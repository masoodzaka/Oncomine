#!/usr/bin/env python3
"""
GEO Database Query for Bulk RNA-seq Cancer Datasets
Queries Gene Expression Omnibus for bulk RNA-seq studies across multiple cancer types
"""

import requests
import pandas as pd
import json
import xml.etree.ElementTree as ET
from datetime import datetime
import logging
from typing import List, Dict, Tuple, Optional
import time
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

try:
    import config
    from utils import (
        retry_with_backoff, setup_logger, validate_csv_structure,
        consolidate_dataframes, safe_divide
    )
except ImportError as e:
    print(f"Error importing configuration modules: {e}")
    print("Please ensure config.py and utils.py are in the same directory")
    sys.exit(1)

# Configure logging
logger = setup_logger(__name__, log_file=str(config.LOG_DIR / "geo_query.log"), level=config.LOG_LEVEL)


class GEOQueryEngine:
    """Query GEO database for bulk RNA-seq cancer datasets"""
    
    def __init__(self):
        """Initialize GEO query engine with validated credentials"""
        self.base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
        self.search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        
        # Use credentials from config (already validated)
        self.email = config.NCBI_EMAIL
        self.api_key = config.NCBI_API_KEY
        self.timeout = config.API_TIMEOUT
        self.max_retries = config.API_RETRIES
        self.datasets = []
        
        logger.info(f"GEO Query Engine initialized (timeout={self.timeout}s, retries={self.max_retries})")
        
    def build_search_query(self, cancer_type: str) -> str:
        """Build GEO search query for specific cancer type"""
        # Use search terms from config
        query = config.GEO_SEARCH_TERMS.get(cancer_type)
        if not query:
            raise ValueError(
                f"Unknown cancer type: {cancer_type}. "
                f"Supported types: {list(config.GEO_SEARCH_TERMS.keys())}"
            )
        return query
    
    @retry_with_backoff(max_attempts=3, initial_delay=2, exceptions=(requests.RequestException,))
    def search_geo(self, cancer_type: str, max_results: int = 100) -> List[str]:
        """
        Search GEO for datasets matching cancer type.
        
        Args:
            cancer_type: Type of cancer to search for
            max_results: Maximum number of results to return
        
        Returns:
            List of GEO dataset IDs
        """
        logger.info(f"Searching GEO for {cancer_type} datasets (max_results={max_results})...")
        
        query = self.build_search_query(cancer_type)
        params = {
            'db': 'gds',
            'term': query,
            'retmax': max_results,
            'rettype': 'json',  # NCBI returns XML regardless, so we parse XML
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        response = requests.get(self.search_url, params=params, timeout=self.timeout)
        response.raise_for_status()
        
        # Check HTTP status and rate limiting
        if response.status_code == 429:
            logger.warning("Rate limited by NCBI - backing off")
            raise requests.RequestException("Rate limited (429)")
        
        # Validate response has content
        if not response.text or not response.text.strip():
            logger.warning(f"Empty response from GEO search for {cancer_type}")
            raise requests.RequestException(f"Empty response from API")
        
        try:
            # NCBI returns XML, not JSON, despite rettype=json parameter
            root = ET.fromstring(response.text)
            
            # Parse XML to extract ID list
            id_list = root.find('.//IdList')
            if id_list is None:
                logger.info(f"No results found for {cancer_type}")
                return []
            
            gse_ids = [id_elem.text for id_elem in id_list.findall('Id')]
            logger.info(f"Found {len(gse_ids)} datasets for {cancer_type}")
            return gse_ids
            
        except ET.ParseError as e:
            logger.error(f"Failed to parse XML response: {e}")
            logger.debug(f"Response content: {response.text[:500]}")
            raise requests.RequestException(f"Invalid XML response: {e}")
        except (ValueError, AttributeError, TypeError) as e:
            logger.error(f"Error extracting IDs from response: {e}")
            logger.debug(f"Response content: {response.text[:500]}")
            raise requests.RequestException(f"Error parsing XML: {e}")
    
    @retry_with_backoff(max_attempts=3, initial_delay=2, exceptions=(requests.RequestException,))
    def fetch_dataset_metadata(self, gse_id: str) -> Optional[Dict]:
        """
        Fetch detailed metadata for a GEO dataset using NCBI esummary API.
        The esummary API returns structured XML data with dataset metadata.
        
        Args:
            gse_id: GEO dataset ID (numerical ID from search results)
        
        Returns:
            Metadata dictionary or None if fetch fails
        """
        logger.debug(f"Fetching metadata for {gse_id}...")
        
        # Use esummary API which returns structured XML with metadata
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            'db': 'gds',
            'id': gse_id,
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        response = requests.get(url, params=params, timeout=self.timeout)
        response.raise_for_status()
        
        # Check rate limiting (HTTP 429)
        if response.status_code == 429:
            logger.warning(f"Rate limited by NCBI for {gse_id} - backing off")
            raise requests.RequestException("Rate limited (429)")
        
        # Validate response has content
        if not response.text or not response.text.strip():
            logger.warning(f"Empty response from NCBI esummary for {gse_id}")
            return None
        
        try:
            # Parse XML response from esummary
            root = ET.fromstring(response.text)
            
            # Initialize metadata with defaults
            metadata = {
                'gse_id': gse_id,
                'title': 'Unknown',
                'summary': '',
                'organism': 'Homo sapiens',
                'sample_count': 0,
                'platform': '',
                'submission_date': '',
                'last_update': '',
            }
            
            # Find the DocSum element which contains the dataset information
            docsum = root.find('.//DocSum')
            if docsum is None:
                logger.warning(f"No DocSum found in esummary response for {gse_id}")
                return None
            
            # Extract data from Item elements
            for item in docsum.findall('Item'):
                name = item.get('Name', '')
                text = item.text or ''
                
                if name == 'Accession':
                    metadata['gse_id'] = text or gse_id
                elif name == 'title':
                    metadata['title'] = text or 'Unknown'
                elif name == 'summary':
                    metadata['summary'] = text or ''
                elif name == 'taxon':
                    metadata['organism'] = text or 'Homo sapiens'
                elif name == 'n_samples':
                    try:
                        metadata['sample_count'] = int(text) if text else 0
                    except (ValueError, TypeError):
                        metadata['sample_count'] = 0
                elif name == 'GPL' or name == 'PlatformTitle':
                    metadata['platform'] = text or ''
                elif name == 'PDAT':  # Publication date
                    metadata['submission_date'] = text or ''
                elif name == 'updatedate':
                    metadata['last_update'] = text or ''
            
            logger.debug(f"Successfully extracted metadata for {gse_id}: {metadata['title'][:80]}")
            return metadata
            
        except ET.ParseError as e:
            logger.error(f"Failed to parse XML for {gse_id}: {e}")
            logger.debug(f"Response content: {response.text[:500]}")
            raise requests.RequestException(f"Invalid XML response: {e}")
        except (KeyError, TypeError, AttributeError) as e:
            logger.error(f"Error extracting metadata for {gse_id}: {e}")
            return None
    
    def query_all_cancers(self, max_results_per_cancer: int = 50) -> pd.DataFrame:
        """
        Query all cancer types and compile results.
        
        Args:
            max_results_per_cancer: Maximum results to retrieve per cancer type
        
        Returns:
            DataFrame with all retrieved datasets
        """
        cancer_types = list(config.CANCER_TYPES.keys())  # Use cancer types from config
        all_datasets = []
        failed_queries = []
        
        for cancer_type in cancer_types:
            try:
                gse_ids = self.search_geo(cancer_type, max_results=max_results_per_cancer)
                logger.info(f"Processing {len(gse_ids)} results for {cancer_type}")
                
                for gse_id in gse_ids:
                    try:
                        metadata = self.fetch_dataset_metadata(gse_id)
                        if metadata is not None:
                            metadata['cancer_type'] = cancer_type
                            all_datasets.append(metadata)
                        else:
                            logger.debug(f"Skipped {gse_id} - no metadata returned")
                    except requests.RequestException as e:
                        logger.warning(f"Failed to fetch metadata for {gse_id}: {e}")
                        continue
                    except Exception as e:
                        logger.error(f"Unexpected error fetching {gse_id}: {e}")
                        continue
                    
                    # Rate limiting - respectful delay between requests
                    time.sleep(0.4)
                    
            except requests.RequestException as e:
                logger.error(f"Failed to query cancer type {cancer_type}: {e}")
                failed_queries.append((cancer_type, str(e)))
                continue
            except Exception as e:
                logger.error(f"Unexpected error querying {cancer_type}: {e}")
                failed_queries.append((cancer_type, str(e)))
                continue
        
        if failed_queries:
            logger.warning(f"Failed to query {len(failed_queries)} cancer types: {failed_queries}")
        
        if not all_datasets:
            logger.warning("No datasets retrieved from GEO")
            return pd.DataFrame()
        
        df = pd.DataFrame(all_datasets)
        logger.info(f"Total datasets retrieved: {len(df)}")
        return df
    
    def filter_bulk_rnaseq(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter for bulk RNA-seq datasets (exclude single-cell).
        
        Args:
            df: Input DataFrame with dataset metadata
        
        Returns:
            Filtered DataFrame containing only bulk RNA-seq datasets
        """
        if df.empty:
            logger.warning("Input DataFrame is empty for bulk RNA-seq filtering")
            return df
        
        # Filter out single-cell studies using keywords from config
        mask = ~df['title'].str.lower().str.contains('|'.join(config.BULK_RNASEQ_EXCLUDE_KEYWORDS), na=False)
        filtered_df = df[mask].copy()
        
        n_removed = len(df) - len(filtered_df)
        logger.info(f"Removed {n_removed} single-cell studies. Remaining: {len(filtered_df)}")
        return filtered_df
    
    def save_results(self, df: pd.DataFrame, output_file: Optional[str] = None) -> str:
        """
        Save query results to CSV.
        
        Args:
            df: DataFrame to save
            output_file: Output file path (uses config default if not specified)
        
        Returns:
            Path to saved file
        """
        if df.empty:
            logger.warning("Attempting to save empty DataFrame - no results to write")
            return ""
        
        if output_file is None:
            output_file = str(config.DATA_DIR / "geo_datasets.csv")
        
        try:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_file, index=False)
            logger.info(f"Results saved to {output_file} ({len(df)} datasets)")
            return output_file
        except IOError as e:
            logger.error(f"Failed to save results to {output_file}: {e}")
            raise


def main():
    """Main execution"""
    try:
        logger.info("=" * 60)
        logger.info("Starting GEO database query...")
        logger.info("=" * 60)
        
        engine = GEOQueryEngine()
        
        # Query all cancer types
        results_df = engine.query_all_cancers(max_results_per_cancer=config.MAX_GEO_RECORDS)
        
        if results_df.empty:
            logger.error("No datasets retrieved from GEO query")
            return None
        
        # Filter for bulk RNA-seq
        bulk_rnaseq_df = engine.filter_bulk_rnaseq(results_df)
        
        if bulk_rnaseq_df.empty:
            logger.error("All datasets filtered out - no bulk RNA-seq data found")
            return None
        
        # Save results
        output_file = engine.save_results(bulk_rnaseq_df)
        
        logger.info("=" * 60)
        logger.info(f"GEO query complete. Found {len(bulk_rnaseq_df)} bulk RNA-seq datasets")
        logger.info("=" * 60)
        
        return bulk_rnaseq_df
    
    except Exception as e:
        logger.error(f"Fatal error during GEO query: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    results = main()
    
    if results is not None and not results.empty:
        print(f"\nSummary:")
        print(f"Total datasets: {len(results)}")
        print(f"\nDatasets by cancer type:")
        print(results['cancer_type'].value_counts())
    else:
        print("\nNo results to display")
        sys.exit(1)