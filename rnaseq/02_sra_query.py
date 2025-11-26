#!/usr/bin/env python3
"""
SRA Database Query for Bulk RNA-seq Cancer Datasets
Queries NCBI Sequence Read Archive for RNA-seq experiments across cancer types
"""

import requests
import pandas as pd
import xml.etree.ElementTree as ET
import json
from typing import List, Dict, Optional
import logging
import time
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

try:
    import config
    from utils import (
        retry_with_backoff, setup_logger, validate_accession_format,
        safe_divide, verify_dependencies
    )
except ImportError as e:
    print(f"Error importing configuration modules: {e}")
    sys.exit(1)

logger = setup_logger(__name__, log_file=str(config.LOG_DIR / "sra_query.log"), level=config.LOG_LEVEL)


class SRAQueryEngine:
    """Query SRA database for bulk RNA-seq cancer datasets"""
    
    def __init__(self):
        """Initialize SRA query engine with validated credentials"""
        self.search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self.email = config.NCBI_EMAIL
        self.api_key = config.NCBI_API_KEY
        self.timeout = config.API_TIMEOUT
        self.max_retries = config.API_RETRIES
        self.datasets = []
        
        logger.info(f"SRA Query Engine initialized (timeout={self.timeout}s, retries={self.max_retries})")
    
    def build_sra_query(self, cancer_type: str) -> str:
        """
        Build SRA search query for specific cancer type.
        
        Args:
            cancer_type: Type of cancer to search for
        
        Returns:
            Query string for SRA search
        """
        # Use search terms from config
        query = config.SRA_SEARCH_TERMS.get(cancer_type)
        if not query:
            raise ValueError(
                f"Unknown cancer type: {cancer_type}. "
                f"Supported types: {list(config.SRA_SEARCH_TERMS.keys())}"
            )
        return query
    
    @retry_with_backoff(max_attempts=3, initial_delay=2, exceptions=(requests.RequestException,))
    def search_sra(self, cancer_type: str, max_results: int = 1000) -> List[str]:
        """
        Search SRA for experiments matching cancer type.
        
        Args:
            cancer_type: Type of cancer to search for
            max_results: Maximum number of results to return
        
        Returns:
            List of SRA experiment IDs
        """
        logger.info(f"Searching SRA for {cancer_type} experiments (max_results={max_results})...")
        
        query = self.build_sra_query(cancer_type)
        params = {
            'db': 'sra',
            'term': query,
            'retmax': max_results,
            'rettype': 'json',
            'retmode': 'json',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        response = requests.get(self.search_url, params=params, timeout=self.timeout)
        response.raise_for_status()
        
        # Validate response has content
        if not response.text or not response.text.strip():
            logger.warning(f"Empty response from SRA search for {cancer_type}")
            raise requests.RequestException(f"Empty response from API")
        
        try:
            data = response.json()
        except (ValueError, json.JSONDecodeError) as e:
            logger.error(f"Error parsing JSON response: {e}")
            logger.error(f"Response content: {response.text[:200]}")
            raise requests.RequestException(f"Invalid JSON response: {e}")
        
        if not data or 'esearchresult' not in data:
            logger.warning(f"Unexpected response format for {cancer_type}: {data}")
            return []
        
        if 'idlist' not in data['esearchresult']:
            logger.debug(f"No results found for {cancer_type}")
            return []
        
        sra_ids = data['esearchresult']['idlist']
        logger.info(f"Found {len(sra_ids)} experiments for {cancer_type}")
        return sra_ids
    
    @retry_with_backoff(max_attempts=3, initial_delay=2, exceptions=(requests.RequestException, ET.ParseError))
    def fetch_experiment_metadata(self, sra_id: str) -> Optional[Dict]:
        """
        Fetch detailed metadata for SRA experiment.
        
        Args:
            sra_id: SRA experiment ID
        
        Returns:
            Metadata dictionary or None if fetch fails
        """
        logger.debug(f"Fetching metadata for SRA ID: {sra_id}...")
        
        params = {
            'db': 'sra',
            'id': sra_id,
            'rettype': 'xml',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        response = requests.get(self.fetch_url, params=params, timeout=self.timeout)
        response.raise_for_status()
        
        if not response.content:
            logger.warning(f"Empty response for {sra_id}")
            return None
        
        # Parse XML response
        try:
            root = ET.fromstring(response.content)
        except ET.ParseError as e:
            logger.error(f"Error parsing XML for {sra_id}: {e}")
            return None
        
        metadata = {
            'sra_id': sra_id,
            'experiment_accession': '',
            'run_accession': '',
            'sample_accession': '',
            'study_accession': '',
            'title': '',
            'organism': 'Homo sapiens',
            'library_strategy': '',
            'library_source': '',
            'library_selection': '',
            'platform': '',
            'instrument_model': '',
            'read_length': 0,
            'base_count': 0,
            'run_count': 0,
            'submission_date': '',
            'publication_date': ''
        }
        
        try:
            # Extract experiment information
            for experiment in root.findall('.//EXPERIMENT'):
                metadata['experiment_accession'] = experiment.get('accession', '')
                title_elem = experiment.find('.//TITLE')
                if title_elem is not None and title_elem.text:
                    metadata['title'] = title_elem.text
                
                # Extract library info
                design = experiment.find('.//LIBRARY_DESCRIPTOR')
                if design is not None:
                    strategy = design.find('LIBRARY_STRATEGY')
                    if strategy is not None and strategy.text:
                        metadata['library_strategy'] = strategy.text
                    source = design.find('LIBRARY_SOURCE')
                    if source is not None and source.text:
                        metadata['library_source'] = source.text
                    selection = design.find('LIBRARY_SELECTION')
                    if selection is not None and selection.text:
                        metadata['library_selection'] = selection.text
                
                # Extract platform info
                platform = experiment.find('.//PLATFORM')
                if platform is not None:
                    for child in platform:
                        metadata['platform'] = child.tag
                        instrument = child.find('INSTRUMENT_MODEL')
                        if instrument is not None and instrument.text:
                            metadata['instrument_model'] = instrument.text
            
            # Extract run information
            for run in root.findall('.//RUN'):
                metadata['run_accession'] = run.get('accession', '')
                metadata['run_count'] += 1
                
                # Extract read length and base count (with safe division)
                total_bases_str = run.get('total_bases', '0')
                total_spots_str = run.get('total_spots', '0')
                
                try:
                    total_bases = int(total_bases_str)
                    total_spots = int(total_spots_str)
                    
                    metadata['base_count'] = total_bases
                    
                    if total_spots > 0:
                        metadata['read_length'] = total_bases // total_spots
                    else:
                        logger.warning(f"Zero total_spots for {sra_id}")
                        
                except ValueError as e:
                    logger.warning(f"Could not parse read metrics for {sra_id}: {e}")
            
            # Extract sample and study info
            for sample in root.findall('.//SAMPLE'):
                metadata['sample_accession'] = sample.get('accession', '')
            
            for study in root.findall('.//STUDY'):
                metadata['study_accession'] = study.get('accession', '')
            
            return metadata
        
        except Exception as e:
            logger.error(f"Error extracting metadata for {sra_id}: {e}")
            return None
    
    def query_all_cancers(self, max_results_per_cancer: int = 50) -> pd.DataFrame:
        """
        Query all cancer types and compile results.
        
        Args:
            max_results_per_cancer: Maximum results to retrieve per cancer type
        
        Returns:
            DataFrame with all retrieved experiments
        """
        cancer_types = list(config.CANCER_TYPES.keys())
        all_experiments = []
        
        for cancer_type in cancer_types:
            try:
                sra_ids = self.search_sra(cancer_type, max_results=max_results_per_cancer)
                
                for sra_id in sra_ids:
                    try:
                        metadata = self.fetch_experiment_metadata(sra_id)
                        if metadata:
                            metadata['cancer_type'] = cancer_type
                            all_experiments.append(metadata)
                    except Exception as e:
                        logger.warning(f"Failed to fetch metadata for {sra_id}: {e}")
                        continue
                    
                    # Rate limiting
                    time.sleep(0.4)
            except Exception as e:
                logger.error(f"Failed to query cancer type {cancer_type}: {e}")
                continue
        
        if not all_experiments:
            logger.warning("No experiments retrieved from SRA")
            return pd.DataFrame()
        
        df = pd.DataFrame(all_experiments)
        logger.info(f"Total experiments retrieved: {len(df)}")
        return df
    
    def filter_bulk_rnaseq(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter for bulk RNA-seq (exclude single-cell and other non-bulk).
        
        Args:
            df: Input DataFrame with experiment metadata
        
        Returns:
            Filtered DataFrame containing only bulk RNA-seq experiments
        """
        if df.empty:
            logger.warning("Input DataFrame is empty for bulk RNA-seq filtering")
            return df
        
        initial_count = len(df)
        
        # Filter for RNA-seq strategy
        mask = df['library_strategy'].str.contains('RNA-Seq', case=False, na=False)
        
        # Exclude single-cell
        exclude_keywords = ['single cell', 'scRNA', '10x', 'droplet']
        mask &= ~df['title'].str.lower().str.contains('|'.join(exclude_keywords), na=False)
        
        # Exclude low-quality runs (with safe value checks)
        mask &= df['read_length'].fillna(0) >= config.MIN_READ_LENGTH
        mask &= df['base_count'].fillna(0) >= int(1e9)  # At least 1 billion bases
        
        filtered_df = df[mask].copy()
        
        n_removed = initial_count - len(filtered_df)
        logger.info(f"Removed {n_removed} experiments during filtering. Remaining: {len(filtered_df)}")
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
            output_file = str(config.DATA_DIR / "sra_experiments.csv")
        
        try:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(output_file, index=False)
            logger.info(f"Results saved to {output_file} ({len(df)} experiments)")
            return output_file
        except IOError as e:
            logger.error(f"Failed to save results to {output_file}: {e}")
            raise


def main():
    """Main execution"""
    try:
        logger.info("=" * 60)
        logger.info("Starting SRA database query...")
        logger.info("=" * 60)
        
        engine = SRAQueryEngine()
        
        # Query all cancer types
        results_df = engine.query_all_cancers(max_results_per_cancer=config.MAX_SRA_RECORDS)
        
        if results_df.empty:
            logger.error("No experiments retrieved from SRA query")
            return None
        
        # Filter for bulk RNA-seq
        bulk_rnaseq_df = engine.filter_bulk_rnaseq(results_df)
        
        if bulk_rnaseq_df.empty:
            logger.error("All experiments filtered out - no bulk RNA-seq data found")
            return None
        
        # Save results
        output_file = engine.save_results(bulk_rnaseq_df)
        
        logger.info("=" * 60)
        logger.info(f"SRA query complete. Found {len(bulk_rnaseq_df)} bulk RNA-seq experiments")
        logger.info("=" * 60)
        
        return bulk_rnaseq_df
    
    except Exception as e:
        logger.error(f"Fatal error during SRA query: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    results = main()
    
    if results is not None and not results.empty:
        print(f"\nSummary:")
        print(f"Total experiments: {len(results)}")
        print(f"\nExperiments by cancer type:")
        print(results['cancer_type'].value_counts())
        print(f"\nPlatform distribution:")
        print(results['platform'].value_counts())
    else:
        print("\nNo results to display")
        sys.exit(1)