#!/usr/bin/env python3
"""
ENA Database Query for Bulk RNA-seq Cancer Datasets
Queries European Nucleotide Archive for RNA-seq studies across cancer types
"""

import requests
import pandas as pd
from typing import List, Dict
import logging
import time
import json
from config import MAX_ENA_RECORDS, DATA_DIR

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('ena_query.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ENAQueryEngine:
    """Query ENA database for bulk RNA-seq cancer datasets"""
    
    def __init__(self):
        self.portal_url = "https://www.ebi.ac.uk/ena/portal/api"
        self.search_url = f"{self.portal_url}/search"
        self.retrieve_url = f"{self.portal_url}/retrieve"
        self.datasets = []
    
    def build_ena_query(self, cancer_type: str) -> str:
        """Build ENA search query for all Homo sapiens RNA-Seq studies using taxid"""
        # ENA portal API uses simpler syntax - just taxid query
        # Studies will be filtered by RNA-seq in the run fetching step
        return 'tax_eq(9606)'
    
    def search_ena(self, cancer_type: str, max_results: int = 100) -> List[Dict]:
        """Search ENA for studies matching cancer type"""
        logger.info(f"Searching ENA for {cancer_type} studies...")
        
        query = self.build_ena_query(cancer_type)
        
        params = {
            'query': query,
            'result': 'study',
            'format': 'json',
            'limit': max_results,
            'offset': 0
        }
        
        try:
            response = requests.get(self.search_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            # ENA API returns list directly, not dict with 'results' key
            if isinstance(data, list):
                logger.info(f"Found {len(data)} studies for {cancer_type}")
                return data
            elif isinstance(data, dict) and 'results' in data:
                studies = data['results']
                logger.info(f"Found {len(studies)} studies for {cancer_type}")
                return studies
            else:
                logger.warning(f"Unexpected response format for {cancer_type}")
                return []
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Error searching ENA: {e}")
            return []
    
    def fetch_study_runs(self, study_accession: str) -> List[Dict]:
        """Fetch all runs associated with a study"""
        logger.info(f"Fetching runs for study {study_accession}...")
        
        params = {
            'query': f'study_accession="{study_accession}"',
            'result': 'read_run',
            'format': 'json',
            'limit': 1000,
            'fields': 'run_accession,sample_accession,experiment_accession,library_strategy,library_source,library_selection,instrument_model,read_count,base_count,first_created,last_updated'
        }
        
        try:
            response = requests.get(self.search_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            # Handle both list and dict responses
            if isinstance(data, list):
                runs = data
            elif isinstance(data, dict) and 'results' in data:
                runs = data['results']
            else:
                runs = []
            
            logger.info(f"Found {len(runs)} runs for study {study_accession}")
            return runs
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching runs for {study_accession}: {e}")
            return []
    
    def fetch_sample_metadata(self, sample_accession: str) -> Dict:
        """Fetch metadata for a specific sample"""
        logger.info(f"Fetching metadata for sample {sample_accession}...")
        
        params = {
            'query': f'sample_accession="{sample_accession}"',
            'result': 'sample',
            'format': 'json',
            'limit': 1
        }
        
        try:
            response = requests.get(self.search_url, params=params, timeout=30)
            response.raise_for_status()
            data = response.json()
            
            # Handle both list and dict responses
            if isinstance(data, list) and len(data) > 0:
                return data[0]
            elif isinstance(data, dict) and 'results' in data and len(data['results']) > 0:
                return data['results'][0]
            else:
                return {}
                
        except requests.exceptions.RequestException as e:
            logger.error(f"Error fetching sample metadata for {sample_accession}: {e}")
            return {}
    
    def query_all_cancers(self, max_results_per_cancer: int = None) -> pd.DataFrame:
        """Query all cancer types and compile results"""
        if max_results_per_cancer is None:
            max_results_per_cancer = MAX_ENA_RECORDS
        
        all_runs = []
        
        # Only query once for all RNA-Seq studies in Homo sapiens
        studies = self.search_ena('all', max_results=max_results_per_cancer)
        logger.info(f"Processing {len(studies)} studies...")
        
        for idx, study in enumerate(studies):
            study_accession = study.get('study_accession', '')
            # ENA returns 'description' field, not 'study_title' or 'study_description'
            study_desc = study.get('description', '').lower()
            
            if not study_accession:
                continue
            
            logger.debug(f"Processing study {idx+1}/{len(studies)}: {study_accession}")
            
            # Fetch runs for this study
            runs = self.fetch_study_runs(study_accession)
            
            for run in runs:
                run_data = {
                    'study_accession': study_accession,
                    'run_accession': run.get('run_accession', ''),
                    'description': study_desc,
                    'raw_data': json.dumps(run)
                }
                
                # Try to extract useful fields if they exist
                for key in ['experiment_accession', 'sample_accession', 'instrument_model', 
                            'library_strategy', 'library_source', 'library_selection',
                            'read_count', 'base_count', 'first_created', 'last_updated']:
                    if key in run:
                        run_data[key] = run[key]
                
                all_runs.append(run_data)
            
            # Rate limiting
            time.sleep(0.3)
        
        df = pd.DataFrame(all_runs)
        logger.info(f"Total runs retrieved: {len(df)}")
        return df
    
    def filter_bulk_rnaseq(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filter for bulk RNA-seq (exclude single-cell and low-quality)"""
        # Handle empty DataFrame
        if df.empty:
            logger.warning("DataFrame is empty - no runs to filter")
            return df
        
        # Start with all rows
        mask = pd.Series([True] * len(df), index=df.index)
        
        # Filter for RNA-Seq strategy if column exists
        if 'library_strategy' in df.columns:
            mask &= df['library_strategy'].fillna('').astype(str).str.contains('RNA-Seq', case=False, na=False)
        
        # Exclude single-cell
        if 'description' in df.columns:
            exclude_keywords = ['single cell', 'scRNA', '10x', 'droplet', 'microfluidic']
            mask &= ~df['description'].fillna('').astype(str).str.lower().str.contains('|'.join(exclude_keywords), na=False)
        
        # Exclude low-quality runs if columns exist
        # Convert read_count and base_count to numeric if they exist
        if 'read_count' in df.columns:
            try:
                df['read_count'] = pd.to_numeric(df['read_count'], errors='coerce')
                mask &= df['read_count'] >= 20e6  # At least 20 million reads
            except Exception as e:
                logger.warning(f"Could not filter by read_count: {e}")
        
        if 'base_count' in df.columns:
            try:
                df['base_count'] = pd.to_numeric(df['base_count'], errors='coerce')
            except Exception as e:
                logger.warning(f"Could not convert base_count: {e}")
        
        # Calculate read length if needed
        if 'read_count' in df.columns and 'base_count' in df.columns:
            try:
                df['read_length'] = (df['base_count'] / df['read_count']).fillna(0).astype(int)
                mask &= df['read_length'] >= 50
            except Exception as e:
                logger.warning(f"Could not calculate read_length: {e}")
        
        filtered_df = df[mask].copy()
        logger.info(f"After bulk RNA-seq filtering: {len(filtered_df)} runs")
        return filtered_df
    
    def save_results(self, df: pd.DataFrame, output_file: str = None):
        """Save query results to CSV"""
        if output_file is None:
            output_file = str(DATA_DIR / 'ena_runs.csv')
        
        df.to_csv(output_file, index=False)
        logger.info(f"Results saved to {output_file}")
        return output_file


def main():
    """Main execution"""
    logger.info("Starting ENA database query...")
    
    try:
        engine = ENAQueryEngine()
        
        # Query all cancer types
        results_df = engine.query_all_cancers()  # Uses MAX_ENA_RECORDS from config
        
        if results_df.empty:
            logger.warning("No results retrieved from ENA")
            return pd.DataFrame()
        
        # Filter for bulk RNA-seq
        bulk_rnaseq_df = engine.filter_bulk_rnaseq(results_df)
        
        if bulk_rnaseq_df.empty:
            logger.warning("No bulk RNA-seq runs after filtering")
            return pd.DataFrame()
        
        # Save results
        output_file = engine.save_results(bulk_rnaseq_df)
        
        logger.info(f"ENA query complete. Found {len(bulk_rnaseq_df)} bulk RNA-seq runs")
        
        return bulk_rnaseq_df
    
    except Exception as e:
        logger.error(f"Fatal error in ENA query: {e}", exc_info=True)
        return pd.DataFrame()


if __name__ == "__main__":
    results = main()
    if not results.empty:
        print(f"\nSummary:")
        print(f"Total runs: {len(results)}")
        if 'cancer_type' in results.columns:
            print(f"\nRuns by cancer type:")
            print(results['cancer_type'].value_counts())
        if 'instrument_model' in results.columns:
            print(f"\nInstrument model distribution:")
            print(results['instrument_model'].value_counts())
        if 'read_length' in results.columns:
            print(f"\nRead length statistics:")
            print(results['read_length'].describe())
    else:
        print("\nNo results to display")
        logger.warning("Pipeline exiting with no results from ENA query")