#!/usr/bin/env python3
"""
Consolidate Metadata from GEO, SRA, and ENA
Merges results from all three databases and identifies duplicates
"""

import pandas as pd
import numpy as np
from typing import Tuple, List
import logging
import hashlib
from pathlib import Path
from config import DATA_DIR

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('consolidate_metadata.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class MetadataConsolidator:
    """Consolidate metadata from multiple databases"""
    
    def __init__(self):
        self.geo_df = None
        self.sra_df = None
        self.ena_df = None
        self.consolidated_df = None
        self.duplicates = []
    
    def load_database_results(self, geo_file: str = None,
                             sra_file: str = None,
                             ena_file: str = None) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Load results from all three databases"""
        logger.info("Loading database results...")
        
        # Use DATA_DIR paths by default
        if geo_file is None:
            geo_file = str(DATA_DIR / 'geo_datasets.csv')
        if sra_file is None:
            sra_file = str(DATA_DIR / 'sra_experiments.csv')
        if ena_file is None:
            ena_file = str(DATA_DIR / 'ena_runs.csv')
        
        try:
            self.geo_df = pd.read_csv(geo_file)
            logger.info(f"Loaded {len(self.geo_df)} GEO datasets from {geo_file}")
        except FileNotFoundError:
            logger.warning(f"GEO file not found: {geo_file}")
            self.geo_df = pd.DataFrame()
        
        try:
            self.sra_df = pd.read_csv(sra_file)
            logger.info(f"Loaded {len(self.sra_df)} SRA experiments from {sra_file}")
        except FileNotFoundError:
            logger.warning(f"SRA file not found: {sra_file}")
            self.sra_df = pd.DataFrame()
        
        try:
            self.ena_df = pd.read_csv(ena_file)
            logger.info(f"Loaded {len(self.ena_df)} ENA runs from {ena_file}")
        except FileNotFoundError:
            logger.warning(f"ENA file not found: {ena_file}")
            self.ena_df = pd.DataFrame()
        
        return self.geo_df, self.sra_df, self.ena_df
    
    def standardize_geo_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize GEO metadata to common schema"""
        logger.info("Standardizing GEO metadata...")
        
        standardized = pd.DataFrame({
            'database': 'GEO',
            'accession': df.get('gse_id', ''),
            'title': df.get('title', ''),
            'organism': df.get('organism', 'Homo sapiens'),
            'sample_count': df.get('sample_count', 0),
            'platform': df.get('platform', ''),
            'cancer_type': df.get('cancer_type', ''),
            'submission_date': df.get('submission_date', ''),
            'last_update': df.get('last_update', ''),
            'data_type': 'bulk_rnaseq'
        })
        
        return standardized
    
    def standardize_sra_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize SRA metadata to common schema"""
        logger.info("Standardizing SRA metadata...")
        
        standardized = pd.DataFrame({
            'database': 'SRA',
            'accession': df.get('run_accession', ''),
            'study_accession': df.get('study_accession', ''),
            'experiment_accession': df.get('experiment_accession', ''),
            'sample_accession': df.get('sample_accession', ''),
            'title': df.get('title', ''),
            'organism': df.get('organism', 'Homo sapiens'),
            'platform': df.get('platform', ''),
            'instrument_model': df.get('instrument_model', ''),
            'read_length': df.get('read_length', 0),
            'base_count': df.get('base_count', 0),
            'cancer_type': df.get('cancer_type', ''),
            'submission_date': df.get('submission_date', ''),
            'data_type': 'bulk_rnaseq'
        })
        
        return standardized
    
    def standardize_ena_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """Standardize ENA metadata to common schema"""
        logger.info("Standardizing ENA metadata...")
        
        standardized = pd.DataFrame({
            'database': 'ENA',
            'accession': df.get('run_accession', ''),
            'study_accession': df.get('study_accession', ''),
            'experiment_accession': df.get('experiment_accession', ''),
            'sample_accession': df.get('sample_accession', ''),
            'title': df.get('study_title', ''),
            'organism': 'Homo sapiens',
            'platform': 'Illumina',
            'instrument_model': df.get('instrument_model', ''),
            'read_length': df.get('read_length', 0),
            'read_count': df.get('read_count', 0),
            'base_count': df.get('base_count', 0),
            'cancer_type': df.get('cancer_type', ''),
            'submission_date': df.get('first_created', ''),
            'data_type': 'bulk_rnaseq'
        })
        
        return standardized
    
    def detect_duplicates(self, df: pd.DataFrame) -> Tuple[pd.DataFrame, List[List[int]]]:
        """Detect potential duplicate samples across databases"""
        logger.info("Detecting duplicate samples...")
        
        duplicates = []
        
        # Create hash of key identifiers for duplicate detection
        for idx, row in df.iterrows():
            # Look for samples with similar titles and similar read characteristics
            similar_mask = (
                (df['title'].str.lower() == row['title'].lower()) &
                (df['organism'] == row['organism']) &
                (df['cancer_type'] == row['cancer_type']) &
                (df.index != idx)
            )
            
            similar_indices = df[similar_mask].index.tolist()
            
            if similar_indices:
                # Check if this combination hasn't been recorded yet
                duplicate_group = sorted([idx] + similar_indices)
                if duplicate_group not in duplicates:
                    duplicates.append(duplicate_group)
                    logger.warning(f"Potential duplicate detected: {duplicate_group}")
        
        return df, duplicates
    
    def consolidate(self) -> pd.DataFrame:
        """Consolidate all metadata into single dataframe"""
        logger.info("Consolidating metadata from all databases...")
        
        dfs_to_concat = []
        
        if not self.geo_df.empty:
            geo_std = self.standardize_geo_metadata(self.geo_df)
            dfs_to_concat.append(geo_std)
        
        if not self.sra_df.empty:
            sra_std = self.standardize_sra_metadata(self.sra_df)
            dfs_to_concat.append(sra_std)
        
        if not self.ena_df.empty:
            ena_std = self.standardize_ena_metadata(self.ena_df)
            dfs_to_concat.append(ena_std)
        
        if dfs_to_concat:
            self.consolidated_df = pd.concat(dfs_to_concat, ignore_index=True)
            logger.info(f"Consolidated {len(self.consolidated_df)} total samples")
        else:
            logger.error("No data to consolidate")
            self.consolidated_df = pd.DataFrame()
        
        return self.consolidated_df
    
    def add_metadata_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        """Add additional metadata columns for curation"""
        logger.info("Adding metadata columns...")
        
        df['sequencing_type'] = df['read_length'].apply(
            lambda x: 'paired-end' if x > 100 else 'single-end' if x > 0 else 'unknown'
        )
        
        df['quality_flag'] = 'pass'
        
        # Flag low-quality samples
        if 'read_length' in df.columns:
            df.loc[df['read_length'] < 50, 'quality_flag'] = 'low_read_length'
        
        if 'base_count' in df.columns:
            df.loc[df['base_count'] < 1e9, 'quality_flag'] = 'low_depth'
        
        df['curation_status'] = 'pending'
        df['notes'] = ''
        
        return df
    
    def generate_summary_statistics(self, df: pd.DataFrame) -> dict:
        """Generate summary statistics for consolidated metadata"""
        logger.info("Generating summary statistics...")
        
        stats = {
            'total_samples': len(df),
            'by_database': df['database'].value_counts().to_dict(),
            'by_cancer_type': df['cancer_type'].value_counts().to_dict(),
            'by_sequencing_type': df['sequencing_type'].value_counts().to_dict() if 'sequencing_type' in df.columns else {},
            'quality_summary': df['quality_flag'].value_counts().to_dict() if 'quality_flag' in df.columns else {},
            'duplicate_groups': len(self.duplicates)
        }
        
        return stats
    
    def save_consolidated_metadata(self, df: pd.DataFrame, output_file: str = None):
        """Save consolidated metadata to CSV"""
        if output_file is None:
            output_file = str(DATA_DIR / 'consolidated_metadata.csv')
        
        df.to_csv(output_file, index=False)
        logger.info(f"Consolidated metadata saved to {output_file}")
        return output_file
    
    def save_duplicate_report(self, output_file: str = None):
        """Save duplicate detection report"""
        if output_file is None:
            output_file = str(DATA_DIR / 'duplicate_report.txt')
        
        with open(output_file, 'w') as f:
            f.write("Duplicate Sample Detection Report\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total duplicate groups detected: {len(self.duplicates)}\n\n")
            
            for i, group in enumerate(self.duplicates, 1):
                f.write(f"Group {i}: {group}\n")
                if self.consolidated_df is not None:
                    for idx in group:
                        if idx < len(self.consolidated_df):
                            row = self.consolidated_df.iloc[idx]
                            f.write(f"  - {row['database']}: {row['accession']} ({row['title']})\n")
                f.write("\n")
        
        logger.info(f"Duplicate report saved to {output_file}")
        return output_file
    
    def save_summary_report(self, stats: dict, output_file: str = 'consolidation_summary.txt'):
        """Save summary statistics report"""
        with open(output_file, 'w') as f:
            f.write("Metadata Consolidation Summary Report\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"Total Samples: {stats['total_samples']}\n\n")
            
            f.write("Samples by Database:\n")
            for db, count in stats['by_database'].items():
                f.write(f"  {db}: {count}\n")
            
            f.write("\nSamples by Cancer Type:\n")
            for cancer, count in stats['by_cancer_type'].items():
                f.write(f"  {cancer}: {count}\n")
            
            if stats['by_sequencing_type']:
                f.write("\nSamples by Sequencing Type:\n")
                for seq_type, count in stats['by_sequencing_type'].items():
                    f.write(f"  {seq_type}: {count}\n")
            
            if stats['quality_summary']:
                f.write("\nQuality Summary:\n")
                for quality, count in stats['quality_summary'].items():
                    f.write(f"  {quality}: {count}\n")
            
            f.write(f"\nDuplicate Groups Detected: {stats['duplicate_groups']}\n")
        
        logger.info(f"Summary report saved to {output_file}")
        return output_file


def main():
    """Main execution"""
    logger.info("Starting metadata consolidation...")
    
    consolidator = MetadataConsolidator()
    
    # Load results from all databases
    consolidator.load_database_results()
    
    # Consolidate metadata
    consolidated_df = consolidator.consolidate()
    
    # Add metadata columns
    consolidated_df = consolidator.add_metadata_columns(consolidated_df)
    
    # Detect duplicates
    consolidated_df, duplicates = consolidator.detect_duplicates(consolidated_df)
    consolidator.duplicates = duplicates
    
    # Generate statistics
    stats = consolidator.generate_summary_statistics(consolidated_df)
    
    # Save results
    consolidator.save_consolidated_metadata(consolidated_df)
    consolidator.save_duplicate_report()
    consolidator.save_summary_report(stats)
    
    logger.info("Metadata consolidation complete!")
    
    return consolidated_df, stats


if __name__ == "__main__":
    consolidated_df, stats = main()
    
    print("\n" + "=" * 60)
    print("CONSOLIDATION SUMMARY")
    print("=" * 60)
    print(f"\nTotal Samples: {stats['total_samples']}")
    print(f"\nBy Database:")
    for db, count in stats['by_database'].items():
        print(f"  {db}: {count}")
    print(f"\nBy Cancer Type:")
    for cancer, count in stats['by_cancer_type'].items():
        print(f"  {cancer}: {count}")
    print(f"\nDuplicate Groups: {stats['duplicate_groups']}")