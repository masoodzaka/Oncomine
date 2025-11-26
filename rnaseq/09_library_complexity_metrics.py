#!/usr/bin/env python3
"""
Library Complexity Metrics Pipeline
Calculates duplication rates, complexity metrics, and library quality indicators
"""

import subprocess
import pandas as pd
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import json

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('library_complexity.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class LibraryComplexityMetrics:
    """Calculate library complexity and duplication metrics"""
    
    def __init__(self, fastq_dir: str = 'fastq_downloads',
                 output_dir: str = 'complexity_reports',
                 max_workers: int = 4):
        self.fastq_dir = Path(fastq_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.results = []
    
    def calculate_sequence_complexity(self, fastq_file: Path) -> Dict:
        """Calculate sequence complexity metrics"""
        logger.info(f"Calculating complexity for {fastq_file.name}")
        
        result_dict = {
            'file': fastq_file.name,
            'total_reads': 0,
            'unique_sequences': 0,
            'duplicate_reads': 0,
            'duplication_rate': 0.0,
            'complexity_score': 0.0,
            'entropy': 0.0,
            'gc_content': 0.0,
            'quality_assessment': 'unknown'
        }
        
        try:
            import gzip
            from collections import Counter
            import math
            
            sequence_counts = Counter()
            total_bases = 0
            gc_count = 0
            quality_scores = []
            
            with gzip.open(fastq_file, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    
                    # FASTQ format: sequence on line 2, quality on line 4
                    if line_num % 4 == 2:
                        sequence = line.strip().upper()
                        sequence_counts[sequence] += 1
                        
                        # Count GC
                        gc_count += sequence.count('G') + sequence.count('C')
                        total_bases += len(sequence)
                        
                    elif line_num % 4 == 0:
                        quality = line.strip()
                        quality_scores.extend([ord(q) - 33 for q in quality])
            
            # Calculate metrics
            total_reads = sum(sequence_counts.values())
            unique_sequences = len(sequence_counts)
            duplicate_reads = total_reads - unique_sequences
            
            result_dict['total_reads'] = total_reads
            result_dict['unique_sequences'] = unique_sequences
            result_dict['duplicate_reads'] = duplicate_reads
            
            # Duplication rate
            if total_reads > 0:
                result_dict['duplication_rate'] = (duplicate_reads / total_reads) * 100
            
            # Complexity score (unique sequences / total reads)
            if total_reads > 0:
                result_dict['complexity_score'] = (unique_sequences / total_reads) * 100
            
            # Shannon entropy (measure of sequence diversity)
            if total_reads > 0:
                entropy = 0.0
                for count in sequence_counts.values():
                    p = count / total_reads
                    if p > 0:
                        entropy -= p * math.log2(p)
                result_dict['entropy'] = entropy
            
            # GC content
            if total_bases > 0:
                result_dict['gc_content'] = (gc_count / total_bases) * 100
            
            # Quality assessment
            if quality_scores:
                mean_quality = sum(quality_scores) / len(quality_scores)
                if mean_quality >= 30:
                    result_dict['quality_assessment'] = 'high'
                elif mean_quality >= 20:
                    result_dict['quality_assessment'] = 'medium'
                else:
                    result_dict['quality_assessment'] = 'low'
            
            logger.info(f"Complexity calculation complete: {result_dict['complexity_score']:.1f}% unique")
            
        except Exception as e:
            logger.error(f"Error calculating complexity: {e}")
        
        return result_dict
    
    def estimate_library_size(self, fastq_file: Path) -> Dict:
        """Estimate library size using Chao1 estimator"""
        logger.info(f"Estimating library size for {fastq_file.name}")
        
        result_dict = {
            'file': fastq_file.name,
            'observed_species': 0,
            'singleton_species': 0,
            'doubleton_species': 0,
            'chao1_estimate': 0.0,
            'coverage_estimate': 0.0
        }
        
        try:
            import gzip
            from collections import Counter
            
            sequence_counts = Counter()
            
            with gzip.open(fastq_file, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    
                    if line_num % 4 == 2:
                        sequence = line.strip().upper()
                        sequence_counts[sequence] += 1
            
            # Count species by frequency
            frequency_counts = Counter(sequence_counts.values())
            
            observed_species = len(sequence_counts)
            singleton_species = frequency_counts.get(1, 0)
            doubleton_species = frequency_counts.get(2, 0)
            
            result_dict['observed_species'] = observed_species
            result_dict['singleton_species'] = singleton_species
            result_dict['doubleton_species'] = doubleton_species
            
            # Chao1 estimator for species richness
            if doubleton_species > 0:
                chao1 = observed_species + (singleton_species ** 2) / (2 * doubleton_species)
            else:
                chao1 = observed_species + (singleton_species * (singleton_species - 1)) / 2
            
            result_dict['chao1_estimate'] = chao1
            
            # Coverage estimate (Good's coverage)
            total_reads = sum(sequence_counts.values())
            if total_reads > 0:
                coverage = 1 - (singleton_species / total_reads)
                result_dict['coverage_estimate'] = coverage * 100
            
            logger.info(f"Library size estimation complete: Chao1={chao1:.0f}")
            
        except Exception as e:
            logger.error(f"Error estimating library size: {e}")
        
        return result_dict
    
    def find_fastq_files(self) -> List[Path]:
        """Find all FASTQ files"""
        logger.info(f"Searching for FASTQ files in {self.fastq_dir}")
        fastq_files = list(self.fastq_dir.rglob('*.fastq.gz'))
        logger.info(f"Found {len(fastq_files)} FASTQ files")
        return fastq_files
    
    def process_sample(self, fastq_file: Path) -> Dict:
        """Process a single sample"""
        result = {
            'file': fastq_file.name,
            'path': str(fastq_file),
            'complexity_metrics': {},
            'library_size_estimate': {}
        }
        
        # Calculate complexity
        complexity = self.calculate_sequence_complexity(fastq_file)
        result['complexity_metrics'] = complexity
        
        # Estimate library size
        library_size = self.estimate_library_size(fastq_file)
        result['library_size_estimate'] = library_size
        
        self.results.append(result)
        return result
    
    def run_complexity_analysis(self):
        """Execute library complexity analysis pipeline"""
        logger.info("Starting library complexity analysis pipeline...")
        
        fastq_files = self.find_fastq_files()
        
        if not fastq_files:
            logger.error("No FASTQ files found!")
            return
        
        logger.info(f"Processing {len(fastq_files)} files with {self.max_workers} workers")
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self.process_sample, fastq_file): fastq_file
                for fastq_file in fastq_files
            }
            
            completed = 0
            for future in as_completed(futures):
                try:
                    result = future.result()
                    completed += 1
                    logger.info(f"[{completed}/{len(fastq_files)}] ✓ {result['file']}")
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
        
        logger.info("Library complexity analysis complete!")
    
    def generate_complexity_report(self) -> pd.DataFrame:
        """Generate library complexity report"""
        logger.info("Generating library complexity report...")
        
        report_data = []
        
        for result in self.results:
            complexity = result['complexity_metrics']
            library = result['library_size_estimate']
            
            report_row = {
                'file': result['file'],
                'total_reads': complexity.get('total_reads', 0),
                'unique_sequences': complexity.get('unique_sequences', 0),
                'duplication_rate': complexity.get('duplication_rate', 0.0),
                'complexity_score': complexity.get('complexity_score', 0.0),
                'entropy': complexity.get('entropy', 0.0),
                'gc_content': complexity.get('gc_content', 0.0),
                'quality_assessment': complexity.get('quality_assessment', 'unknown'),
                'chao1_estimate': library.get('chao1_estimate', 0.0),
                'coverage_estimate': library.get('coverage_estimate', 0.0)
            }
            report_data.append(report_row)
        
        report_df = pd.DataFrame(report_data)
        report_df.to_csv(self.output_dir / 'library_complexity_report.csv', index=False)
        logger.info("Library complexity report saved")
        
        return report_df
    
    def generate_summary(self, report_df: pd.DataFrame):
        """Generate summary statistics"""
        logger.info("Generating library complexity summary...")
        
        summary_file = self.output_dir / 'LIBRARY_COMPLEXITY_SUMMARY.txt'
        
        with open(summary_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("LIBRARY COMPLEXITY METRICS REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Total samples analyzed: {len(report_df)}\n\n")
            
            f.write("Read Statistics:\n")
            f.write(f"  Total reads (mean): {report_df['total_reads'].mean():.0f}\n")
            f.write(f"  Total reads (median): {report_df['total_reads'].median():.0f}\n")
            f.write(f"  Total reads (range): [{report_df['total_reads'].min():.0f}, {report_df['total_reads'].max():.0f}]\n\n")
            
            f.write("Duplication Metrics:\n")
            f.write(f"  Mean duplication rate: {report_df['duplication_rate'].mean():.1f}%\n")
            f.write(f"  Median duplication rate: {report_df['duplication_rate'].median():.1f}%\n")
            f.write(f"  Samples with <5% duplication: {(report_df['duplication_rate'] < 5).sum()}\n")
            f.write(f"  Samples with 5-10% duplication: {((report_df['duplication_rate'] >= 5) & (report_df['duplication_rate'] < 10)).sum()}\n")
            f.write(f"  Samples with >10% duplication: {(report_df['duplication_rate'] >= 10).sum()}\n\n")
            
            f.write("Complexity Metrics:\n")
            f.write(f"  Mean complexity score: {report_df['complexity_score'].mean():.1f}%\n")
            f.write(f"  Median complexity score: {report_df['complexity_score'].median():.1f}%\n")
            f.write(f"  Mean entropy: {report_df['entropy'].mean():.2f}\n")
            f.write(f"  Median entropy: {report_df['entropy'].median():.2f}\n\n")
            
            f.write("GC Content:\n")
            f.write(f"  Mean GC%: {report_df['gc_content'].mean():.1f}%\n")
            f.write(f"  Median GC%: {report_df['gc_content'].median():.1f}%\n")
            f.write(f"  Std Dev: {report_df['gc_content'].std():.1f}%\n\n")
            
            f.write("Library Size Estimates (Chao1):\n")
            f.write(f"  Mean: {report_df['chao1_estimate'].mean():.0f}\n")
            f.write(f"  Median: {report_df['chao1_estimate'].median():.0f}\n")
            f.write(f"  Range: [{report_df['chao1_estimate'].min():.0f}, {report_df['chao1_estimate'].max():.0f}]\n\n")
            
            f.write("Coverage Estimates:\n")
            f.write(f"  Mean coverage: {report_df['coverage_estimate'].mean():.1f}%\n")
            f.write(f"  Median coverage: {report_df['coverage_estimate'].median():.1f}%\n\n")
            
            f.write("Quality Assessment Distribution:\n")
            quality_counts = report_df['quality_assessment'].value_counts()
            for quality, count in quality_counts.items():
                pct = (count / len(report_df)) * 100
                f.write(f"  {quality}: {count} ({pct:.1f}%)\n")
        
        logger.info(f"Summary saved to {summary_file}")
    
    def generate_quality_flags(self, report_df: pd.DataFrame) -> pd.DataFrame:
        """Generate quality flags for each sample"""
        logger.info("Generating quality flags...")
        
        flags = []
        
        for idx, row in report_df.iterrows():
            sample_flags = []
            
            # Check duplication rate
            if row['duplication_rate'] > 20:
                sample_flags.append('HIGH_DUPLICATION')
            elif row['duplication_rate'] > 10:
                sample_flags.append('MODERATE_DUPLICATION')
            
            # Check complexity
            if row['complexity_score'] < 50:
                sample_flags.append('LOW_COMPLEXITY')
            
            # Check GC content (should be ~50% for human)
            if row['gc_content'] < 40 or row['gc_content'] > 60:
                sample_flags.append('ABNORMAL_GC_CONTENT')
            
            # Check coverage
            if row['coverage_estimate'] < 80:
                sample_flags.append('LOW_COVERAGE')
            
            flags.append('|'.join(sample_flags) if sample_flags else 'PASS')
        
        report_df['quality_flags'] = flags
        report_df.to_csv(self.output_dir / 'library_complexity_with_flags.csv', index=False)
        
        logger.info("Quality flags generated")
        return report_df


def main():
    """Main execution"""
    logger.info("Starting library complexity analysis...")
    
    metrics = LibraryComplexityMetrics(
        fastq_dir='fastq_downloads',
        output_dir='complexity_reports',
        max_workers=4
    )
    
    # Run analysis
    metrics.run_complexity_analysis()
    
    # Generate reports
    report_df = metrics.generate_complexity_report()
    metrics.generate_summary(report_df)
    report_df = metrics.generate_quality_flags(report_df)
    
    logger.info("Library complexity analysis complete!")
    
    return metrics, report_df


if __name__ == "__main__":
    metrics, report_df = main()
    
    print("\n" + "=" * 80)
    print("LIBRARY COMPLEXITY SUMMARY")
    print("=" * 80)
    print(f"\nTotal samples: {len(report_df)}")
    print(f"\nDuplication Rate: {report_df['duplication_rate'].mean():.1f}% ± {report_df['duplication_rate'].std():.1f}%")
    print(f"Complexity Score: {report_df['complexity_score'].mean():.1f}% ± {report_df['complexity_score'].std():.1f}%")
    print(f"GC Content: {report_df['gc_content'].mean():.1f}% ± {report_df['gc_content'].std():.1f}%")
    print(f"Entropy: {report_df['entropy'].mean():.2f} ± {report_df['entropy'].std():.2f}")
    print(f"\nLibrary Size (Chao1): {report_df['chao1_estimate'].mean():.0f} ± {report_df['chao1_estimate'].std():.0f}")
    print(f"Coverage: {report_df['coverage_estimate'].mean():.1f}% ± {report_df['coverage_estimate'].std():.1f}%")
    print(f"\nQuality Flags:")
    print(report_df['quality_flags'].value_counts())