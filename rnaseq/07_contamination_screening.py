#!/usr/bin/env python3
"""
Contamination Screening Pipeline
Kraken2-based species composition analysis and k-mer contamination detection
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
        logging.FileHandler('contamination_screening.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class ContaminationScreening:
    """Contamination detection and species composition analysis"""
    
    def __init__(self, fastq_dir: str = 'fastq_downloads',
                 output_dir: str = 'contamination_reports',
                 kraken_db: str = None,
                 max_workers: int = 4):
        self.fastq_dir = Path(fastq_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.kraken_db = kraken_db or '/path/to/kraken2/db'
        self.max_workers = max_workers
        self.results = []
        self.contaminated_samples = []
    
    def check_kraken2_installation(self) -> bool:
        """Check if Kraken2 is installed"""
        logger.info("Checking Kraken2 installation...")
        
        try:
            result = subprocess.run(['kraken2', '--version'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"Kraken2 found: {result.stdout.strip()}")
                return True
            else:
                logger.error("Kraken2 not found")
                return False
        except FileNotFoundError:
            logger.error("Kraken2 not installed. Install with: conda install -c bioconda kraken2")
            return False
        except Exception as e:
            logger.error(f"Error checking Kraken2: {e}")
            return False
    
    def run_kraken2_analysis(self, fastq_file: Path, output_dir: Path) -> Tuple[bool, Dict]:
        """Run Kraken2 species composition analysis"""
        logger.info(f"Running Kraken2 on {fastq_file.name}")
        
        result_dict = {
            'file': fastq_file.name,
            'success': False,
            'species_composition': {},
            'contamination_level': 'unknown',
            'human_percentage': 0.0,
            'contaminant_species': []
        }
        
        try:
            # Kraken2 output files
            kraken_output = output_dir / f"{fastq_file.stem}_kraken.txt"
            kraken_report = output_dir / f"{fastq_file.stem}_kraken_report.txt"
            
            # Run Kraken2
            cmd = [
                'kraken2',
                '--db', self.kraken_db,
                '--output', str(kraken_output),
                '--report', str(kraken_report),
                '--threads', '4',
                '--gzip-compressed',
                str(fastq_file)
            ]
            
            process = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            
            if process.returncode != 0:
                logger.error(f"Kraken2 failed for {fastq_file.name}: {process.stderr}")
                return False, result_dict
            
            # Parse Kraken2 report
            species_data = self._parse_kraken_report(kraken_report)
            result_dict['species_composition'] = species_data
            
            # Calculate contamination metrics
            human_pct = species_data.get('Homo sapiens', 0.0)
            result_dict['human_percentage'] = human_pct
            
            # Identify contaminants (non-human species >1%)
            contaminants = {sp: pct for sp, pct in species_data.items() 
                          if sp != 'Homo sapiens' and pct > 1.0}
            result_dict['contaminant_species'] = contaminants
            
            # Classify contamination level
            if human_pct >= 95.0:
                result_dict['contamination_level'] = 'pass'
            elif human_pct >= 90.0:
                result_dict['contamination_level'] = 'warn'
            else:
                result_dict['contamination_level'] = 'fail'
            
            result_dict['success'] = True
            logger.info(f"Kraken2 analysis complete for {fastq_file.name}: {human_pct:.1f}% human")
            
            return True, result_dict
            
        except subprocess.TimeoutExpired:
            logger.error(f"Kraken2 timeout for {fastq_file.name}")
            return False, result_dict
        except Exception as e:
            logger.error(f"Error running Kraken2: {e}")
            return False, result_dict
    
    def _parse_kraken_report(self, report_file: Path) -> Dict[str, float]:
        """Parse Kraken2 report file"""
        species_data = {}
        
        try:
            with open(report_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        # Format: percentage, reads, reads_assigned, rank, taxid, name
                        percentage = float(parts[0])
                        rank = parts[3].strip()
                        name = parts[5].strip()
                        
                        # Only keep species-level classifications
                        if rank == 'S' and percentage > 0.1:
                            species_data[name] = percentage
            
            logger.info(f"Parsed {len(species_data)} species from {report_file.name}")
            
        except Exception as e:
            logger.error(f"Error parsing Kraken report: {e}")
        
        return species_data
    
    def kmer_contamination_detection(self, fastq_file: Path, output_dir: Path,
                                     kmer_size: int = 31) -> Dict:
        """Detect contamination using k-mer analysis"""
        logger.info(f"Running k-mer analysis on {fastq_file.name}")
        
        result_dict = {
            'file': fastq_file.name,
            'kmer_size': kmer_size,
            'unique_kmers': 0,
            'duplicate_kmers': 0,
            'kmer_complexity': 0.0,
            'potential_contamination': False
        }
        
        try:
            import gzip
            from collections import Counter
            
            kmer_counts = Counter()
            total_kmers = 0
            
            # Extract k-mers from FASTQ
            with gzip.open(fastq_file, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    
                    # FASTQ sequence is on line 2, 6, 10, etc. (every 4th line starting from 2)
                    if line_num % 4 == 2:
                        sequence = line.strip()
                        
                        # Extract k-mers
                        for i in range(len(sequence) - kmer_size + 1):
                            kmer = sequence[i:i+kmer_size]
                            kmer_counts[kmer] += 1
                            total_kmers += 1
            
            # Calculate metrics
            unique_kmers = len(kmer_counts)
            duplicate_kmers = sum(1 for count in kmer_counts.values() if count > 1)
            
            # K-mer complexity: unique k-mers / total k-mers
            kmer_complexity = unique_kmers / total_kmers if total_kmers > 0 else 0.0
            
            result_dict['unique_kmers'] = unique_kmers
            result_dict['duplicate_kmers'] = duplicate_kmers
            result_dict['kmer_complexity'] = kmer_complexity
            
            # Flag potential contamination if complexity is too low
            # (indicates repetitive/low-complexity sequences)
            if kmer_complexity < 0.7:
                result_dict['potential_contamination'] = True
                logger.warning(f"Low k-mer complexity detected in {fastq_file.name}: {kmer_complexity:.3f}")
            
            logger.info(f"K-mer analysis complete: {unique_kmers} unique, {kmer_complexity:.3f} complexity")
            
        except Exception as e:
            logger.error(f"Error in k-mer analysis: {e}")
        
        return result_dict
    
    def find_fastq_files(self) -> List[Path]:
        """Find all FASTQ files"""
        logger.info(f"Searching for FASTQ files in {self.fastq_dir}")
        fastq_files = list(self.fastq_dir.rglob('*.fastq.gz'))
        logger.info(f"Found {len(fastq_files)} FASTQ files")
        return fastq_files
    
    def process_sample(self, fastq_file: Path) -> Dict:
        """Process a single sample"""
        output_dir = self.output_dir / fastq_file.parent.name
        output_dir.mkdir(parents=True, exist_ok=True)
        
        result = {
            'file': fastq_file.name,
            'path': str(fastq_file),
            'kraken2_analysis': {},
            'kmer_analysis': {}
        }
        
        # Run Kraken2 if available
        if self.check_kraken2_installation():
            success, kraken_result = self.run_kraken2_analysis(fastq_file, output_dir)
            result['kraken2_analysis'] = kraken_result
            
            if kraken_result['contamination_level'] in ['warn', 'fail']:
                self.contaminated_samples.append(result)
        
        # Run k-mer analysis
        kmer_result = self.kmer_contamination_detection(fastq_file, output_dir)
        result['kmer_analysis'] = kmer_result
        
        self.results.append(result)
        return result
    
    def run_contamination_screening(self):
        """Execute contamination screening pipeline"""
        logger.info("Starting contamination screening pipeline...")
        
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
                    status = "✓" if result['kraken2_analysis'].get('success', False) else "✗"
                    logger.info(f"[{completed}/{len(fastq_files)}] {status} {result['file']}")
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
        
        logger.info("Contamination screening complete!")
    
    def generate_contamination_report(self) -> pd.DataFrame:
        """Generate contamination report"""
        logger.info("Generating contamination report...")
        
        report_data = []
        
        for result in self.results:
            kraken = result['kraken2_analysis']
            kmer = result['kmer_analysis']
            
            report_row = {
                'file': result['file'],
                'human_percentage': kraken.get('human_percentage', 0.0),
                'contamination_level': kraken.get('contamination_level', 'unknown'),
                'contaminant_species': str(kraken.get('contaminant_species', {})),
                'kmer_complexity': kmer.get('kmer_complexity', 0.0),
                'potential_kmer_contamination': kmer.get('potential_contamination', False)
            }
            report_data.append(report_row)
        
        report_df = pd.DataFrame(report_data)
        report_df.to_csv(self.output_dir / 'contamination_report.csv', index=False)
        logger.info(f"Contamination report saved")
        
        return report_df
    
    def generate_summary(self, report_df: pd.DataFrame):
        """Generate summary statistics"""
        logger.info("Generating contamination summary...")
        
        summary_file = self.output_dir / 'CONTAMINATION_SUMMARY.txt'
        
        with open(summary_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CONTAMINATION SCREENING REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Total samples analyzed: {len(report_df)}\n")
            f.write(f"Samples with PASS status: {(report_df['contamination_level'] == 'pass').sum()}\n")
            f.write(f"Samples with WARN status: {(report_df['contamination_level'] == 'warn').sum()}\n")
            f.write(f"Samples with FAIL status: {(report_df['contamination_level'] == 'fail').sum()}\n\n")
            
            f.write("Human Content Statistics:\n")
            f.write(f"  Mean: {report_df['human_percentage'].mean():.1f}%\n")
            f.write(f"  Median: {report_df['human_percentage'].median():.1f}%\n")
            f.write(f"  Min: {report_df['human_percentage'].min():.1f}%\n")
            f.write(f"  Max: {report_df['human_percentage'].max():.1f}%\n\n")
            
            f.write("K-mer Complexity Statistics:\n")
            f.write(f"  Mean: {report_df['kmer_complexity'].mean():.3f}\n")
            f.write(f"  Median: {report_df['kmer_complexity'].median():.3f}\n")
            f.write(f"  Samples with low complexity: {report_df['potential_kmer_contamination'].sum()}\n\n")
            
            if self.contaminated_samples:
                f.write("Contaminated Samples (Human <90%):\n")
                f.write("-" * 80 + "\n")
                for sample in self.contaminated_samples:
                    f.write(f"  {sample['file']}\n")
        
        logger.info(f"Summary saved to {summary_file}")


def main():
    """Main execution"""
    logger.info("Starting contamination screening...")
    
    screening = ContaminationScreening(
        fastq_dir='fastq_downloads',
        output_dir='contamination_reports',
        kraken_db='/path/to/kraken2/db',  # Update with actual path
        max_workers=4
    )
    
    # Run screening
    screening.run_contamination_screening()
    
    # Generate reports
    report_df = screening.generate_contamination_report()
    screening.generate_summary(report_df)
    
    logger.info("Contamination screening complete!")
    
    return screening, report_df


if __name__ == "__main__":
    screening, report_df = main()
    
    print("\n" + "=" * 80)
    print("CONTAMINATION SCREENING SUMMARY")
    print("=" * 80)
    print(f"\nTotal samples: {len(report_df)}")
    print(f"Pass: {(report_df['contamination_level'] == 'pass').sum()}")
    print(f"Warn: {(report_df['contamination_level'] == 'warn').sum()}")
    print(f"Fail: {(report_df['contamination_level'] == 'fail').sum()}")
    print(f"\nHuman content: {report_df['human_percentage'].mean():.1f}% ± {report_df['human_percentage'].std():.1f}%")