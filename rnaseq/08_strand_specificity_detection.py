#!/usr/bin/env python3
"""
Strand-Specificity Detection Pipeline
Determines library strandedness (unstranded, forward-stranded, reverse-stranded)
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
        logging.FileHandler('strand_specificity.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class StrandSpecificityDetection:
    """Detect library strand-specificity"""
    
    def __init__(self, fastq_dir: str = 'fastq_downloads',
                 output_dir: str = 'strand_reports',
                 reference_gtf: str = None,
                 reference_fasta: str = None,
                 max_workers: int = 4):
        self.fastq_dir = Path(fastq_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.reference_gtf = reference_gtf
        self.reference_fasta = reference_fasta
        self.max_workers = max_workers
        self.results = []
    
    def check_rseqc_installation(self) -> bool:
        """Check if RSeQC is installed"""
        logger.info("Checking RSeQC installation...")
        
        try:
            result = subprocess.run(['infer_experiment.py', '--version'],
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info("RSeQC found")
                return True
            else:
                logger.error("RSeQC not found")
                return False
        except FileNotFoundError:
            logger.error("RSeQC not installed. Install with: conda install -c bioconda rseqc")
            return False
        except Exception as e:
            logger.error(f"Error checking RSeQC: {e}")
            return False
    
    def run_rseqc_analysis(self, bam_file: Path, output_dir: Path) -> Dict:
        """Run RSeQC infer_experiment to detect strandedness"""
        logger.info(f"Running RSeQC on {bam_file.name}")
        
        result_dict = {
            'file': bam_file.name,
            'success': False,
            'strandedness': 'unknown',
            'forward_reads': 0.0,
            'reverse_reads': 0.0,
            'undetermined_reads': 0.0,
            'confidence': 'low'
        }
        
        try:
            output_file = output_dir / f"{bam_file.stem}_rseqc.txt"
            
            cmd = [
                'infer_experiment.py',
                '-r', self.reference_gtf,
                '-i', str(bam_file),
                '-o', str(output_file)
            ]
            
            process = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
            
            if process.returncode != 0:
                logger.error(f"RSeQC failed for {bam_file.name}: {process.stderr}")
                return result_dict
            
            # Parse RSeQC output
            result_dict = self._parse_rseqc_output(output_file)
            result_dict['success'] = True
            
            logger.info(f"RSeQC analysis complete: {result_dict['strandedness']}")
            
        except subprocess.TimeoutExpired:
            logger.error(f"RSeQC timeout for {bam_file.name}")
        except Exception as e:
            logger.error(f"Error running RSeQC: {e}")
        
        return result_dict
    
    def _parse_rseqc_output(self, output_file: Path) -> Dict:
        """Parse RSeQC output file"""
        result_dict = {
            'file': output_file.stem,
            'strandedness': 'unknown',
            'forward_reads': 0.0,
            'reverse_reads': 0.0,
            'undetermined_reads': 0.0,
            'confidence': 'low'
        }
        
        try:
            with open(output_file, 'r') as f:
                content = f.read()
            
            # Parse percentages
            for line in content.split('\n'):
                if 'Fraction of reads explained by "++,--"' in line:
                    # Unstranded
                    value = float(line.split(':')[1].strip())
                    result_dict['undetermined_reads'] = value
                elif 'Fraction of reads explained by "+-,-+"' in line:
                    # Reverse-stranded
                    value = float(line.split(':')[1].strip())
                    result_dict['reverse_reads'] = value
                elif 'Fraction of reads explained by "++-+--"' in line:
                    # Forward-stranded
                    value = float(line.split(':')[1].strip())
                    result_dict['forward_reads'] = value
            
            # Determine strandedness
            max_fraction = max(result_dict['forward_reads'],
                             result_dict['reverse_reads'],
                             result_dict['undetermined_reads'])
            
            if result_dict['forward_reads'] == max_fraction and result_dict['forward_reads'] > 0.7:
                result_dict['strandedness'] = 'forward-stranded'
                result_dict['confidence'] = 'high' if result_dict['forward_reads'] > 0.9 else 'medium'
            elif result_dict['reverse_reads'] == max_fraction and result_dict['reverse_reads'] > 0.7:
                result_dict['strandedness'] = 'reverse-stranded'
                result_dict['confidence'] = 'high' if result_dict['reverse_reads'] > 0.9 else 'medium'
            elif result_dict['undetermined_reads'] > 0.7:
                result_dict['strandedness'] = 'unstranded'
                result_dict['confidence'] = 'high' if result_dict['undetermined_reads'] > 0.9 else 'medium'
            else:
                result_dict['strandedness'] = 'ambiguous'
                result_dict['confidence'] = 'low'
            
        except Exception as e:
            logger.error(f"Error parsing RSeQC output: {e}")
        
        return result_dict
    
    def infer_strandedness_from_fastq(self, fastq_file: Path) -> Dict:
        """Infer strandedness from FASTQ using k-mer analysis"""
        logger.info(f"Inferring strandedness from {fastq_file.name}")
        
        result_dict = {
            'file': fastq_file.name,
            'method': 'kmer_analysis',
            'strandedness': 'unknown',
            'gc_skew': 0.0,
            'at_skew': 0.0,
            'confidence': 'low'
        }
        
        try:
            import gzip
            
            forward_gc = 0
            reverse_gc = 0
            forward_at = 0
            reverse_at = 0
            total_reads = 0
            
            with gzip.open(fastq_file, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    
                    # FASTQ sequence is on line 2, 6, 10, etc.
                    if line_num % 4 == 2:
                        sequence = line.strip().upper()
                        total_reads += 1
                        
                        # Count GC and AT
                        gc_count = sequence.count('G') + sequence.count('C')
                        at_count = sequence.count('A') + sequence.count('T')
                        
                        # Assume first half is forward, second half is reverse
                        if total_reads % 2 == 0:
                            forward_gc += gc_count
                            forward_at += at_count
                        else:
                            reverse_gc += gc_count
                            reverse_at += at_count
                        
                        if total_reads >= 100000:  # Sample first 100k reads
                            break
            
            # Calculate skew
            if forward_gc + reverse_gc > 0:
                gc_skew = (forward_gc - reverse_gc) / (forward_gc + reverse_gc)
                result_dict['gc_skew'] = gc_skew
            
            if forward_at + reverse_at > 0:
                at_skew = (forward_at - reverse_at) / (forward_at + reverse_at)
                result_dict['at_skew'] = at_skew
            
            # Infer strandedness based on skew
            skew_magnitude = abs(gc_skew)
            
            if skew_magnitude > 0.1:
                result_dict['strandedness'] = 'stranded'
                result_dict['confidence'] = 'high' if skew_magnitude > 0.2 else 'medium'
            else:
                result_dict['strandedness'] = 'unstranded'
                result_dict['confidence'] = 'high' if skew_magnitude < 0.05 else 'medium'
            
            logger.info(f"Strandedness inference complete: {result_dict['strandedness']} (skew: {gc_skew:.3f})")
            
        except Exception as e:
            logger.error(f"Error inferring strandedness: {e}")
        
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
            'rseqc_analysis': {},
            'kmer_analysis': {}
        }
        
        # Try RSeQC if BAM file available
        # (In practice, would need to align FASTQ to BAM first)
        
        # Run k-mer based inference
        kmer_result = self.infer_strandedness_from_fastq(fastq_file)
        result['kmer_analysis'] = kmer_result
        
        self.results.append(result)
        return result
    
    def run_strand_detection(self):
        """Execute strand-specificity detection pipeline"""
        logger.info("Starting strand-specificity detection pipeline...")
        
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
        
        logger.info("Strand-specificity detection complete!")
    
    def generate_strand_report(self) -> pd.DataFrame:
        """Generate strand-specificity report"""
        logger.info("Generating strand-specificity report...")
        
        report_data = []
        
        for result in self.results:
            kmer = result['kmer_analysis']
            
            report_row = {
                'file': result['file'],
                'strandedness': kmer.get('strandedness', 'unknown'),
                'gc_skew': kmer.get('gc_skew', 0.0),
                'at_skew': kmer.get('at_skew', 0.0),
                'confidence': kmer.get('confidence', 'low')
            }
            report_data.append(report_row)
        
        report_df = pd.DataFrame(report_data)
        report_df.to_csv(self.output_dir / 'strand_specificity_report.csv', index=False)
        logger.info("Strand-specificity report saved")
        
        return report_df
    
    def generate_summary(self, report_df: pd.DataFrame):
        """Generate summary statistics"""
        logger.info("Generating strand-specificity summary...")
        
        summary_file = self.output_dir / 'STRAND_SPECIFICITY_SUMMARY.txt'
        
        with open(summary_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("STRAND-SPECIFICITY DETECTION REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Total samples analyzed: {len(report_df)}\n\n")
            
            f.write("Strandedness Distribution:\n")
            strand_counts = report_df['strandedness'].value_counts()
            for strand, count in strand_counts.items():
                pct = (count / len(report_df)) * 100
                f.write(f"  {strand}: {count} ({pct:.1f}%)\n")
            
            f.write("\nConfidence Distribution:\n")
            conf_counts = report_df['confidence'].value_counts()
            for conf, count in conf_counts.items():
                pct = (count / len(report_df)) * 100
                f.write(f"  {conf}: {count} ({pct:.1f}%)\n")
            
            f.write("\nGC Skew Statistics:\n")
            f.write(f"  Mean: {report_df['gc_skew'].mean():.4f}\n")
            f.write(f"  Median: {report_df['gc_skew'].median():.4f}\n")
            f.write(f"  Std Dev: {report_df['gc_skew'].std():.4f}\n")
            f.write(f"  Range: [{report_df['gc_skew'].min():.4f}, {report_df['gc_skew'].max():.4f}]\n")
        
        logger.info(f"Summary saved to {summary_file}")


def main():
    """Main execution"""
    logger.info("Starting strand-specificity detection...")
    
    detection = StrandSpecificityDetection(
        fastq_dir='fastq_downloads',
        output_dir='strand_reports',
        reference_gtf=None,  # Update with actual GTF path
        reference_fasta=None,  # Update with actual FASTA path
        max_workers=4
    )
    
    # Run detection
    detection.run_strand_detection()
    
    # Generate reports
    report_df = detection.generate_strand_report()
    detection.generate_summary(report_df)
    
    logger.info("Strand-specificity detection complete!")
    
    return detection, report_df


if __name__ == "__main__":
    detection, report_df = main()
    
    print("\n" + "=" * 80)
    print("STRAND-SPECIFICITY SUMMARY")
    print("=" * 80)
    print(f"\nTotal samples: {len(report_df)}")
    print(f"\nStrandedness distribution:")
    print(report_df['strandedness'].value_counts())
    print(f"\nConfidence distribution:")
    print(report_df['confidence'].value_counts())
    print(f"\nGC Skew: {report_df['gc_skew'].mean():.4f} ± {report_df['gc_skew'].std():.4f}")