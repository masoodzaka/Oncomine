#!/usr/bin/env python3
"""
FastQC Quality Control Pipeline
Generates FastQC reports, MultiQC aggregation, and quality metrics for all FASTQ files
"""

import subprocess
import pandas as pd
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('fastqc_qc.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class FastQCQualityControl:
    """Execute FastQC and quality control analysis"""
    
    def __init__(self, fastq_dir: str = 'fastq_downloads',
                 output_dir: str = 'qc_reports',
                 max_workers: int = 4):
        self.fastq_dir = Path(fastq_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.qc_results = []
        self.failed_qc = []
    
    def find_fastq_files(self) -> List[Path]:
        """Find all FASTQ files in directory"""
        logger.info(f"Searching for FASTQ files in {self.fastq_dir}")
        
        fastq_files = list(self.fastq_dir.rglob('*.fastq.gz'))
        logger.info(f"Found {len(fastq_files)} FASTQ files")
        
        return fastq_files
    
    def run_fastqc(self, fastq_file: Path, output_dir: Path) -> Tuple[bool, str]:
        """Run FastQC on a single FASTQ file"""
        logger.info(f"Running FastQC on {fastq_file.name}")
        
        try:
            cmd = [
                'fastqc',
                '--outdir', str(output_dir),
                '--threads', '2',
                '--nogroup',
                str(fastq_file)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0:
                logger.info(f"FastQC completed for {fastq_file.name}")
                return True, f"FastQC completed for {fastq_file.name}"
            else:
                logger.error(f"FastQC failed for {fastq_file.name}: {result.stderr}")
                return False, result.stderr
                
        except FileNotFoundError:
            logger.error("FastQC not found. Install with: conda install -c bioconda fastqc")
            return False, "FastQC not installed"
        except subprocess.TimeoutExpired:
            logger.error(f"FastQC timeout for {fastq_file.name}")
            return False, "Timeout"
        except Exception as e:
            logger.error(f"Error running FastQC: {e}")
            return False, str(e)
    
    def extract_fastqc_metrics(self, fastqc_html: Path) -> Dict:
        """Extract key metrics from FastQC HTML report"""
        logger.info(f"Extracting metrics from {fastqc_html.name}")
        
        metrics = {
            'file': fastqc_html.stem.replace('_fastqc', ''),
            'total_sequences': 0,
            'sequence_length': '',
            'gc_content': 0,
            'per_base_quality': 'unknown',
            'per_sequence_quality': 'unknown',
            'adapter_content': 'unknown',
            'overrepresented_sequences': 'unknown',
            'status': 'unknown'
        }
        
        try:
            with open(fastqc_html, 'r') as f:
                content = f.read()
            
            # Extract total sequences
            match = re.search(r'Total Sequences</td><td>([0-9,]+)</td>', content)
            if match:
                metrics['total_sequences'] = int(match.group(1).replace(',', ''))
            
            # Extract sequence length
            match = re.search(r'Sequence length</td><td>([0-9\-]+)</td>', content)
            if match:
                metrics['sequence_length'] = match.group(1)
            
            # Extract GC content
            match = re.search(r'%GC</td><td>([0-9]+)</td>', content)
            if match:
                metrics['gc_content'] = int(match.group(1))
            
            # Extract status indicators
            if 'Per base sequence quality' in content:
                if '<span class="pass">' in content.split('Per base sequence quality')[1].split('</tr>')[0]:
                    metrics['per_base_quality'] = 'pass'
                elif '<span class="warn">' in content.split('Per base sequence quality')[1].split('</tr>')[0]:
                    metrics['per_base_quality'] = 'warn'
                else:
                    metrics['per_base_quality'] = 'fail'
            
            if 'Per sequence quality scores' in content:
                if '<span class="pass">' in content.split('Per sequence quality scores')[1].split('</tr>')[0]:
                    metrics['per_sequence_quality'] = 'pass'
                elif '<span class="warn">' in content.split('Per sequence quality scores')[1].split('</tr>')[0]:
                    metrics['per_sequence_quality'] = 'warn'
                else:
                    metrics['per_sequence_quality'] = 'fail'
            
            if 'Adapter Content' in content:
                if '<span class="pass">' in content.split('Adapter Content')[1].split('</tr>')[0]:
                    metrics['adapter_content'] = 'pass'
                elif '<span class="warn">' in content.split('Adapter Content')[1].split('</tr>')[0]:
                    metrics['adapter_content'] = 'warn'
                else:
                    metrics['adapter_content'] = 'fail'
            
            if 'Overrepresented sequences' in content:
                if '<span class="pass">' in content.split('Overrepresented sequences')[1].split('</tr>')[0]:
                    metrics['overrepresented_sequences'] = 'pass'
                elif '<span class="warn">' in content.split('Overrepresented sequences')[1].split('</tr>')[0]:
                    metrics['overrepresented_sequences'] = 'warn'
                else:
                    metrics['overrepresented_sequences'] = 'fail'
            
            # Determine overall status
            if metrics['per_base_quality'] == 'fail' or metrics['per_sequence_quality'] == 'fail':
                metrics['status'] = 'fail'
            elif metrics['per_base_quality'] == 'warn' or metrics['per_sequence_quality'] == 'warn':
                metrics['status'] = 'warn'
            else:
                metrics['status'] = 'pass'
                
        except Exception as e:
            logger.error(f"Error extracting metrics from {fastqc_html}: {e}")
        
        return metrics
    
    def run_multiqc(self, qc_dir: Path) -> Tuple[bool, str]:
        """Run MultiQC to aggregate FastQC reports"""
        logger.info(f"Running MultiQC on {qc_dir}")
        
        try:
            cmd = [
                'multiqc',
                '--outdir', str(self.output_dir),
                '--force',
                '--interactive',
                str(qc_dir)
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            
            if result.returncode == 0:
                logger.info("MultiQC completed successfully")
                return True, "MultiQC completed"
            else:
                logger.error(f"MultiQC failed: {result.stderr}")
                return False, result.stderr
                
        except FileNotFoundError:
            logger.error("MultiQC not found. Install with: conda install -c bioconda multiqc")
            return False, "MultiQC not installed"
        except Exception as e:
            logger.error(f"Error running MultiQC: {e}")
            return False, str(e)
    
    def validate_fastq_format(self, fastq_file: Path) -> Tuple[bool, Dict]:
        """Validate FASTQ file format integrity"""
        logger.info(f"Validating FASTQ format for {fastq_file.name}")
        
        validation = {
            'file': fastq_file.name,
            'valid': True,
            'total_reads': 0,
            'errors': []
        }
        
        try:
            import gzip
            
            read_count = 0
            with gzip.open(fastq_file, 'rt') as f:
                line_num = 0
                for line in f:
                    line_num += 1
                    
                    # FASTQ format: 4 lines per read
                    if line_num % 4 == 1:
                        if not line.startswith('@'):
                            validation['errors'].append(f"Line {line_num}: Invalid header (should start with @)")
                            validation['valid'] = False
                    elif line_num % 4 == 3:
                        if not line.startswith('+'):
                            validation['errors'].append(f"Line {line_num}: Invalid separator (should start with +)")
                            validation['valid'] = False
                    elif line_num % 4 == 0:
                        read_count += 1
            
            validation['total_reads'] = read_count
            
            if validation['valid']:
                logger.info(f"FASTQ validation passed for {fastq_file.name}: {read_count} reads")
            else:
                logger.warning(f"FASTQ validation failed for {fastq_file.name}")
                
        except Exception as e:
            validation['valid'] = False
            validation['errors'].append(str(e))
            logger.error(f"Error validating {fastq_file.name}: {e}")
        
        return validation['valid'], validation
    
    def process_fastq_file(self, fastq_file: Path) -> Dict:
        """Process a single FASTQ file"""
        result = {
            'file': fastq_file.name,
            'path': str(fastq_file),
            'fastqc_success': False,
            'format_valid': False,
            'metrics': {}
        }
        
        # Validate format
        valid, validation = self.validate_fastq_format(fastq_file)
        result['format_valid'] = valid
        result['validation'] = validation
        
        if not valid:
            logger.warning(f"Skipping FastQC for invalid file: {fastq_file.name}")
            self.failed_qc.append(result)
            return result
        
        # Run FastQC
        qc_output_dir = self.output_dir / 'fastqc_reports'
        qc_output_dir.mkdir(parents=True, exist_ok=True)
        
        success, message = self.run_fastqc(fastq_file, qc_output_dir)
        result['fastqc_success'] = success
        result['fastqc_message'] = message
        
        if success:
            # Extract metrics
            fastqc_html = qc_output_dir / f"{fastq_file.stem}_fastqc.html"
            if fastqc_html.exists():
                metrics = self.extract_fastqc_metrics(fastqc_html)
                result['metrics'] = metrics
        else:
            self.failed_qc.append(result)
        
        self.qc_results.append(result)
        return result
    
    def run_quality_control(self):
        """Execute quality control pipeline"""
        logger.info("Starting quality control pipeline...")
        
        # Find FASTQ files
        fastq_files = self.find_fastq_files()
        
        if not fastq_files:
            logger.warning("No FASTQ files found - skipping FastQC execution")
            logger.info("QC pipeline skipped (no input files)")
            return
        
        logger.info(f"Processing {len(fastq_files)} FASTQ files with {self.max_workers} workers")
        
        # Process files in parallel
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self.process_fastq_file, fastq_file): fastq_file
                for fastq_file in fastq_files
            }
            
            completed = 0
            for future in as_completed(futures):
                try:
                    result = future.result()
                    completed += 1
                    status = "✓" if result['fastqc_success'] else "✗"
                    logger.info(f"[{completed}/{len(fastq_files)}] {status} {result['file']}")
                except Exception as e:
                    logger.error(f"Error processing file: {e}")
        
        # Run MultiQC
        qc_dir = self.output_dir / 'fastqc_reports'
        if qc_dir.exists():
            self.run_multiqc(qc_dir)
        
        logger.info("Quality control pipeline complete!")
    
    def generate_qc_summary(self) -> pd.DataFrame:
        """Generate QC summary report"""
        logger.info("Generating QC summary report...")
        
        summary_data = []
        
        for result in self.qc_results:
            summary_row = {
                'file': result['file'],
                'format_valid': result['format_valid'],
                'total_reads': result['validation'].get('total_reads', 0) if result['format_valid'] else 0,
                'fastqc_success': result['fastqc_success'],
                'gc_content': result['metrics'].get('gc_content', 0),
                'sequence_length': result['metrics'].get('sequence_length', ''),
                'per_base_quality': result['metrics'].get('per_base_quality', 'unknown'),
                'per_sequence_quality': result['metrics'].get('per_sequence_quality', 'unknown'),
                'adapter_content': result['metrics'].get('adapter_content', 'unknown'),
                'overall_status': result['metrics'].get('status', 'unknown')
            }
            summary_data.append(summary_row)
        
        # Handle empty results - create empty DataFrame with expected columns
        if not summary_data:
            summary_data = [{
                'file': '',
                'format_valid': False,
                'total_reads': 0,
                'fastqc_success': False,
                'gc_content': 0,
                'sequence_length': '',
                'per_base_quality': 'unknown',
                'per_sequence_quality': 'unknown',
                'adapter_content': 'unknown',
                'overall_status': 'unknown'
            }]
            logger.warning("No QC results available - creating empty summary")
        
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(self.output_dir / 'qc_summary.csv', index=False)
        logger.info(f"QC summary saved to {self.output_dir / 'qc_summary.csv'}")
        
        return summary_df
    
    def generate_qc_report(self, summary_df: pd.DataFrame):
        """Generate comprehensive QC report"""
        logger.info("Generating comprehensive QC report...")
        
        report_file = self.output_dir / 'QC_REPORT.txt'
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("FASTQ QUALITY CONTROL REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            # Handle empty DataFrame (no files processed)
            if summary_df.empty or (len(summary_df) == 1 and summary_df.iloc[0]['file'] == ''):
                f.write("No FASTQ files were found for quality control analysis.\n")
                f.write("Please ensure FASTQ files have been downloaded to the fastq_downloads directory.\n\n")
                f.write("Expected directory structure:\n")
                f.write("  fastq_downloads/\n")
                f.write("    ├── cancer_type_1/\n")
                f.write("    │   ├── *.fastq.gz\n")
                f.write("    └── cancer_type_2/\n")
                f.write("        ├── *.fastq.gz\n")
                logger.warning("No FASTQ files found for QC analysis")
            else:
                # Normal report generation
                f.write(f"Total files processed: {len(summary_df)}\n")
                
                # Only access these columns if they exist and DataFrame is not empty
                if 'format_valid' in summary_df.columns:
                    f.write(f"Files with valid format: {summary_df['format_valid'].sum()}\n")
                if 'fastqc_success' in summary_df.columns:
                    f.write(f"Files with successful FastQC: {summary_df['fastqc_success'].sum()}\n")
                if 'overall_status' in summary_df.columns:
                    f.write(f"Files with PASS status: {(summary_df['overall_status'] == 'pass').sum()}\n")
                    f.write(f"Files with WARN status: {(summary_df['overall_status'] == 'warn').sum()}\n")
                    f.write(f"Files with FAIL status: {(summary_df['overall_status'] == 'fail').sum()}\n\n")
                
                f.write("Quality Metrics Summary:\n")
                f.write("-" * 80 + "\n")
                
                if 'gc_content' in summary_df.columns and not summary_df['gc_content'].empty:
                    f.write(f"Average GC content: {summary_df['gc_content'].mean():.1f}%\n")
                if 'total_reads' in summary_df.columns and not summary_df['total_reads'].empty:
                    f.write(f"Average read count: {summary_df['total_reads'].mean():.0f}\n")
                    f.write(f"Median read count: {summary_df['total_reads'].median():.0f}\n\n")
                
                f.write("Per-Base Quality Distribution:\n")
                if 'per_base_quality' in summary_df.columns:
                    f.write(f"  PASS: {(summary_df['per_base_quality'] == 'pass').sum()}\n")
                    f.write(f"  WARN: {(summary_df['per_base_quality'] == 'warn').sum()}\n")
                    f.write(f"  FAIL: {(summary_df['per_base_quality'] == 'fail').sum()}\n\n")
                
                f.write("Adapter Content Distribution:\n")
                if 'adapter_content' in summary_df.columns:
                    f.write(f"  PASS: {(summary_df['adapter_content'] == 'pass').sum()}\n")
                    f.write(f"  WARN: {(summary_df['adapter_content'] == 'warn').sum()}\n")
                    f.write(f"  FAIL: {(summary_df['adapter_content'] == 'fail').sum()}\n\n")
                
                if self.failed_qc:
                    f.write("Failed QC Files:\n")
                    f.write("-" * 80 + "\n")
                    for failed in self.failed_qc:
                        f.write(f"  {failed['file']}\n")
                        if failed['validation'].get('errors'):
                            for error in failed['validation']['errors'][:3]:
                                f.write(f"    - {error}\n")
        
        logger.info(f"QC report saved to {report_file}")


def main():
    """Main execution"""
    logger.info("Starting FastQC quality control...")
    
    qc = FastQCQualityControl(
        fastq_dir='fastq_downloads',
        output_dir='qc_reports',
        max_workers=4
    )
    
    # Run quality control
    qc.run_quality_control()
    
    # Generate reports
    summary_df = qc.generate_qc_summary()
    qc.generate_qc_report(summary_df)
    
    logger.info("Quality control complete!")
    
    return qc, summary_df


if __name__ == "__main__":
    qc, summary_df = main()
    
    print("\n" + "=" * 80)
    print("QUALITY CONTROL SUMMARY")
    print("=" * 80)
    
    # Handle empty results
    if summary_df.empty or (len(summary_df) == 1 and summary_df.iloc[0]['file'] == ''):
        print("\nNo FASTQ files were processed.")
        print("Status: SKIPPED - No input files found")
    else:
        print(f"\nTotal files: {len(summary_df)}")
        if 'format_valid' in summary_df.columns:
            print(f"Valid format: {summary_df['format_valid'].sum()}")
        if 'fastqc_success' in summary_df.columns:
            print(f"FastQC successful: {summary_df['fastqc_success'].sum()}")
        if 'overall_status' in summary_df.columns:
            print(f"\nStatus distribution:")
            print(summary_df['overall_status'].value_counts())
        if 'gc_content' in summary_df.columns:
            print(f"\nGC content: {summary_df['gc_content'].mean():.1f}% ± {summary_df['gc_content'].std():.1f}%")
    
    print("=" * 80)