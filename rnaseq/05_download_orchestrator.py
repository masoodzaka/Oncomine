#!/usr/bin/env python3
"""
Download Orchestrator for FASTQ Files
Manages downloads from SRA, ENA, and GEO with parallel processing, retry logic, and checksum verification
"""

import subprocess
import pandas as pd
import os
import hashlib
import logging
from typing import List, Dict, Tuple
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import requests
from urllib.parse import urljoin

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('download_orchestrator.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

class FASTQDownloadOrchestrator:
    """Orchestrate FASTQ downloads from multiple sources"""
    
    def __init__(self, metadata_file: str = 'consolidated_metadata.csv',
                 output_dir: str = 'fastq_downloads',
                 max_workers: int = 4):
        self.metadata_df = pd.read_csv(metadata_file)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.download_log = []
        self.failed_downloads = []
    
    def create_directory_structure(self):
        """Create organized directory structure by cancer type"""
        logger.info("Creating directory structure...")
        
        cancer_types = self.metadata_df['cancer_type'].dropna().unique()
        
        for cancer_type in cancer_types:
            # Skip NaN or non-string values
            if not isinstance(cancer_type, str):
                logger.warning(f"Skipping invalid cancer type: {cancer_type} (type: {type(cancer_type)})")
                continue
            
            cancer_dir = self.output_dir / str(cancer_type)
            cancer_dir.mkdir(parents=True, exist_ok=True)
            
            # Create subdirectories for PE and SE
            (cancer_dir / 'paired-end').mkdir(exist_ok=True)
            (cancer_dir / 'single-end').mkdir(exist_ok=True)
            
            logger.info(f"Created directory structure for {cancer_type}")
    
    def download_from_sra(self, run_accession: str, output_path: Path, 
                         max_retries: int = 3) -> Tuple[bool, str]:
        """Download FASTQ from SRA using parallel-fastq-dump"""
        logger.info(f"Downloading SRA run: {run_accession}")
        
        for attempt in range(max_retries):
            try:
                # Use parallel-fastq-dump for faster downloads
                cmd = [
                    'parallel-fastq-dump',
                    '--sra-id', run_accession,
                    '--outdir', str(output_path),
                    '--split-files',
                    '--gzip',
                    '--threads', '4',
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                
                if result.returncode == 0:
                    logger.info(f"Successfully downloaded {run_accession}")
                    return True, f"Downloaded {run_accession}"
                else:
                    logger.warning(f"Attempt {attempt + 1} failed for {run_accession}: {result.stderr}")
                    
            except subprocess.TimeoutExpired:
                logger.error(f"Timeout downloading {run_accession}")
            except FileNotFoundError:
                logger.error("parallel-fastq-dump not found. Install with: conda install -c bioconda parallel-fastq-dump")
                return False, "parallel-fastq-dump not installed"
            except Exception as e:
                logger.error(f"Error downloading {run_accession}: {e}")
            
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt  # Exponential backoff
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
        
        return False, f"Failed to download {run_accession} after {max_retries} attempts"
    
    def download_from_ena(self, run_accession: str, output_path: Path,
                         max_retries: int = 3) -> Tuple[bool, str]:
        """Download FASTQ from ENA FTP"""
        logger.info(f"Downloading ENA run: {run_accession}")
        
        # ENA FTP URL pattern
        ena_ftp_base = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
        
        # Construct ENA FTP path
        run_prefix = run_accession[:6]
        ena_ftp_path = f"{ena_ftp_base}/{run_prefix}/{run_accession}"
        
        for attempt in range(max_retries):
            try:
                # Use wget for FTP downloads
                cmd = [
                    'wget',
                    '-r',
                    '-nd',
                    '-P', str(output_path),
                    '--accept', '*.fastq.gz',
                    ena_ftp_path
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                
                if result.returncode == 0:
                    logger.info(f"Successfully downloaded {run_accession} from ENA")
                    return True, f"Downloaded {run_accession} from ENA"
                else:
                    logger.warning(f"Attempt {attempt + 1} failed for {run_accession}: {result.stderr}")
                    
            except subprocess.TimeoutExpired:
                logger.error(f"Timeout downloading {run_accession} from ENA")
            except Exception as e:
                logger.error(f"Error downloading {run_accession} from ENA: {e}")
            
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                logger.info(f"Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
        
        return False, f"Failed to download {run_accession} from ENA after {max_retries} attempts"
    
    def download_from_geo(self, geo_accession: str, output_path: Path,
                         max_retries: int = 3) -> Tuple[bool, str]:
        """Download supplementary files from GEO"""
        logger.info(f"Downloading GEO dataset: {geo_accession}")
        
        geo_url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={geo_accession}&format=file"
        
        for attempt in range(max_retries):
            try:
                cmd = [
                    'wget',
                    '-P', str(output_path),
                    geo_url
                ]
                
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                
                if result.returncode == 0:
                    logger.info(f"Successfully downloaded {geo_accession} from GEO")
                    return True, f"Downloaded {geo_accession} from GEO"
                else:
                    logger.warning(f"Attempt {attempt + 1} failed for {geo_accession}")
                    
            except Exception as e:
                logger.error(f"Error downloading {geo_accession}: {e}")
            
            if attempt < max_retries - 1:
                wait_time = 2 ** attempt
                time.sleep(wait_time)
        
        return False, f"Failed to download {geo_accession} after {max_retries} attempts"
    
    def calculate_checksum(self, file_path: Path, algorithm: str = 'md5') -> str:
        """Calculate file checksum"""
        hash_obj = hashlib.new(algorithm)
        
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b''):
                hash_obj.update(chunk)
        
        return hash_obj.hexdigest()
    
    def verify_downloads(self, output_path: Path) -> Dict[str, bool]:
        """Verify downloaded FASTQ files"""
        logger.info(f"Verifying downloads in {output_path}")
        
        verification_results = {}
        
        for fastq_file in output_path.glob('*.fastq.gz'):
            try:
                # Check file size (should be > 0)
                file_size = fastq_file.stat().st_size
                
                if file_size > 0:
                    # Calculate checksum
                    checksum = self.calculate_checksum(fastq_file)
                    verification_results[fastq_file.name] = {
                        'status': 'pass',
                        'size': file_size,
                        'checksum': checksum
                    }
                    logger.info(f"Verified {fastq_file.name}: {file_size} bytes")
                else:
                    verification_results[fastq_file.name] = {
                        'status': 'fail',
                        'reason': 'empty_file'
                    }
                    logger.warning(f"Empty file: {fastq_file.name}")
                    
            except Exception as e:
                verification_results[fastq_file.name] = {
                    'status': 'fail',
                    'reason': str(e)
                }
                logger.error(f"Error verifying {fastq_file.name}: {e}")
        
        return verification_results
    
    def download_sample(self, row: pd.Series) -> Dict:
        """Download a single sample"""
        database = row['database']
        accession = row['accession']
        cancer_type = row['cancer_type']
        seq_type = row.get('sequencing_type', 'unknown')
        
        # Determine output directory
        output_path = self.output_dir / cancer_type / seq_type
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Download based on database
        success = False
        message = ""
        
        try:
            if database == 'SRA':
                success, message = self.download_from_sra(accession, output_path)
            elif database == 'ENA':
                success, message = self.download_from_ena(accession, output_path)
            elif database == 'GEO':
                success, message = self.download_from_geo(accession, output_path)
            else:
                message = f"Unknown database: {database}"
        
        except Exception as e:
            success = False
            message = str(e)
        
        # Verify downloads
        verification = self.verify_downloads(output_path) if success else {}
        
        result = {
            'accession': accession,
            'database': database,
            'cancer_type': cancer_type,
            'success': success,
            'message': message,
            'verification': verification
        }
        
        if not success:
            self.failed_downloads.append(result)
        
        self.download_log.append(result)
        return result
    
    def orchestrate_downloads(self, sample_limit: int = None):
        """Orchestrate parallel downloads"""
        logger.info(f"Starting download orchestration with {self.max_workers} workers")
        
        # Create directory structure
        self.create_directory_structure()
        
        # Prepare download queue
        download_queue = self.metadata_df.copy()
        
        if sample_limit:
            download_queue = download_queue.head(sample_limit)
        
        logger.info(f"Queued {len(download_queue)} samples for download")
        
        # Execute parallel downloads
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self.download_sample, row): idx
                for idx, (_, row) in enumerate(download_queue.iterrows())
            }
            
            completed = 0
            for future in as_completed(futures):
                try:
                    result = future.result()
                    completed += 1
                    status = "✓" if result['success'] else "✗"
                    logger.info(f"[{completed}/{len(download_queue)}] {status} {result['accession']}")
                except Exception as e:
                    logger.error(f"Error in download task: {e}")
        
        logger.info(f"Download orchestration complete. {len(self.download_log)} total, {len(self.failed_downloads)} failed")
    
    def save_download_report(self, output_file: str = 'download_report.csv'):
        """Save download report"""
        report_df = pd.DataFrame(self.download_log)
        report_df.to_csv(output_file, index=False)
        logger.info(f"Download report saved to {output_file}")
        
        return report_df
    
    def save_failed_downloads(self, output_file: str = 'failed_downloads.txt'):
        """Save list of failed downloads for retry"""
        with open(output_file, 'w') as f:
            f.write("Failed Downloads - Retry List\n")
            f.write("=" * 50 + "\n\n")
            
            for item in self.failed_downloads:
                f.write(f"Accession: {item['accession']}\n")
                f.write(f"Database: {item['database']}\n")
                f.write(f"Message: {item['message']}\n")
                f.write("-" * 50 + "\n")
        
        logger.info(f"Failed downloads list saved to {output_file}")


def main():
    """Main execution"""
    logger.info("Starting FASTQ download orchestration...")
    
    orchestrator = FASTQDownloadOrchestrator(
        metadata_file='consolidated_metadata.csv',
        output_dir='fastq_downloads',
        max_workers=4
    )
    
    # Start downloads (limit to 10 for testing)
    orchestrator.orchestrate_downloads(sample_limit=10)
    
    # Save reports
    report_df = orchestrator.save_download_report()
    orchestrator.save_failed_downloads()
    
    logger.info("Download orchestration complete!")
    
    return orchestrator


if __name__ == "__main__":
    orchestrator = main()
    
    print("\n" + "=" * 60)
    print("DOWNLOAD SUMMARY")
    print("=" * 60)
    print(f"Total downloads: {len(orchestrator.download_log)}")
    print(f"Successful: {sum(1 for d in orchestrator.download_log if d['success'])}")
    print(f"Failed: {len(orchestrator.failed_downloads)}")