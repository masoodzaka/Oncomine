# Cancer RNA-seq FASTQ Curation Pipeline

## Overview

This comprehensive pipeline automates the curation of raw FASTQ files for bulk RNA-seq cancer datasets from public databases (GEO, SRA, ENA). It handles both paired-end and single-end sequencing data across multiple cancer indications including breast cancer subtypes, colorectal cancer (COAD), brain cancer, and lung cancer.

### Key Features

- **Multi-Database Integration**: Queries GEO, SRA, and ENA simultaneously
- **Intelligent Metadata Consolidation**: Merges metadata from all sources with duplicate detection
- **Parallel Download Orchestration**: Efficient multi-threaded downloads with retry logic
- **Comprehensive Quality Control**: FastQC + MultiQC analysis with detailed metrics
- **Cancer Subtype Stratification**: Automatic classification by cancer type and subtype
- **Reproducibility**: Full logging, checksum verification, and documentation
- **Production-Ready**: Error handling, rate limiting, and bandwidth optimization

## System Requirements

### Hardware
- **Minimum**: 8 CPU cores, 32GB RAM, 2TB storage
- **Recommended**: 16+ CPU cores, 64GB RAM, 5TB+ storage
- **Network**: High-speed internet connection (downloads can be 100GB+)

### Software Dependencies

```bash
# Core bioinformatics tools
conda install -c bioconda fastqc multiqc parallel-fastq-dump

# Python packages
pip install pandas numpy requests

# System utilities
# wget, curl, gzip (usually pre-installed)
```

### Installation

```bash
# 1. Clone or download the pipeline
cd rnaseq

# 2. Create conda environment
conda create -n rnaseq-curation python=3.11
conda activate rnaseq-curation

# 3. Install dependencies
conda install -c bioconda fastqc multiqc parallel-fastq-dump
pip install pandas numpy requests

# 4. Make scripts executable
chmod +x *.py *.sh

# 5. Configure NCBI API (optional but recommended)
# Add your NCBI API key to scripts for higher rate limits
```

## Configuration Management

### Centralized Configuration (`config.py`)

All pipeline parameters, thresholds, and search terms are managed in a single `config.py` file. This provides:

- **Environment Variable Support**: Override settings without editing code
- **Unified Query Terms**: Centralized database search queries
- **Consistent Filtering**: Standardized exclusion keywords across all databases
- **Easy Customization**: Change thresholds and parameters in one place

#### Key Configuration Parameters

**Cancer Type Search Terms**:
```python
CANCER_TYPES = {
    "breast": ["breast cancer", "brca", "mammary cancer"],
    "lung": ["lung cancer", "luad", "lusc", "nsclc", "sclc"],
    "colorectal": ["colorectal cancer", "crc", "colon cancer", "rectal cancer"],
    "prostate": ["prostate cancer"],
    "melanoma": ["melanoma", "skin cancer"],
    "pancreatic": ["pancreatic cancer"],
    "ovarian": ["ovarian cancer"],
}
```

**Database-Specific Search Terms**:
- `GEO_SEARCH_TERMS`: Optimized queries for Gene Expression Omnibus
- `SRA_SEARCH_TERMS`: Optimized queries for Sequence Read Archive
- `EXCLUDED_KEYWORDS`: Keywords for filtering single-cell and non-bulk studies

**Query Limits** (configurable via environment variables):
```bash
export MAX_GEO_RECORDS=5          # Maximum GEO datasets per cancer type
export MAX_SRA_RECORDS=5          # Maximum SRA experiments per cancer type
export MAX_ENA_RECORDS=5          # Maximum ENA runs per cancer type
```

**Quality Control Thresholds**:
```python
MIN_READ_LENGTH = 50              # Minimum read length (bp)
MIN_READS_PER_SAMPLE = 100000     # Minimum reads per sample
MIN_COMPLEXITY_SCORE = 0.7        # Library complexity threshold
MAX_DUPLICATION_RATE = 0.20       # Maximum duplication rate (20%)
```

---

## Pipeline Architecture

### Phase 1: Database Query & Metadata Extraction
**Scripts**: `01_geo_query.py`, `02_sra_query.py`, `03_ena_query.py`

Queries three major public databases for bulk RNA-seq cancer datasets using search terms from `config.py`:

- **GEO (Gene Expression Omnibus)**
  - Searches for GSE datasets with RNA-seq data using `GEO_SEARCH_TERMS`
  - Extracts sample metadata and platform information
  - Filters for bulk RNA-seq using `EXCLUDED_KEYWORDS`

- **SRA (Sequence Read Archive)**
  - Queries for RNA-seq experiments using `SRA_SEARCH_TERMS`
  - Extracts run information, read lengths, and base counts
  - Filters for Illumina platform, ≥50bp reads, ≥1GB data

- **ENA (European Nucleotide Archive)**
  - REST API queries for RNA-seq studies with cancer type classification
  - Retrieves run-level metadata with read statistics
  - Filters for quality thresholds (≥20M reads, ≥50bp)

**Centralized Features**:
- All three scripts use `CANCER_TYPES` from config for consistent categorization
- Search queries managed in config for easy updates without code changes
- Exclude keywords (`EXCLUDED_KEYWORDS`) applied consistently across all databases

**Output Files** (saved to `output/data/`):
- `geo_datasets.csv` - GEO dataset metadata
- `sra_experiments.csv` - SRA experiment metadata
- `ena_runs.csv` - ENA run metadata

### Phase 2: Metadata Consolidation
**Script**: `04_consolidate_metadata.py`

Merges metadata from all three databases into unified format:

- Standardizes field names and data types across all sources
- Automatically loads metadata from `output/data/` directory (configured in `config.py`)
- Detects and flags potential duplicate samples
- Adds quality control flags
- Stratifies samples by cancer type using consistent `CANCER_TYPES` from config
- Classifies sequencing type (paired-end vs single-end)

**Output Files** (saved to `output/data/`):
- `consolidated_metadata.csv` - Master metadata file
- `duplicate_report.txt` - Identified duplicate samples
- `consolidation_summary.txt` - Summary statistics

### Phase 3: Download Orchestration
**Script**: `05_download_orchestrator.py`

Manages parallel downloads from multiple sources:

- **SRA Downloads**: Uses `parallel-fastq-dump` for speed
- **ENA Downloads**: Direct FTP access via wget
- **GEO Downloads**: Supplementary file retrieval
- **Retry Logic**: Exponential backoff for failed downloads
- **Checksum Verification**: MD5/SHA256 validation
- **Directory Organization**: Organized by cancer type and sequencing type

**Features**:
- Configurable parallel workers (default: 4)
- Automatic retry with exponential backoff
- File integrity verification
- Detailed download logging
- Failed download tracking for retry

**Output Structure**:
```
fastq_downloads/
├── breast_cancer/
│   ├── paired-end/
│   │   ├── SRR*.fastq.gz
│   │   └── ...
│   └── single-end/
├── coad/
├── brain/
└── lung/
```

### Phase 4: Quality Control & Validation
**Script**: `06_fastqc_quality_control.py`

Comprehensive quality assessment of all FASTQ files:

- **FASTQ Format Validation**: Checks file integrity and format compliance
- **FastQC Analysis**: Per-base quality, adapter content, overrepresented sequences
- **MultiQC Aggregation**: Interactive HTML reports across all samples
- **Quality Metrics**: GC content, read length distribution, quality scores
- **QC Flagging**: Automatic pass/warn/fail classification

**Quality Thresholds**:
- Per-base quality: ≥Q30 (99.9% accuracy)
- Sequence length: ≥50bp
- Adapter content: <5% contamination
- Overrepresented sequences: <1% of reads

**Output Files**:
- `qc_reports/qc_summary.csv` - QC metrics for all samples
- `qc_reports/QC_REPORT.txt` - Comprehensive QC report
- `qc_reports/multiqc_report.html` - Interactive MultiQC report
- `qc_reports/fastqc_reports/` - Individual FastQC reports

## Usage

### Quick Start

```bash
# Run complete pipeline (requires NCBI_EMAIL and NCBI_API_KEY environment variables)
export NCBI_EMAIL="your.email@example.com"
export NCBI_API_KEY="your_api_key"
bash 00_master_orchestrator.sh

# Run with specific options
bash 00_master_orchestrator.sh --skip-query      # Skip database queries
bash 00_master_orchestrator.sh --skip-download   # Skip downloads
bash 00_master_orchestrator.sh --skip-qc         # Skip QC analysis
bash 00_master_orchestrator.sh --skip-advanced   # Skip advanced analysis
```

### Customizing Configuration

All pipeline behavior is controlled through `config.py` and environment variables:

```bash
# Adjust record limits for database queries
export MAX_GEO_RECORDS=10          # Increase from default 5
export MAX_SRA_RECORDS=10
export MAX_ENA_RECORDS=10

# Adjust quality thresholds
export MIN_READ_LENGTH=100         # Increase from default 50
export MAX_DUPLICATION_RATE=0.3    # Increase from default 0.20

# Adjust performance settings
export MAX_PARALLEL_DOWNLOADS=8    # Increase from default 4
export FASTQC_THREADS=8            # Increase from default 4

# Run pipeline with custom settings
bash 00_master_orchestrator.sh
```

### Individual Phase Execution

```bash
# Phase 1: Query databases (with custom record limits)
export MAX_GEO_RECORDS=20
python3 01_geo_query.py
python3 02_sra_query.py
python3 03_ena_query.py

# Phase 2: Consolidate metadata
python3 04_consolidate_metadata.py

# Phase 3: Download FASTQ files
python3 05_download_orchestrator.py

# Phase 4: Quality control
python3 06_fastqc_quality_control.py

# Phase 5: Advanced analysis
python3 07_contamination_screening.py
python3 08_strand_specificity_detection.py
python3 09_library_complexity_metrics.py
```

### Modifying Search Terms and Filters

To customize which studies are retrieved and filtered:

**Edit `config.py`** to modify:

```python
# Add new cancer types or keywords
CANCER_TYPES = {
    "breast": ["breast cancer", "brca", "mammary cancer"],
    # ... add more cancer types
}

# Customize GEO search queries
GEO_SEARCH_TERMS = {
    "breast": "(your custom search terms)",
    # ... customize per cancer type
}

# Customize SRA search queries
SRA_SEARCH_TERMS = {
    "breast": "(your custom search terms)",
    # ... customize per cancer type
}

# Modify exclusion keywords (filters out single-cell, etc.)
EXCLUDED_KEYWORDS = ['single cell', 'your keywords', '10x', ...]
```

After editing `config.py`, all scripts will automatically use the new settings on the next run.

## Output Directory Structure

All pipeline outputs are organized hierarchically under `output/` directory (configured in `config.py`):

```
output/
├── data/                                    # All CSV metadata files
│   ├── geo_datasets.csv                    # GEO query results
│   ├── sra_experiments.csv                 # SRA query results
│   ├── ena_runs.csv                        # ENA query results
│   ├── consolidated_metadata.csv           # Master metadata
│   └── duplicate_report.txt                # Duplicate detection results
│
├── logs/                                    # Pipeline execution logs
│   ├── master_orchestration_*.log          # Main orchestrator log
│   ├── geo_query.log
│   ├── sra_query.log
│   ├── ena_query.log
│   └── ...
│
└── results/                                 # Analysis results
    ├── qc_reports/                         # Quality control reports
    ├── contamination_reports/              # Contamination screening results
    ├── strand_reports/                     # Strand specificity analysis
    └── complexity_reports/                 # Library complexity metrics
```

**Key Configuration** (in `config.py`):
```python
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR", "./output"))
LOG_DIR = OUTPUT_DIR / "logs"
DATA_DIR = OUTPUT_DIR / "data"              # Where query results are saved
RESULTS_DIR = OUTPUT_DIR / "results"        # Where analysis results go
```

---

## Output Files & Interpretation

### Metadata Files

**consolidated_metadata.csv**
- `accession`: Database accession ID
- `database`: Source database (GEO/SRA/ENA)
- `cancer_type`: Cancer indication (from `CANCER_TYPES` config)
- `sequencing_type`: paired-end or single-end
- `read_length`: Average read length
- `base_count`: Total bases in sample
- `quality_flag`: QC status (pass/warn/fail)

### Quality Control Reports

**qc_summary.csv**
- `file`: FASTQ filename
- `total_reads`: Number of reads
- `gc_content`: GC percentage
- `per_base_quality`: Quality score status
- `adapter_content`: Adapter contamination status
- `overall_status`: Pass/Warn/Fail

**QC_REPORT.txt**
- Summary statistics across all samples
- Quality metric distributions
- Failed sample list with error details

**multiqc_report.html**
- Interactive visualization of all QC metrics
- Comparison across samples and cancer types
- Detailed quality score distributions

## Cancer Type Stratification

### Breast Cancer Subtypes
- ER+ (Estrogen Receptor Positive)
- PR+ (Progesterone Receptor Positive)
- HER2+ (Human Epidermal Growth Factor Receptor 2 Positive)
- TNBC (Triple Negative Breast Cancer)

### Colorectal Cancer (COAD)
- MSI-H (Microsatellite Instability High)
- MSS (Microsatellite Stable)

### Brain Cancer
- Glioblastoma (Grade IV)
- LGG (Lower Grade Glioma)
- Other brain tumors

### Lung Cancer
- LUAD (Lung Adenocarcinoma)
- LUSC (Lung Squamous Cell Carcinoma)

## Advanced Features

### Duplicate Detection

The pipeline identifies potential duplicate samples across databases using:
- Title/description matching
- Organism and cancer type matching
- Read characteristic similarity

Duplicates are flagged in `duplicate_report.txt` for manual review.

### Retry Logic

Failed downloads are automatically retried with exponential backoff:
- Attempt 1: Immediate retry
- Attempt 2: Wait 2 seconds
- Attempt 3: Wait 4 seconds

Failed downloads are logged in `failed_downloads.txt` for manual retry.

### Checksum Verification

All downloaded files are verified using MD5 checksums:
- Checksums stored in download report
- Corrupted files automatically flagged
- Enables data integrity validation

### Rate Limiting

Pipeline respects database rate limits:
- NCBI: ~3 requests/second
- ENA: Configurable rate limiting
- Automatic backoff on 429 (Too Many Requests) errors

## Troubleshooting

### Common Issues

**Issue**: "parallel-fastq-dump not found"
```bash
# Solution: Install SRA Toolkit
conda install -c bioconda sra-tools
```

**Issue**: "FastQC not found"
```bash
# Solution: Install FastQC
conda install -c bioconda fastqc
```

**Issue**: Network timeouts during downloads
```bash
# Solution: Increase timeout in download_orchestrator.py
timeout=7200  # 2 hours instead of 1 hour
```

**Issue**: Out of disk space
```bash
# Solution: Monitor disk usage
df -h
# Delete intermediate files if needed
rm -rf qc_reports/fastqc_reports/*.zip
```

### Performance Optimization

**For faster downloads**:
```python
max_workers=8  # Increase parallel workers (if network allows)
```

**For faster QC**:
```python
# Run FastQC on subset first
fastq_files = fastq_files[:100]  # Test on 100 files
```

**For memory efficiency**:
```python
# Process files in batches
batch_size = 50
for i in range(0, len(fastq_files), batch_size):
    batch = fastq_files[i:i+batch_size]
    # Process batch
```

## Data Organization Best Practices

### Directory Structure
```
project_root/
├── rnaseq/
│   ├── *.py (scripts)
│   ├── *.sh (orchestrators)
│   └── README.md
├── consolidated_metadata.csv
├── download_report.csv
├── fastq_downloads/
│   ├── breast_cancer/
│   ├── coad/
│   ├── brain/
│   └── lung/
├── qc_reports/
│   ├── qc_summary.csv
│   ├── QC_REPORT.txt
│   ├── multiqc_report.html
│   └── fastqc_reports/
└── logs/
    └── *.log
```

### Backup Strategy
```bash
# Backup metadata and reports
tar -czf metadata_backup_$(date +%Y%m%d).tar.gz \
    consolidated_metadata.csv \
    download_report.csv \
    qc_reports/

# Backup FASTQ files (optional - very large)
rsync -av fastq_downloads/ /backup/location/
```

## Citation & Attribution

When using this pipeline, please cite:
- GEO: Barrett et al. (2013) Nucleic Acids Res.
- SRA: Leinonen et al. (2011) Nucleic Acids Res.
- ENA: Cochrane et al. (2016) Nucleic Acids Res.
- FastQC: Andrews, S. (2010)
- MultiQC: Ewels et al. (2016) Bioinformatics

## License & Support

This pipeline is provided as-is for research purposes. For issues or questions:
1. Check the troubleshooting section
2. Review log files for detailed error messages
3. Verify all dependencies are installed
4. Test with a small subset of data first

## Version History

- **v1.0** (2025-11-24): Initial release
  - GEO, SRA, ENA integration
  - Parallel download orchestration
  - FastQC + MultiQC quality control
  - Comprehensive metadata curation

## Future Enhancements

- [ ] Kraken2 contamination screening
- [ ] Strand-specificity detection
- [ ] Library complexity metrics (Picard)
- [ ] Automated subtype classification
- [ ] Integration with Nextflow/Snakemake
- [ ] Cloud storage support (AWS S3, GCS)
- [ ] Real-time monitoring dashboard
- [ ] Machine learning-based QC prediction

---

**Last Updated**: November 24, 2025
**Maintainer**: Bioinformatics Team
**Status**: Production Ready