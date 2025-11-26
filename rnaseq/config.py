"""
Centralized configuration for RNA-seq pipeline.
Manages all hardcoded constants, thresholds, and environment-based settings.
"""

import os
from pathlib import Path


# ===== API CONFIGURATION =====
NCBI_EMAIL = os.getenv("NCBI_EMAIL", None)
NCBI_API_KEY = os.getenv("NCBI_API_KEY", None)
API_TIMEOUT = int(os.getenv("API_TIMEOUT", "30"))
API_RETRIES = int(os.getenv("API_RETRIES", "3"))
API_RETRY_DELAY = int(os.getenv("API_RETRY_DELAY", "5"))  # seconds

# Validate NCBI credentials - email is required, API key is optional
if not NCBI_EMAIL:
    raise ValueError(
        "NCBI_EMAIL environment variable must be set. "
        "Required by NCBI Entrez. "
        "For better performance, also set NCBI_API_KEY (3 requests/sec vs 100 requests/sec). "
        "Register at: https://www.ncbi.nlm.nih.gov/account/register/"
    )

if NCBI_API_KEY:
    import logging
    logging.getLogger(__name__).info("NCBI API key configured - using high rate limit (100 req/sec)")
else:
    import logging
    logging.getLogger(__name__).warning("NCBI API key not set - using standard rate limit (3 req/sec). Performance will be slower.")


# ===== DATABASE PATHS =====
KRAKEN2_DB_PATH = os.getenv(
    "KRAKEN2_DB_PATH",
    "/opt/kraken2/standard"
)
REFERENCE_GTF = os.getenv("REFERENCE_GTF", None)
REFERENCE_GENOME = os.getenv("REFERENCE_GENOME", None)


# ===== QUALITY CONTROL THRESHOLDS =====
# Complexity metrics (library complexity detection)
MIN_COMPLEXITY_SCORE = float(os.getenv("MIN_COMPLEXITY_SCORE", "0.7"))
MIN_GC_CONTENT_PERCENT = float(os.getenv("MIN_GC_CONTENT_PERCENT", "40"))
MAX_GC_CONTENT_PERCENT = float(os.getenv("MAX_GC_CONTENT_PERCENT", "60"))
MAX_DUPLICATION_RATE = float(os.getenv("MAX_DUPLICATION_RATE", "0.20"))

# Strand detection
STRAND_DETECT_SAMPLE_SIZE = int(os.getenv("STRAND_DETECT_SAMPLE_SIZE", "100000"))
FORWARD_STRAND_THRESHOLD = float(os.getenv("FORWARD_STRAND_THRESHOLD", "0.1"))
REVERSE_STRAND_THRESHOLD = float(os.getenv("REVERSE_STRAND_THRESHOLD", "0.2"))

# Contamination thresholds
CONTAMINATION_THRESHOLD = float(os.getenv("CONTAMINATION_THRESHOLD", "0.05"))  # 5%

# FastQC metrics
MIN_READ_LENGTH = int(os.getenv("MIN_READ_LENGTH", "50"))
MIN_READS_PER_SAMPLE = int(os.getenv("MIN_READS_PER_SAMPLE", "100000"))

# K-mer complexity
KMER_COMPLEXITY_THRESHOLD = float(os.getenv("KMER_COMPLEXITY_THRESHOLD", "0.7"))


# ===== QUERY CONFIGURATION =====
# Cancer types and search terms
CANCER_TYPES = {
    "breast": ["breast cancer", "brca", "mammary cancer"],
    "lung": ["lung cancer", "luad", "lusc", "nsclc", "sclc"],
    "colorectal": ["colorectal cancer", "crc", "colon cancer", "rectal cancer"],
    "prostate": ["prostate cancer"],
    "melanoma": ["melanoma", "skin cancer"],
    "pancreatic": ["pancreatic cancer"],
    "ovarian": ["ovarian cancer"],
}

# GEO database search query terms for each cancer type
GEO_SEARCH_TERMS = {
    "breast": '(breast cancer) AND (RNA-seq OR "RNA seq") AND (homo sapiens) AND (GSE)',
    "colorectal": '(colorectal cancer OR colon cancer OR COAD) AND (RNA-seq OR "RNA seq") AND (homo sapiens) AND (GSE)',
    "lung": '(lung cancer OR LUAD OR LUSC) AND (RNA-seq OR "RNA seq") AND (homo sapiens) AND (GSE)',
    "pancreatic": '(pancreatic cancer) AND (RNA-seq OR "RNA seq") AND (homo sapiens) AND (GSE)',
    "ovarian": '(ovarian cancer) AND (RNA-seq OR "RNA seq") AND (homo sapiens) AND (GSE)',
}

# SRA database search query terms for each cancer type
SRA_SEARCH_TERMS = {
    "breast": '(breast cancer) AND (RNA-seq) AND (homo sapiens) AND (ILLUMINA)',
    "colorectal": '(colorectal cancer OR colon cancer OR COAD) AND (RNA-seq) AND (homo sapiens) AND (ILLUMINA)',
    "lung": '(lung cancer OR LUAD OR LUSC) AND (RNA-seq) AND (homo sapiens) AND (ILLUMINA)',
    "pancreatic": '(pancreatic cancer) AND (RNA-seq) AND (homo sapiens) AND (ILLUMINA)',
    "ovarian": '(ovarian cancer) AND (RNA-seq) AND (homo sapiens) AND (ILLUMINA)',
}

# GEO/SRA/ENA query parameters
MAX_GEO_RECORDS = int(os.getenv("MAX_GEO_RECORDS", "5"))
MAX_SRA_RECORDS = int(os.getenv("MAX_SRA_RECORDS", "5"))
MAX_ENA_RECORDS = int(os.getenv("MAX_ENA_RECORDS", "5"))
MIN_SAMPLE_COUNT = int(os.getenv("MIN_SAMPLE_COUNT", "3"))
MAX_SAMPLE_COUNT = int(os.getenv("MAX_SAMPLE_COUNT", "1000"))


# ===== DOWNLOAD CONFIGURATION =====
MAX_DOWNLOAD_ATTEMPTS = int(os.getenv("MAX_DOWNLOAD_ATTEMPTS", "5"))
DOWNLOAD_TIMEOUT = int(os.getenv("DOWNLOAD_TIMEOUT", "3600"))  # 1 hour
MAX_PARALLEL_DOWNLOADS = int(os.getenv("MAX_PARALLEL_DOWNLOADS", "4"))
FASTQC_THREADS = int(os.getenv("FASTQC_THREADS", "4"))


# ===== PROCESSING CONFIGURATION =====
SAMPLE_LIMIT_TESTING = int(os.getenv("SAMPLE_LIMIT_TESTING", "10"))
ENABLE_SAMPLE_LIMIT = os.getenv("ENABLE_SAMPLE_LIMIT", "false").lower() == "true"


# ===== FILE PATHS =====
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR", "./output"))
LOG_DIR = OUTPUT_DIR / "logs"
DATA_DIR = OUTPUT_DIR / "data"
RESULTS_DIR = OUTPUT_DIR / "results"

# Create directories if they don't exist
for directory in [OUTPUT_DIR, LOG_DIR, DATA_DIR, RESULTS_DIR]:
    directory.mkdir(parents=True, exist_ok=True)


# ===== QUALITY FLAGS =====
QUALITY_FLAGS = {
    "PASS": "Sample passed all quality checks",
    "LOW_COMPLEXITY": "Library complexity below threshold",
    "HIGH_DUPLICATION": "Duplication rate exceeds threshold",
    "GC_BIAS": "GC content outside acceptable range",
    "CONTAMINATION": "Contamination detected above threshold",
    "LOW_READS": "Insufficient number of reads",
    "STRAND_AMBIGUOUS": "Strand orientation cannot be determined",
    "TRUNCATED": "FASTQ file appears truncated",
}


# ===== LOGGING CONFIGURATION =====
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO")
LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
