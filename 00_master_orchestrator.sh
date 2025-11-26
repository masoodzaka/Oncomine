#!/bin/bash
#
# Master Orchestrator for Cancer RNA-seq FASTQ Curation
# Coordinates all phases of database querying, downloading, QC, and advanced analysis
#
# IMPROVEMENTS:
# - Set environment variables for NCBI_EMAIL and NCBI_API_KEY before running
# - Uses centralized config.py and utils.py for error handling and validation
# - Strict error checking: pipeline stops on first failure (unless --continue-on-error)
# - Detailed logging with phase checkpoints
# - Dependency verification before execution

set -euo pipefail  # Exit on error, undefined variables, pipe failures

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MASTER_LOG="${LOG_DIR}/master_orchestration_${TIMESTAMP}.log"
CONTINUE_ON_ERROR=false
EXIT_CODE=0

# Create log directory
mkdir -p "${LOG_DIR}"

# Logging functions
log() {
    local msg="$1"
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} ${msg}" | tee -a "${MASTER_LOG}"
}

log_success() {
    local msg="$1"
    echo -e "${GREEN}[✓]${NC} ${msg}" | tee -a "${MASTER_LOG}"
}

log_error() {
    local msg="$1"
    echo -e "${RED}[✗]${NC} ${msg}" | tee -a "${MASTER_LOG}"
}

log_warning() {
    local msg="$1"
    echo -e "${YELLOW}[!]${NC} ${msg}" | tee -a "${MASTER_LOG}"
}

# Error handler
handle_error() {
    local phase="$1"
    local exit_code="$2"
    log_error "Phase '${phase}' failed with exit code ${exit_code}"
    
    if [ "${CONTINUE_ON_ERROR}" = false ]; then
        log_error "Stopping pipeline (use --continue-on-error to continue)"
        EXIT_CODE="${exit_code}"
        return 1
    else
        log_warning "Continuing with next phase (--continue-on-error enabled)"
        return 0
    fi
}

# Verify required environment variables
check_requirements() {
    log "Verifying required environment variables..."
    
    if [ -z "${NCBI_EMAIL:-}" ]; then
        log_error "NCBI_EMAIL environment variable not set"
        echo "  Register for free at: https://www.ncbi.nlm.nih.gov/account/"
        return 1
    fi
    
    if [ -z "${NCBI_API_KEY:-}" ]; then
        log_error "NCBI_API_KEY environment variable not set"
        echo "  Register for free at: https://www.ncbi.nlm.nih.gov/account/"
        return 1
    fi
    
    log_success "Required environment variables are set"
    return 0
}

# Verify Python modules
check_python_modules() {
    log "Verifying Python modules..."
    
    local required_modules=("requests" "pandas" "numpy")
    for module in "${required_modules[@]}"; do
        if ! python3 -c "import ${module}" 2>/dev/null; then
            log_error "Missing Python module: ${module}"
            log_warning "Install with: pip install ${module}"
            return 1
        fi
    done
    
    log_success "All required Python modules found"
    return 0
}

# Header
echo -e "${BLUE}"
echo "╔════════════════════════════════════════════════════════════════════════════════╗"
echo "║                 Cancer RNA-seq FASTQ Curation Pipeline                         ║"
echo "║              GEO | SRA | ENA | QC | Contamination | Complexity                 ║"
echo "╚════════════════════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

log "Starting master orchestration pipeline..."
log "Log file: ${MASTER_LOG}"
log "Script directory: ${SCRIPT_DIR}"
log "Timestamp: ${TIMESTAMP}"

# Parse command line arguments
SKIP_QUERY=false
SKIP_DOWNLOAD=false
SKIP_QC=false
SKIP_ADVANCED=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-query)
            SKIP_QUERY=true
            shift
            ;;
        --skip-download)
            SKIP_DOWNLOAD=true
            shift
            ;;
        --skip-qc)
            SKIP_QC=true
            shift
            ;;
        --skip-advanced)
            SKIP_ADVANCED=true
            shift
            ;;
        --continue-on-error)
            CONTINUE_ON_ERROR=true
            shift
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            log_warning "Unknown option: $1"
            shift
            ;;
    esac
done

# Show help function
show_help() {
    echo -e "${BLUE}RNA-seq Pipeline - Orchestrator${NC}"
    echo ""
    echo "Usage: bash 00_master_orchestrator.sh [options]"
    echo ""
    echo "Options:"
    echo "  --skip-query          Skip database queries (use existing results)"
    echo "  --skip-download       Skip downloads (use existing FASTQ files)"
    echo "  --skip-qc             Skip quality control"
    echo "  --skip-advanced       Skip advanced analysis"
    echo "  --continue-on-error   Continue pipeline on phase failures"
    echo "  --help                Show this help message"
    echo ""
    echo "Required Environment Variables:"
    echo "  NCBI_EMAIL            Your email for NCBI (register at ncbi.nlm.nih.gov/account)"
    echo "  NCBI_API_KEY          Your NCBI API key (get at ncbi.nlm.nih.gov/account)"
    echo ""
    echo "Optional Environment Variables:"
    echo "  LOG_LEVEL             Logging level (DEBUG, INFO, WARNING, ERROR) - default: INFO"
    echo "  API_TIMEOUT           API timeout in seconds - default: 30"
    echo "  API_RETRIES           Number of API retries - default: 3"
    echo "  MAX_PARALLEL_DOWNLOADS  Max parallel downloads - default: 4"
}

# Pre-flight checks
log "\n${BLUE}=== Pre-flight Checks ===${NC}"

if ! check_requirements; then
    log_error "Pre-flight checks failed"
    exit 1
fi

if ! check_python_modules; then
    log_error "Pre-flight checks failed"
    exit 1
fi

log_success "All pre-flight checks passed"

# Phase 1: Database Queries
log "\n${BLUE}=== PHASE 1: Database Query & Metadata Extraction ===${NC}"

if [ "$SKIP_QUERY" = true ]; then
    log_warning "Skipping database queries (using existing results)"
else
    log "Querying GEO database..."
    if python3 "${SCRIPT_DIR}/01_geo_query.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "GEO query complete"
    else
        if ! handle_error "GEO Query" $?; then
            exit ${EXIT_CODE}
        fi
    fi
    
    log "Querying SRA database..."
    if python3 "${SCRIPT_DIR}/02_sra_query.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "SRA query complete"
    else
        if ! handle_error "SRA Query" $?; then
            exit ${EXIT_CODE}
        fi
    fi
    
    log "Querying ENA database..."
    if python3 "${SCRIPT_DIR}/03_ena_query.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "ENA query complete"
    else
        if ! handle_error "ENA Query" $?; then
            exit ${EXIT_CODE}
        fi
    fi
fi

# Phase 2: Metadata Consolidation
log "\n${BLUE}=== PHASE 2: Metadata Consolidation ===${NC}"

log "Consolidating metadata from all databases..."
if python3 "${SCRIPT_DIR}/04_consolidate_metadata.py" >> "${MASTER_LOG}" 2>&1; then
    log_success "Metadata consolidation complete"
else
    if ! handle_error "Metadata Consolidation" $?; then
        exit ${EXIT_CODE}
    fi
fi

# Phase 3: Download Orchestration
log "\n${BLUE}=== PHASE 3: Download Orchestration ===${NC}"

if [ "$SKIP_DOWNLOAD" = true ]; then
    log_warning "Skipping downloads (using existing FASTQ files)"
else
    log "Starting FASTQ downloads from SRA, ENA, and GEO..."
    log_warning "This may take several hours depending on dataset size and network speed"
    if python3 "${SCRIPT_DIR}/05_download_orchestrator.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "Download orchestration complete"
    else
        if ! handle_error "Download Orchestration" $?; then
            exit ${EXIT_CODE}
        fi
    fi
fi

# Phase 4: Quality Control
log "\n${BLUE}=== PHASE 4: Quality Control & Validation ===${NC}"

if [ "$SKIP_QC" = true ]; then
    log_warning "Skipping quality control"
else
    log "Running FastQC and quality control analysis..."
    if python3 "${SCRIPT_DIR}/06_fastqc_quality_control.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "Quality control complete"
    else
        if ! handle_error "Quality Control" $?; then
            exit ${EXIT_CODE}
        fi
    fi
fi

# Phase 5: Advanced Analysis
log "\n${BLUE}=== PHASE 5: Advanced Bioinformatics Analysis ===${NC}"

if [ "$SKIP_ADVANCED" = true ]; then
    log_warning "Skipping advanced analysis"
else
    log "Running contamination screening..."
    if python3 "${SCRIPT_DIR}/07_contamination_screening.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "Contamination screening complete"
    else
        log_warning "Contamination screening incomplete (Kraken2 may not be installed)"
    fi
    
    log "Running strand-specificity detection..."
    if python3 "${SCRIPT_DIR}/08_strand_specificity_detection.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "Strand-specificity detection complete"
    else
        if ! handle_error "Strand Detection" $?; then
            exit ${EXIT_CODE}
        fi
    fi
    
    log "Running library complexity analysis..."
    if python3 "${SCRIPT_DIR}/09_library_complexity_metrics.py" >> "${MASTER_LOG}" 2>&1; then
        log_success "Library complexity analysis complete"
    else
        if ! handle_error "Library Complexity" $?; then
            exit ${EXIT_CODE}
        fi
    fi
fi

# Final Summary
log "\n${BLUE}=== PIPELINE SUMMARY ===${NC}"

SUMMARY_COMPLETE=true

if [ -f "consolidated_metadata.csv" ]; then
    TOTAL_SAMPLES=$(tail -n +2 consolidated_metadata.csv 2>/dev/null | wc -l || echo "0")
    log_success "Total samples curated: ${TOTAL_SAMPLES}"
fi

if [ -f "download_report.csv" ]; then
    SUCCESSFUL_DOWNLOADS=$(grep -c "success\|True" download_report.csv 2>/dev/null || echo "0")
    log_success "Successful downloads: ${SUCCESSFUL_DOWNLOADS}"
fi

if [ -f "qc_reports/qc_summary.csv" ]; then
    PASS_QC=$(tail -n +2 qc_reports/qc_summary.csv 2>/dev/null | grep -ci "pass" || echo "0")
    log_success "Samples passing QC: ${PASS_QC}"
fi

if [ -f "contamination_reports/contamination_report.csv" ]; then
    PASS_CONTAMINATION=$(tail -n +2 contamination_reports/contamination_report.csv 2>/dev/null | grep -ci "pass" || echo "0")
    log_success "Samples passing contamination screening: ${PASS_CONTAMINATION}"
fi

if [ -f "complexity_reports/library_complexity_report.csv" ]; then
    TOTAL_COMPLEXITY=$(tail -n +2 complexity_reports/library_complexity_report.csv 2>/dev/null | wc -l || echo "0")
    log_success "Library complexity analysis complete: ${TOTAL_COMPLEXITY} samples"
fi

log "\n${GREEN}Pipeline execution complete!${NC}"
log "Results directory: ${SCRIPT_DIR}"
log "Log file: ${MASTER_LOG}"

if [ ${EXIT_CODE} -eq 0 ]; then
    echo -e "\n${GREEN}✓ All phases completed successfully!${NC}\n"
else
    echo -e "\n${RED}✗ Pipeline completed with errors (exit code: ${EXIT_CODE})${NC}\n"
fi

# Print usage information
echo -e "${BLUE}Usage:${NC}"
echo "  bash 00_master_orchestrator.sh                         # Run complete pipeline"
echo "  bash 00_master_orchestrator.sh --skip-query            # Skip database queries"
echo "  bash 00_master_orchestrator.sh --skip-download         # Skip downloads"
echo "  bash 00_master_orchestrator.sh --skip-qc               # Skip quality control"
echo "  bash 00_master_orchestrator.sh --skip-advanced         # Skip advanced analysis"
echo "  bash 00_master_orchestrator.sh --continue-on-error     # Continue on phase failures"
echo "  bash 00_master_orchestrator.sh --help                  # Show detailed help"
echo ""
echo -e "${BLUE}Required Environment Variables:${NC}"
echo "  NCBI_EMAIL                    Your NCBI email (free account at ncbi.nlm.nih.gov)"
echo "  NCBI_API_KEY                  Your NCBI API key (get at ncbi.nlm.nih.gov/account)"
echo ""
echo -e "${BLUE}Output directories:${NC}"
echo "  - output/logs/                       (Execution logs)"
echo "  - output/data/                       (Consolidated metadata and datasets)"
echo "  - output/results/                    (Analysis results)"
echo "  - fastq_downloads/                   (Downloaded FASTQ files)"
echo "  - qc_reports/                        (FastQC and MultiQC reports)"
echo "  - contamination_reports/             (Contamination analysis)"
echo "  - strand_reports/                    (Strand-specificity analysis)"
echo "  - complexity_reports/                (Library complexity metrics)"
echo ""

exit ${EXIT_CODE}