"""
Utility functions for RNA-seq pipeline.
Includes validation, retry logic, and error handling utilities.
"""

import logging
import time
from functools import wraps
from typing import Any, Callable, Optional, TypeVar
from pathlib import Path
import pandas as pd


logger = logging.getLogger(__name__)

T = TypeVar("T")


# ===== RETRY DECORATOR =====
def retry_with_backoff(
    max_attempts: int = 3,
    initial_delay: int = 1,
    backoff_factor: float = 2.0,
    exceptions: tuple = (Exception,)
) -> Callable:
    """
    Decorator to retry a function with exponential backoff.
    
    Args:
        max_attempts: Maximum number of retry attempts
        initial_delay: Initial delay between retries (seconds)
        backoff_factor: Multiply delay by this factor after each retry
        exceptions: Tuple of exceptions to catch and retry on
    
    Example:
        @retry_with_backoff(max_attempts=3, initial_delay=2, exceptions=(requests.RequestException,))
        def fetch_data(url):
            return requests.get(url)
    """
    def decorator(func: Callable[..., T]) -> Callable[..., T]:
        @wraps(func)
        def wrapper(*args, **kwargs) -> T:
            delay = initial_delay
            last_exception = None
            
            for attempt in range(1, max_attempts + 1):
                try:
                    logger.debug(f"Attempt {attempt}/{max_attempts} for {func.__name__}")
                    return func(*args, **kwargs)
                except exceptions as e:
                    last_exception = e
                    if attempt < max_attempts:
                        logger.warning(
                            f"{func.__name__} failed (attempt {attempt}/{max_attempts}): {str(e)}. "
                            f"Retrying in {delay}s..."
                        )
                        time.sleep(delay)
                        delay *= backoff_factor
                    else:
                        logger.error(
                            f"{func.__name__} failed after {max_attempts} attempts: {str(e)}"
                        )
            
            raise last_exception
        
        return wrapper
    return decorator


# ===== VALIDATION FUNCTIONS =====
def validate_file_exists(file_path: Optional[str], description: str = "File") -> Optional[Path]:
    """
    Validate that a file exists.
    
    Args:
        file_path: Path to file (can be None)
        description: Human-readable description for error message
    
    Returns:
        Path object if file exists, None if file_path is None
    
    Raises:
        FileNotFoundError: If file_path is provided but doesn't exist
    """
    if file_path is None:
        return None
    
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {file_path}")
    
    return path


def validate_csv_structure(csv_path: str, required_columns: list) -> pd.DataFrame:
    """
    Validate CSV file exists and contains required columns.
    
    Args:
        csv_path: Path to CSV file
        required_columns: List of column names that must exist
    
    Returns:
        Loaded DataFrame
    
    Raises:
        FileNotFoundError: If CSV doesn't exist
        ValueError: If required columns are missing or CSV is empty
    """
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"CSV file not found: {csv_path}")
    except pd.errors.EmptyDataError:
        raise ValueError(f"CSV file is empty: {csv_path}")
    
    if df.empty:
        raise ValueError(f"CSV file contains no data rows: {csv_path}")
    
    missing_columns = set(required_columns) - set(df.columns)
    if missing_columns:
        raise ValueError(
            f"CSV file {csv_path} missing required columns: {missing_columns}. "
            f"Found columns: {list(df.columns)}"
        )
    
    return df


def validate_positive_integer(value: Any, name: str = "Value") -> int:
    """
    Validate that value is a positive integer.
    
    Args:
        value: Value to validate
        name: Name for error message
    
    Returns:
        Integer value
    
    Raises:
        ValueError: If value is not a positive integer
    """
    try:
        int_val = int(value)
        if int_val <= 0:
            raise ValueError(f"{name} must be positive, got {int_val}")
        return int_val
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a positive integer, got {value}")


def validate_percentage(value: Any, name: str = "Value") -> float:
    """
    Validate that value is between 0 and 100.
    
    Args:
        value: Value to validate
        name: Name for error message
    
    Returns:
        Float value between 0-100
    
    Raises:
        ValueError: If value is not a valid percentage
    """
    try:
        float_val = float(value)
        if not 0 <= float_val <= 100:
            raise ValueError(f"{name} must be between 0-100, got {float_val}")
        return float_val
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a valid percentage (0-100), got {value}")


def validate_accession_format(accession: str, prefix: Optional[str] = None) -> str:
    """
    Validate that accession follows SRA/ENA format.
    
    Args:
        accession: Accession string to validate
        prefix: Optional expected prefix (e.g., 'SRR', 'ERR', 'GEO')
    
    Returns:
        Validated accession
    
    Raises:
        ValueError: If accession format is invalid
    """
    if not isinstance(accession, str) or not accession.strip():
        raise ValueError("Accession must be a non-empty string")
    
    accession = accession.strip().upper()
    
    # Check format: 3 letters + digits (SRR, ERR, etc.)
    if not (len(accession) >= 4 and accession[:3].isalpha() and accession[3:].isdigit()):
        raise ValueError(
            f"Invalid accession format: {accession}. "
            f"Expected format like SRR123456 or ERR987654"
        )
    
    if prefix and not accession.startswith(prefix):
        raise ValueError(
            f"Accession {accession} does not match expected prefix {prefix}"
        )
    
    return accession


def safe_divide(numerator: float, denominator: float, default: float = 0.0) -> float:
    """
    Safe division that handles zero denominator.
    
    Args:
        numerator: Dividend
        denominator: Divisor
        default: Value to return if denominator is zero
    
    Returns:
        numerator / denominator, or default if denominator is zero
    """
    if denominator == 0:
        logger.warning(f"Division by zero: {numerator}/{denominator}, returning {default}")
        return default
    return numerator / denominator


# ===== DATAFRAME UTILITIES =====
def safe_get_column(df: pd.DataFrame, column: str, default: Any = None) -> Any:
    """
    Safely get a column from DataFrame, returning default if not found.
    
    Args:
        df: DataFrame
        column: Column name
        default: Default value if column doesn't exist
    
    Returns:
        Column Series or default value
    """
    if column not in df.columns:
        logger.warning(f"Column '{column}' not found in DataFrame. Available: {list(df.columns)}")
        return default
    return df[column]


def consolidate_dataframes(
    *dataframes: pd.DataFrame,
    how: str = "outer",
    on_duplicates: str = "keep_first"
) -> pd.DataFrame:
    """
    Safely consolidate multiple DataFrames.
    
    Args:
        dataframes: Variable number of DataFrames to concatenate
        how: How to merge ('outer', 'inner', etc.)
        on_duplicates: How to handle duplicates ('keep_first', 'drop', etc.)
    
    Returns:
        Consolidated DataFrame
    
    Raises:
        ValueError: If no DataFrames provided or all are empty
    """
    if not dataframes:
        raise ValueError("No DataFrames provided for consolidation")
    
    # Filter out empty DataFrames
    valid_dfs = [df for df in dataframes if not df.empty]
    
    if not valid_dfs:
        raise ValueError("All provided DataFrames are empty")
    
    if len(valid_dfs) == 1:
        return valid_dfs[0].copy()
    
    # Concatenate with outer join to preserve all rows
    consolidated = pd.concat(valid_dfs, axis=0, ignore_index=True)
    
    # Remove duplicates if requested
    if on_duplicates == "drop":
        consolidated = consolidated.drop_duplicates()
        logger.info(f"Removed {len(valid_dfs) - len(consolidated)} duplicate rows")
    
    return consolidated


# ===== LOGGING UTILITIES =====
def setup_logger(name: str, log_file: Optional[str] = None, level: str = "INFO") -> logging.Logger:
    """
    Set up a logger with both console and optional file output.
    
    Args:
        name: Logger name
        log_file: Optional path to log file
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
    
    Returns:
        Configured logger
    """
    from config import LOG_FORMAT, LOG_DATE_FORMAT
    
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, level.upper(), logging.INFO))
    console_formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # File handler
    if log_file:
        try:
            Path(log_file).parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(log_file)
            file_handler.setLevel(getattr(logging, level.upper(), logging.INFO))
            file_formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        except Exception as e:
            logger.warning(f"Failed to create file handler for {log_file}: {e}")
    
    return logger


# ===== PROCESS MONITORING =====
def check_command_available(command: str) -> bool:
    """
    Check if a command is available in PATH.
    
    Args:
        command: Command name to check (e.g., 'fastqc', 'kraken2')
    
    Returns:
        True if command is available, False otherwise
    """
    import shutil
    result = shutil.which(command) is not None
    
    if not result:
        logger.error(f"Command not found in PATH: {command}")
    
    return result


def verify_dependencies(required_commands: list) -> bool:
    """
    Verify all required commands are available.
    
    Args:
        required_commands: List of command names to check
    
    Returns:
        True if all commands available, raises error if any missing
    
    Raises:
        RuntimeError: If any required command is not available
    """
    missing = []
    for cmd in required_commands:
        if not check_command_available(cmd):
            missing.append(cmd)
    
    if missing:
        raise RuntimeError(
            f"Missing required commands: {', '.join(missing)}. "
            f"Please install these tools before running the pipeline."
        )
    
    return True
