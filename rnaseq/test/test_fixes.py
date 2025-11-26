#!/usr/bin/env python3
"""
Test script to verify all critical fixes are working correctly.
Tests: response validation, JSON parsing, return statements, rate limiting.
"""

import os
import sys
import json
import time
import requests
from unittest.mock import Mock, patch, MagicMock
from io import StringIO

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from config import setup_logger, GEO_EMAIL, GEO_API_KEY
    from utils import retry_with_backoff, validate_accession_format
    from rnaseq.one_geo_query import GEOQueryEngine  # Import after path setup
except ImportError as e:
    print(f"ERROR: Failed to import required modules: {e}")
    print("Make sure config.py and utils.py are in the same directory")
    sys.exit(1)

logger = setup_logger(__name__)

class TestGEOFixes:
    """Test suite for GEO query fixes"""
    
    def __init__(self):
        self.passed = 0
        self.failed = 0
        self.engine = GEOQueryEngine()
    
    def test_empty_response_handling(self):
        """Test that empty API responses are handled gracefully"""
        test_name = "Empty Response Handling"
        try:
            with patch('requests.get') as mock_get:
                # Simulate empty response
                mock_response = Mock()
                mock_response.text = ''
                mock_response.status_code = 200
                mock_response.raise_for_status.return_value = None
                mock_get.return_value = mock_response
                
                # Should return empty list, not crash
                result = self.engine.search_geo('breast cancer', max_results=10)
                
                if result == []:
                    logger.info(f"✓ PASS: {test_name}")
                    self.passed += 1
                else:
                    logger.error(f"✗ FAIL: {test_name} - Expected empty list, got {result}")
                    self.failed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Exception: {e}")
            self.failed += 1
    
    def test_malformed_json_handling(self):
        """Test that malformed JSON responses are handled gracefully"""
        test_name = "Malformed JSON Handling"
        try:
            with patch('requests.get') as mock_get:
                # Simulate malformed JSON response
                mock_response = Mock()
                mock_response.text = '{invalid json'
                mock_response.status_code = 200
                mock_response.json.side_effect = json.JSONDecodeError('msg', 'doc', 0)
                mock_response.raise_for_status.return_value = None
                mock_get.return_value = mock_response
                
                # Should return empty list, not crash
                result = self.engine.search_geo('breast cancer', max_results=10)
                
                if result == []:
                    logger.info(f"✓ PASS: {test_name}")
                    self.passed += 1
                else:
                    logger.error(f"✗ FAIL: {test_name}")
                    self.failed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Exception: {e}")
            self.failed += 1
    
    def test_rate_limiting_detection(self):
        """Test that HTTP 429 rate limiting is detected and triggers retry"""
        test_name = "Rate Limiting Detection (429)"
        try:
            with patch('requests.get') as mock_get:
                # Simulate rate limiting response
                mock_response = Mock()
                mock_response.status_code = 429
                mock_response.raise_for_status.return_value = None
                mock_get.return_value = mock_response
                
                # Should raise RequestException on 429
                try:
                    self.engine.search_geo('breast cancer', max_results=10)
                    logger.error(f"✗ FAIL: {test_name} - Should have raised exception for 429")
                    self.failed += 1
                except requests.RequestException:
                    logger.info(f"✓ PASS: {test_name}")
                    self.passed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Unexpected exception: {e}")
            self.failed += 1
    
    def test_valid_response_parsing(self):
        """Test that valid responses are parsed correctly"""
        test_name = "Valid Response Parsing"
        try:
            with patch('requests.get') as mock_get:
                # Simulate valid GEO response
                valid_response = {
                    'esearchresult': {
                        'idlist': ['GSE12345', 'GSE67890']
                    }
                }
                
                mock_response = Mock()
                mock_response.text = json.dumps(valid_response)
                mock_response.status_code = 200
                mock_response.json.return_value = valid_response
                mock_response.raise_for_status.return_value = None
                mock_get.return_value = mock_response
                
                result = self.engine.search_geo('breast cancer', max_results=10)
                
                if result == ['GSE12345', 'GSE67890']:
                    logger.info(f"✓ PASS: {test_name}")
                    self.passed += 1
                else:
                    logger.error(f"✗ FAIL: {test_name} - Got {result}")
                    self.failed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Exception: {e}")
            self.failed += 1
    
    def test_fetch_metadata_return_statement(self):
        """Test that fetch_dataset_metadata returns metadata dict correctly"""
        test_name = "Fetch Metadata Return Statement"
        try:
            with patch('requests.get') as mock_get:
                # Simulate valid metadata response
                valid_response = {
                    'result': {
                        'GSE12345': {
                            'title': 'Test Dataset',
                            'summary': 'Test summary',
                            'organism': 'Homo sapiens',
                            'n_samples': 10,
                            'platform': 'GPL570'
                        }
                    }
                }
                
                mock_response = Mock()
                mock_response.text = json.dumps(valid_response)
                mock_response.status_code = 200
                mock_response.json.return_value = valid_response
                mock_response.raise_for_status.return_value = None
                mock_get.return_value = mock_response
                
                result = self.engine.fetch_dataset_metadata('GSE12345')
                
                if result and isinstance(result, dict) and 'title' in result:
                    logger.info(f"✓ PASS: {test_name}")
                    self.passed += 1
                else:
                    logger.error(f"✗ FAIL: {test_name} - Got {result}")
                    self.failed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Exception: {e}")
            self.failed += 1
    
    def test_accession_format_validation(self):
        """Test that accession formats are validated correctly"""
        test_name = "Accession Format Validation"
        try:
            # Test valid formats
            valid_accessions = ['GSE123456', 'SRX987654', 'PRJNA123456']
            all_valid = all(validate_accession_format(acc) for acc in valid_accessions)
            
            # Test invalid formats
            invalid_accessions = ['INVALID', 'GSE', 'SRX-123']
            all_invalid = all(not validate_accession_format(acc) for acc in invalid_accessions)
            
            if all_valid and all_invalid:
                logger.info(f"✓ PASS: {test_name}")
                self.passed += 1
            else:
                logger.error(f"✗ FAIL: {test_name}")
                self.failed += 1
        except Exception as e:
            logger.error(f"✗ FAIL: {test_name} - Exception: {e}")
            self.failed += 1
    
    def run_all_tests(self):
        """Run all tests and report results"""
        logger.info("=" * 60)
        logger.info("Running GEO Query Fixes Test Suite")
        logger.info("=" * 60)
        
        self.test_empty_response_handling()
        self.test_malformed_json_handling()
        self.test_rate_limiting_detection()
        self.test_valid_response_parsing()
        self.test_fetch_metadata_return_statement()
        self.test_accession_format_validation()
        
        logger.info("=" * 60)
        logger.info(f"Test Results: {self.passed} PASSED, {self.failed} FAILED")
        logger.info("=" * 60)
        
        return self.failed == 0


def main():
    """Main test runner"""
    tester = TestGEOFixes()
    success = tester.run_all_tests()
    
    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
