#!/usr/bin/env python3
"""
Diagnostic script to test NCBI API responses and identify JSON parsing issues.
"""

import requests
import json
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    import config
    from utils import setup_logger
except ImportError as e:
    print(f"Error importing modules: {e}")
    sys.exit(1)

logger = setup_logger(__name__)

def test_geo_search_api():
    """Test GEO search API with diagnostic output"""
    logger.info("=" * 60)
    logger.info("Testing GEO Search API")
    logger.info("=" * 60)
    
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        'db': 'gds',
        'term': '(breast cancer) AND (RNA-seq) AND (homo sapiens) AND (GSE)',
        'retmax': 5,
        'rettype': 'json',
        'email': config.NCBI_EMAIL
    }
    
    if config.NCBI_API_KEY:
        params['api_key'] = config.NCBI_API_KEY
        logger.info("Using API key for increased rate limits")
    
    try:
        logger.info(f"Request URL: {search_url}")
        logger.info(f"Parameters: {params}")
        
        response = requests.get(search_url, params=params, timeout=config.API_TIMEOUT)
        
        logger.info(f"Status Code: {response.status_code}")
        logger.info(f"Response Headers: {dict(response.headers)}")
        logger.info(f"Response Text Length: {len(response.text)} characters")
        
        if not response.text or not response.text.strip():
            logger.error("❌ EMPTY RESPONSE - This causes 'Expecting value' error!")
            return False
        
        logger.info(f"First 500 chars of response:\n{response.text[:500]}")
        
        try:
            data = response.json()
            logger.info("✓ JSON parsed successfully")
            logger.info(f"Response keys: {list(data.keys())}")
            
            if 'esearchresult' in data:
                logger.info(f"esearchresult keys: {list(data['esearchresult'].keys())}")
                if 'idlist' in data['esearchresult']:
                    ids = data['esearchresult']['idlist']
                    logger.info(f"✓ Found {len(ids)} GEO IDs: {ids[:3]}...")
            
            return True
            
        except json.JSONDecodeError as e:
            logger.error(f"❌ JSON DECODE ERROR: {e}")
            logger.error(f"Response text: {response.text}")
            return False
            
    except requests.RequestException as e:
        logger.error(f"❌ REQUEST ERROR: {e}")
        return False


def test_geo_fetch_api():
    """Test GEO fetch API with diagnostic output"""
    logger.info("=" * 60)
    logger.info("Testing GEO Fetch API")
    logger.info("=" * 60)
    
    fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        'db': 'gds',
        'id': 'GSE12345',  # Test with a common ID
        'rettype': 'json',
        'email': config.NCBI_EMAIL
    }
    
    if config.NCBI_API_KEY:
        params['api_key'] = config.NCBI_API_KEY
    
    try:
        logger.info(f"Request URL: {fetch_url}")
        logger.info(f"Parameters: {params}")
        
        response = requests.get(fetch_url, params=params, timeout=config.API_TIMEOUT)
        
        logger.info(f"Status Code: {response.status_code}")
        logger.info(f"Response Headers: {dict(response.headers)}")
        logger.info(f"Response Text Length: {len(response.text)} characters")
        
        if not response.text or not response.text.strip():
            logger.error("❌ EMPTY RESPONSE - This causes 'Expecting value' error!")
            return False
        
        logger.info(f"First 500 chars of response:\n{response.text[:500]}")
        
        try:
            data = response.json()
            logger.info("✓ JSON parsed successfully")
            logger.info(f"Response keys: {list(data.keys())}")
            
            if 'result' in data:
                logger.info(f"Number of results: {len(data['result'])}")
                if 'GSE12345' in data['result']:
                    logger.info(f"✓ Found GSE12345 in results")
            
            return True
            
        except json.JSONDecodeError as e:
            logger.error(f"❌ JSON DECODE ERROR: {e}")
            logger.error(f"Response text: {response.text}")
            return False
            
    except requests.RequestException as e:
        logger.error(f"❌ REQUEST ERROR: {e}")
        return False


def main():
    """Run all diagnostics"""
    logger.info("NCBI API Diagnostic Test")
    logger.info(f"NCBI Email: {config.NCBI_EMAIL}")
    logger.info(f"API Key Configured: {bool(config.NCBI_API_KEY)}")
    
    results = []
    
    results.append(("GEO Search API", test_geo_search_api()))
    logger.info("")
    
    results.append(("GEO Fetch API", test_geo_fetch_api()))
    
    logger.info("=" * 60)
    logger.info("Diagnostic Summary")
    logger.info("=" * 60)
    
    for test_name, passed in results:
        status = "✓ PASS" if passed else "❌ FAIL"
        logger.info(f"{test_name}: {status}")
    
    all_passed = all(passed for _, passed in results)
    
    if all_passed:
        logger.info("\n✓ All API tests passed! Pipeline should work correctly.")
        return 0
    else:
        logger.error("\n❌ Some API tests failed. See details above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
