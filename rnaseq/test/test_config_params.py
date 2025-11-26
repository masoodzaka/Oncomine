from config import MAX_GEO_RECORDS, MAX_SRA_RECORDS, MAX_ENA_RECORDS

print("=" * 60)
print("Database Query Configuration Parameters")
print("=" * 60)
print(f"MAX_GEO_RECORDS: {MAX_GEO_RECORDS}")
print(f"MAX_SRA_RECORDS: {MAX_SRA_RECORDS}")
print(f"MAX_ENA_RECORDS: {MAX_ENA_RECORDS}")
print("=" * 60)

import os
print("\nEnvironment variable overrides:")
print(f"MAX_GEO_RECORDS env:  {os.getenv('MAX_GEO_RECORDS', 'not set')}")
print(f"MAX_SRA_RECORDS env:  {os.getenv('MAX_SRA_RECORDS', 'not set')}")
print(f"MAX_ENA_RECORDS env:  {os.getenv('MAX_ENA_RECORDS', 'not set')}")

print("\nTo override, set environment variables:")
print("  export MAX_GEO_RECORDS=100")
print("  export MAX_SRA_RECORDS=50")
print("  export MAX_ENA_RECORDS=200")
