"""
Phase 4b-helper: Query GDC API to map file UUIDs -> TCGA sample barcodes
Saves: data_processed/gdc_uuid_barcode_map.tsv
"""

import os
import sys
import json
import time
import requests
import pandas as pd

# Edit scripts/config.py to set ROOT for your machine, or change the fallback below:
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import ROOT
except ImportError:
    ROOT = "/Users/nam.tnguyen2022/KIRP_tRF_project"
GDC_DIR = os.path.join(ROOT, "data_raw", "gdc_rnaseq",
                        "gdc_download_20260506_091015.911939")
OUT_DIR = os.path.join(ROOT, "data_processed")
OUT_MAP = os.path.join(OUT_DIR, "gdc_uuid_barcode_map.tsv")

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"

# Get all folder UUIDs (these are the GDC case/file bundle UUIDs)
uuid_list = [d for d in os.listdir(GDC_DIR)
             if os.path.isdir(os.path.join(GDC_DIR, d))]
print(f"Querying GDC API for {len(uuid_list)} UUIDs ...")

records = []
batch_size = 50  # GDC API max per request

for i in range(0, len(uuid_list), batch_size):
    batch = uuid_list[i:i+batch_size]
    payload = {
        "filters": {
            "op": "in",
            "content": {
                "field": "file_id",
                "value": batch
            }
        },
        "fields": "file_id,cases.submitter_id,cases.samples.submitter_id,"
                  "cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id",
        "format": "json",
        "size": str(batch_size)
    }

    try:
        resp = requests.post(GDC_FILES_ENDPOINT, json=payload, timeout=30)
        resp.raise_for_status()
        hits = resp.json().get("data", {}).get("hits", [])
        for hit in hits:
            file_id = hit.get("file_id", "")
            for case in hit.get("cases", []):
                case_barcode = case.get("submitter_id", "")
                for sample in case.get("samples", []):
                    sample_barcode = sample.get("submitter_id", "")
                    sample_type    = sample.get("sample_type", "")
                    for portion in sample.get("portions", []):
                        for analyte in portion.get("analytes", []):
                            for aliquot in analyte.get("aliquots", []):
                                aliquot_barcode = aliquot.get("submitter_id", "")
                                records.append({
                                    "file_uuid":      file_id,
                                    "case_barcode":   case_barcode,
                                    "sample_barcode": sample_barcode,
                                    "sample_type":    sample_type,
                                    "aliquot_barcode":aliquot_barcode,
                                })
        print(f"  Batch {i//batch_size + 1}: {len(hits)} hits")
        time.sleep(0.3)  # be polite to the API
    except Exception as e:
        print(f"  ⚠️  Batch {i//batch_size + 1} failed: {e}")

df = pd.DataFrame(records)
print(f"\nTotal records retrieved: {len(df)}")
print(df.head())

# Deduplicate: one row per file_uuid (pick first sample barcode)
df_dedup = df.drop_duplicates(subset="file_uuid")[
    ["file_uuid", "case_barcode", "sample_barcode", "sample_type", "aliquot_barcode"]
]

# Parse sample type code
def get_stc(barcode):
    parts = str(barcode).split("-")
    return parts[3][:2] if len(parts) >= 4 else None

df_dedup = df_dedup.copy()
df_dedup["sample_type_code"] = df_dedup["sample_barcode"].apply(get_stc)

df_dedup.to_csv(OUT_MAP, sep="\t", index=False)
print(f"\nSaved UUID→barcode map → {OUT_MAP}")
print(f"  Tumor  (01): {(df_dedup['sample_type_code']=='01').sum()}")
print(f"  Normal (11): {(df_dedup['sample_type_code']=='11').sum()}")
