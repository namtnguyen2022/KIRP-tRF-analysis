"""
Phase 4a: Build tRF Raw Count Matrix from MINTbase Exclusive-tRF Expression Files
KIRP tRF Project

For each sample, reads the *-exclusive-tRFs.expression.txt file and extracts:
  - Column 1: MINTbase Unique ID  (e.g. tRF-19-6SM83OJX)
  - Column 4: Unnormalized read counts  ← used for edgeR

Produces:
  data_processed/trf_count_matrix.tsv   (tRFs × samples, raw integer counts)
  data_processed/trf_sample_metadata.tsv (sample_id, patient_id, sample_type_code, sample_type_label)
"""

import os
import re
import sys
import pandas as pd

# Edit scripts/config.py to set ROOT for your machine, or change the fallback below:
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import ROOT
except ImportError:
    ROOT = "/Users/nam.tnguyen2022/KIRP_tRF_project"
MINTBASE_DIR = os.path.join(ROOT, "data_raw", "mintbase_trf", "KIRP")
OUT_DIR      = os.path.join(ROOT, "data_processed")
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLE_TYPE_LABELS = {"01": "Primary_Tumor", "11": "Solid_Tissue_Normal"}

def parse_barcode(barcode):
    parts = barcode.strip().split("-")
    if len(parts) < 4:
        return None, None
    return "-".join(parts[:3]), parts[3][:2].zfill(2)

# ── Collect per-sample count vectors ─────────────────────────────────────────
all_counts   = {}   # sample_id → {trf_id: count}
sample_meta  = []

subdirs = sorted([
    d for d in os.listdir(MINTBASE_DIR)
    if os.path.isdir(os.path.join(MINTBASE_DIR, d)) and d.startswith("TCGA")
])

print(f"Found {len(subdirs)} extracted sample directories.")

skipped = 0
for sample_dir in subdirs:
    sample_path = os.path.join(MINTBASE_DIR, sample_dir)
    barcode     = sample_dir  # e.g. TCGA-BQ-5875-01A-11R-1591-13

    # Find the exclusive-tRFs expression file
    expr_file = None
    for fname in os.listdir(sample_path):
        if "exclusive-tRFs.expression.txt" in fname:
            expr_file = os.path.join(sample_path, fname)
            break

    if expr_file is None:
        print(f"  ⚠️  No exclusive expression file for {sample_dir}, skipping.")
        skipped += 1
        continue

    # Parse expression file
    try:
        df = pd.read_csv(expr_file, sep="\t", comment=None, low_memory=False)
    except Exception as e:
        print(f"  ⚠️  Error reading {expr_file}: {e}")
        skipped += 1
        continue

    # Column names may vary slightly — find them robustly
    df.columns = df.columns.str.strip()
    id_col    = df.columns[0]   # MINTbase Unique ID
    count_col = df.columns[3]   # Unnormalized read counts

    counts = dict(zip(df[id_col].astype(str).str.strip(),
                      pd.to_numeric(df[count_col], errors="coerce").fillna(0).astype(int)))

    all_counts[barcode] = counts

    patient_id, stc = parse_barcode(barcode)
    label = SAMPLE_TYPE_LABELS.get(stc, f"Other_{stc}")
    sample_meta.append({
        "sample_id":         barcode,
        "patient_id":        patient_id,
        "sample_type_code":  stc,
        "sample_type_label": label,
    })

print(f"Loaded {len(all_counts)} samples  ({skipped} skipped).")

# ── Build matrix ──────────────────────────────────────────────────────────────
count_df = pd.DataFrame(all_counts).fillna(0).astype(int)
count_df.index.name = "trf_id"

# Keep only tumor (01) and normal (11) samples
meta_df = pd.DataFrame(sample_meta)
keep_samples = meta_df[meta_df["sample_type_code"].isin(["01","11"])]["sample_id"].tolist()
count_df = count_df[[c for c in keep_samples if c in count_df.columns]]
meta_df  = meta_df[meta_df["sample_id"].isin(count_df.columns)].reset_index(drop=True)

print(f"\nMatrix shape (tRFs × samples): {count_df.shape}")
print(f"  Tumor samples  (01): {(meta_df['sample_type_code']=='01').sum()}")
print(f"  Normal samples (11): {(meta_df['sample_type_code']=='11').sum()}")

# ── Save ──────────────────────────────────────────────────────────────────────
matrix_out = os.path.join(OUT_DIR, "trf_count_matrix.tsv")
meta_out   = os.path.join(OUT_DIR, "trf_sample_metadata.tsv")

count_df.to_csv(matrix_out, sep="\t")
meta_df.to_csv(meta_out, sep="\t", index=False)

print(f"\nSaved count matrix  → {matrix_out}")
print(f"Saved sample metadata → {meta_out}")

# ── Quick QC summary ─────────────────────────────────────────────────────────
total_counts_per_sample = count_df.sum(axis=0)
print(f"\nPer-sample total count range: "
      f"{total_counts_per_sample.min():,} – {total_counts_per_sample.max():,}")
print(f"Total unique tRF IDs: {len(count_df)}")
print(f"tRFs with all-zero counts: {(count_df.sum(axis=1) == 0).sum()}")
