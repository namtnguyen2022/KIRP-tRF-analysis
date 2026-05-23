"""
Phase 4b: Build mRNA Raw Count Matrix from GDC STAR-Counts TSV files
KIRP tRF Project

Reads unstranded counts from each UUID folder's .rna_seq.augmented_star_gene_counts.tsv
Maps UUID → TCGA sample barcode via the GDC files-table or biospecimen aliquot.tsv

Produces:
  data_processed/mrna_count_matrix.tsv     (genes x samples, raw unstranded counts)
  data_processed/mrna_sample_metadata.tsv  (sample_id, patient_id, sample_type_code, label)
"""

import os
import sys
import pandas as pd

# Edit scripts/config.py to set ROOT for your machine, or change the fallback below:
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import ROOT
except ImportError:
    ROOT = "/Users/nam.tnguyen2022/KIRP_tRF_project"
GDC_DIR      = os.path.join(ROOT, "data_raw", "gdc_rnaseq",
                             "gdc_download_20260506_091015.911939")
UUID_MAP_TSV = os.path.join(ROOT, "data_processed", "gdc_uuid_barcode_map.tsv")
OUT_DIR      = os.path.join(ROOT, "data_processed")
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLE_TYPE_LABELS = {"01": "Primary_Tumor", "11": "Solid_Tissue_Normal"}

def parse_barcode(barcode):
    parts = str(barcode).strip().split("-")
    if len(parts) < 4:
        return None, None
    return "-".join(parts[:3]), parts[3][:2]

# ── Load GDC API UUID → sample barcode map ────────────────────────────────────
print("Loading UUID → barcode map from GDC API results ...")
uuid_map = pd.read_csv(UUID_MAP_TSV, sep="\t", low_memory=False)
# file_uuid = folder UUID downloaded by gdc-client
uuid_to_barcode = dict(zip(
    uuid_map["file_uuid"].astype(str).str.strip(),
    uuid_map["sample_barcode"].astype(str).str.strip()
))

uuid_dirs = [d for d in os.listdir(GDC_DIR)
             if os.path.isdir(os.path.join(GDC_DIR, d))]

print(f"Found {len(uuid_dirs)} UUID folders.")
print(f"Mapped {sum(1 for u in uuid_dirs if u in uuid_to_barcode)}/{len(uuid_dirs)} via API map.")

# ── Load counts for each UUID ────────────────────────────────────────────────
print("\nLoading STAR counts ...")
all_counts  = {}
sample_meta = []
skipped     = 0
SKIP_ROWS   = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}

for uuid_dir in sorted(uuid_dirs):
    uuid_path = os.path.join(GDC_DIR, uuid_dir)
    barcode   = uuid_to_barcode.get(uuid_dir)

    # Find the TSV file
    tsv_file = None
    for fname in os.listdir(uuid_path):
        if fname.endswith(".tsv"):
            tsv_file = os.path.join(uuid_path, fname)
            break

    if tsv_file is None:
        skipped += 1
        continue

    try:
        df = pd.read_csv(tsv_file, sep="\t", comment="#", low_memory=False)
    except Exception as e:
        print(f"  ⚠️  Error reading {tsv_file}: {e}")
        skipped += 1
        continue

    # Remove summary rows
    df = df[~df["gene_id"].isin(SKIP_ROWS)]
    # Use gene_id as index, unstranded counts
    counts = dict(zip(df["gene_id"].astype(str), df["unstranded"].fillna(0).astype(int)))

    # Use barcode as column name, or UUID if no barcode found
    col_name = barcode if barcode else uuid_dir
    all_counts[col_name] = counts

    patient_id, stc = parse_barcode(col_name) if barcode else (None, None)
    label = SAMPLE_TYPE_LABELS.get(stc, f"Other_{stc}") if stc else "Unknown"
    sample_meta.append({
        "sample_id":         col_name,
        "patient_id":        patient_id,
        "sample_type_code":  stc,
        "sample_type_label": label,
        "uuid":              uuid_dir,
    })

print(f"Loaded {len(all_counts)} samples  ({skipped} skipped).")

# ── Build matrix ──────────────────────────────────────────────────────────────
count_df = pd.DataFrame(all_counts).fillna(0).astype(int)
count_df.index.name = "gene_id"

# Keep only tumor (01) and normal (11)
meta_df      = pd.DataFrame(sample_meta)
keep_samples = meta_df[meta_df["sample_type_code"].isin(["01", "11"])]["sample_id"].tolist()
count_df     = count_df[[c for c in keep_samples if c in count_df.columns]]
meta_df      = meta_df[meta_df["sample_id"].isin(count_df.columns)].reset_index(drop=True)

print(f"\nMatrix shape (genes x samples): {count_df.shape}")
print(f"  Tumor  samples (01): {(meta_df['sample_type_code']=='01').sum()}")
print(f"  Normal samples (11): {(meta_df['sample_type_code']=='11').sum()}")

# ── Save ──────────────────────────────────────────────────────────────────────
matrix_out = os.path.join(OUT_DIR, "mrna_count_matrix.tsv")
meta_out   = os.path.join(OUT_DIR, "mrna_sample_metadata.tsv")
count_df.to_csv(matrix_out, sep="\t")
meta_df.to_csv(meta_out, sep="\t", index=False)

print(f"\nSaved mRNA count matrix   → {matrix_out}")
print(f"Saved mRNA sample metadata → {meta_out}")

total_per_sample = count_df.sum(axis=0)
print(f"\nPer-sample total count range: "
      f"{total_per_sample.min():,} – {total_per_sample.max():,}")
print(f"Total genes: {len(count_df)}")

# Report unmapped samples
unmapped_final = meta_df[meta_df["sample_type_code"].isna()]
if len(unmapped_final):
    print(f"\n⚠️  {len(unmapped_final)} samples could not be mapped to a TCGA barcode.")
    print("   These were excluded from the matrix.")
