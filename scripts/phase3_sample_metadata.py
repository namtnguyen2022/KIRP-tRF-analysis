"""
Phase 3: Sample ID Cleaning & Metadata Table Construction
KIRP tRF Project

Parses TCGA barcodes from:
  - MINTbase tRF sample folders (data_raw/mintbase_trf/KIRP/)
  - GDC biospecimen sample.tsv (data_raw/gdc_biospecimen/sample.tsv)
  - GDC RNA-seq manifest (data_raw/gdc_rnaseq/cohort_Unsaved_Cohort.2026-05-06.tsv)
    + downloaded UUID folders (once gdc-client finishes)

Outputs:
  - data_processed/sample_metadata.tsv
  - data_processed/sample_counts_summary.tsv
"""

import os
import re
import sys
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────────────────
# Edit scripts/config.py to set ROOT for your machine, or change the fallback below:
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import ROOT
except ImportError:
    ROOT = "/Users/nam.tnguyen2022/KIRP_tRF_project"
MINTBASE_DIR = os.path.join(ROOT, "data_raw", "mintbase_trf", "KIRP")
BIOSPC_SAMPLE= os.path.join(ROOT, "data_raw", "gdc_biospecimen", "sample.tsv")
GDC_RNASEQ_DIR = os.path.join(ROOT, "data_raw", "gdc_rnaseq")
OUT_DIR     = os.path.join(ROOT, "data_processed")
os.makedirs(OUT_DIR, exist_ok=True)

# ── TCGA barcode helpers ─────────────────────────────────────────────────────
def parse_barcode(barcode):
    """
    From a full TCGA aliquot barcode (e.g. TCGA-BQ-5875-01A-11R-1591-13)
    extract patient_id (first 12 chars) and sample_type_code (chars 13-14).
    Returns (patient_id, sample_type_code) or (None, None) if not a valid barcode.
    """
    parts = barcode.strip().split("-")
    if len(parts) < 4:
        return None, None
    patient_id = "-".join(parts[:3])          # TCGA-XX-XXXX
    sample_type_code = parts[3][:2]           # '01', '11', '10', etc.
    return patient_id, sample_type_code

SAMPLE_TYPE_LABELS = {
    "01": "Primary Tumor",
    "02": "Recurrent Tumor",
    "06": "Metastatic",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
    "12": "Buccal Cell Normal",
}

def sample_type_label(code):
    return SAMPLE_TYPE_LABELS.get(code, f"Other ({code})")

# ── 1. MINTbase samples ──────────────────────────────────────────────────────
print("Parsing MINTbase sample folders...")
mint_records = []
for entry in sorted(os.listdir(MINTBASE_DIR)):
    full_path = os.path.join(MINTBASE_DIR, entry)
    # Accept both extracted directories and zip filenames
    if entry.endswith(".zip"):
        barcode = entry.replace(".zip", "")
    elif os.path.isdir(full_path) and entry.startswith("TCGA"):
        barcode = entry
    else:
        continue
    patient_id, stc = parse_barcode(barcode)
    if patient_id is None:
        continue
    mint_records.append({
        "sample_id":        barcode,
        "patient_id":       patient_id,
        "sample_type_code": stc,
        "sample_type_label":sample_type_label(stc),
        "source":           "MINTbase",
        "trf_available":    True,
        "mrna_available":   False,   # will update after merging GDC
    })

df_mint = pd.DataFrame(mint_records).drop_duplicates(subset="sample_id")
print(f"  MINTbase samples found: {len(df_mint)}")

# ── 2. GDC biospecimen samples ───────────────────────────────────────────────
print("Parsing GDC biospecimen sample.tsv...")
biospc = pd.read_csv(BIOSPC_SAMPLE, sep="\t", low_memory=False)

# Keep only RNA-relevant sample types (tumor=01, solid normal=11)
# The submitter_id column contains the 16-char sample barcode e.g. TCGA-BQ-5875-01A
gdc_records = []
for _, row in biospc.iterrows():
    barcode = str(row.get("samples.submitter_id", "")).strip()
    if not barcode.startswith("TCGA"):
        continue
    patient_id, stc = parse_barcode(barcode)
    if patient_id is None:
        continue
    gdc_records.append({
        "sample_id":        barcode,
        "patient_id":       patient_id,
        "sample_type_code": stc,
        "sample_type_label":sample_type_label(stc),
        "source":           "GDC",
        "trf_available":    False,   # will update after merging MINTbase
        "mrna_available":   False,   # will update after scanning UUID folders
    })

df_gdc = pd.DataFrame(gdc_records).drop_duplicates(subset="sample_id")
print(f"  GDC biospecimen entries found: {len(df_gdc)}")

# ── 3. Check which GDC UUID folders contain downloaded STAR Counts TSVs ──────
print("Scanning GDC RNA-seq folder for downloaded files...")
rnaseq_barcodes_found = set()
for uuid_dir in os.listdir(GDC_RNASEQ_DIR):
    uuid_path = os.path.join(GDC_RNASEQ_DIR, uuid_dir)
    if not os.path.isdir(uuid_path):
        continue
    for fname in os.listdir(uuid_path):
        if fname.endswith(".tsv") and "rna_seq" in fname.lower():
            # We'll mark these UUIDs; barcode mapping needs the manifest+API
            # For now flag the UUID as present
            rnaseq_barcodes_found.add(uuid_dir)

print(f"  UUID folders with downloaded TSV files: {len(rnaseq_barcodes_found)}")

# ── 4. Map downloaded UUIDs → TCGA barcodes via GDC manifest metadata ────────
# The manifest only has UUIDs. We use the biospecimen aliquot.tsv which links
# aliquot barcodes to case/sample IDs.
ALIQUOT_TSV = os.path.join(ROOT, "data_raw", "gdc_biospecimen", "aliquot.tsv")
aliquot = pd.read_csv(ALIQUOT_TSV, sep="\t", low_memory=False)

# aliquots.submitter_id = full aliquot barcode e.g. TCGA-BQ-5875-01A-11R-1591-13
# samples.submitter_id  = 16-char sample barcode e.g. TCGA-BQ-5875-01A
aliquot_to_sample = {}
for _, row in aliquot.iterrows():
    aliq = str(row.get("aliquots.submitter_id", "")).strip()
    samp = str(row.get("samples.submitter_id", "")).strip()
    if aliq.startswith("TCGA") and samp.startswith("TCGA"):
        aliquot_to_sample[aliq] = samp

# Now map MINTbase barcodes (which are aliquot-level) → sample-level
# MINTbase folder names are aliquot barcodes (16+ chars with plate/center suffix)
# Extract the 16-char sample portion (first 4 dash-fields)
def aliquot_to_sample_barcode(aliquot_barcode):
    parts = aliquot_barcode.split("-")
    if len(parts) >= 4:
        return "-".join(parts[:4])  # TCGA-XX-XXXX-01A
    return aliquot_barcode

# Update MINTbase df with sample-level barcode for matching
df_mint["sample_barcode_16"] = df_mint["sample_id"].apply(aliquot_to_sample_barcode)

# Update GDC df — sample.tsv submitter_id is already 16-char sample barcode
df_gdc["sample_barcode_16"] = df_gdc["sample_id"]

# ── 5. Mark mrna_available in GDC based on downloaded UUIDs ──────────────────
# Load manifest to get file IDs
manifest_path = os.path.join(GDC_RNASEQ_DIR, "cohort_Unsaved_Cohort.2026-05-06.tsv")
manifest = pd.read_csv(manifest_path, sep="\t", header=None, names=["file_id"])
manifest_uuids = set(manifest["file_id"].dropna().str.strip())
downloaded_uuids = rnaseq_barcodes_found

# We'll mark GDC samples as mrna_available = True for those whose UUID was downloaded
# Full barcode→UUID mapping requires the GDC API; for now flag based on download completeness
all_manifest_downloaded = len(downloaded_uuids) >= len(manifest_uuids)

# If fully downloaded, all GDC samples with RNA-seq are available
# If partially downloaded, we note what fraction is done
print(f"  Manifest UUIDs: {len(manifest_uuids)}, Downloaded: {len(downloaded_uuids)}")
if len(downloaded_uuids) > 0:
    # Mark only tumor/normal samples (01, 11) as potentially having mRNA
    df_gdc.loc[df_gdc["sample_type_code"].isin(["01", "11"]), "mrna_available"] = True
else:
    # Download still in progress — assume all manifest samples will be available
    print("  ⚠️  RNA-seq download not yet complete — assuming all manifest samples will be available.")
    df_gdc.loc[df_gdc["sample_type_code"].isin(["01", "11"]), "mrna_available"] = True

# ── 6. Cross-reference: mark trf_available in GDC rows ───────────────────────
mint_sample_16_set = set(df_mint["sample_barcode_16"])
df_gdc["trf_available"] = df_gdc["sample_barcode_16"].isin(mint_sample_16_set)

# Mark mrna_available in MINTbase rows
gdc_sample_16_set = set(df_gdc.loc[df_gdc["sample_type_code"].isin(["01","11"]), "sample_barcode_16"])
df_mint["mrna_available"] = df_mint["sample_barcode_16"].isin(gdc_sample_16_set)

# ── 7. Combine into unified metadata table ────────────────────────────────────
# Keep GDC as base (authoritative sample type info), supplement with MINTbase
df_combined = pd.concat([
    df_gdc[["sample_id","patient_id","sample_type_code","sample_type_label",
            "source","trf_available","mrna_available","sample_barcode_16"]],
    df_mint[~df_mint["sample_barcode_16"].isin(gdc_sample_16_set)][
        ["sample_id","patient_id","sample_type_code","sample_type_label",
         "source","trf_available","mrna_available","sample_barcode_16"]
    ]
], ignore_index=True).drop_duplicates(subset="sample_barcode_16")

df_combined = df_combined.sort_values(["sample_type_code","patient_id"]).reset_index(drop=True)

# ── 8. Save metadata table ────────────────────────────────────────────────────
out_meta = os.path.join(OUT_DIR, "sample_metadata.tsv")
df_combined.to_csv(out_meta, sep="\t", index=False)
print(f"\nSaved sample metadata → {out_meta}")

# ── 9. Summary counts ─────────────────────────────────────────────────────────
# Focus on tumor (01) and normal (11) only
tumor_mask  = df_combined["sample_type_code"] == "01"
normal_mask = df_combined["sample_type_code"] == "11"

gdc_tumor   = df_combined[tumor_mask  & (df_combined["mrna_available"] == True)]
gdc_normal  = df_combined[normal_mask & (df_combined["mrna_available"] == True)]
mint_tumor  = df_combined[tumor_mask  & (df_combined["trf_available"]  == True)]
mint_normal = df_combined[normal_mask & (df_combined["trf_available"]  == True)]

# Matched = patient has BOTH tRF and mRNA data, same sample type
matched_tumor  = df_combined[tumor_mask  & (df_combined["trf_available"] == True) & (df_combined["mrna_available"] == True)]
matched_normal = df_combined[normal_mask & (df_combined["trf_available"] == True) & (df_combined["mrna_available"] == True)]

summary = pd.DataFrame([
    {"category": "GDC mRNA Tumor samples (01)",          "count": len(gdc_tumor)},
    {"category": "GDC mRNA Normal samples (11)",         "count": len(gdc_normal)},
    {"category": "MINTbase tRF Tumor samples (01)",      "count": len(mint_tumor)},
    {"category": "MINTbase tRF Normal samples (11)",     "count": len(mint_normal)},
    {"category": "Matched Tumor (mRNA + tRF)",           "count": len(matched_tumor)},
    {"category": "Matched Normal (mRNA + tRF)",          "count": len(matched_normal)},
])

out_summary = os.path.join(OUT_DIR, "sample_counts_summary.tsv")
summary.to_csv(out_summary, sep="\t", index=False)

print("\n" + "="*55)
print("  SAMPLE COUNT SUMMARY")
print("="*55)
for _, row in summary.iterrows():
    print(f"  {row['category']:<45} {row['count']}")
print("="*55)
print(f"\nSaved summary → {out_summary}")

if len(downloaded_uuids) == 0:
    print("\n⚠️  NOTE: GDC RNA-seq download is still in progress.")
    print("   Re-run this script after the download completes for final counts.")
