"""
Phase 7b: Parse RNA22 output files and produce a clean binding table.
Input:  data_raw/data_raw_rna22_*_results - Sheet1.tsv
Output: data_processed/rna22_binding_results.tsv
        results/tables/rna22_binding_significant.tsv
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

DATA_RAW  = os.path.join(ROOT, "data_raw")
DATA_PROC = os.path.join(ROOT, "data_processed")
RESULTS   = os.path.join(ROOT, "results", "tables")

genes = ["MET", "TERT", "CDKN2A"]

all_rows = []

for gene in genes:
    fname = f"data_raw_rna22_{gene}_results - Sheet1.tsv"
    fpath = os.path.join(DATA_RAW, fname)
    if not os.path.exists(fpath):
        print(f"WARNING: {fname} not found, skipping.")
        continue

    df = pd.read_csv(fpath, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    # Rename columns to standard names
    df = df.rename(columns={
        "miR Name": "tRF_id",
        "transcript name": "transcript",
        "leftmost position of predicted target site": "target_position",
        "folding energy (in -Kcal/mol)": "folding_energy_kcal",
        "heteroduplex": "heteroduplex",
        "p value": "p_value"
    })

    df["gene"] = gene
    df["folding_energy_kcal"] = pd.to_numeric(df["folding_energy_kcal"], errors="coerce")
    df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")

    # Clean up tRF ID: replace _ with - for MINTbase style
    df["tRF_id"] = df["tRF_id"].str.replace("_", "-", n=2)  # only first 2 underscores (type-length-id)

    all_rows.append(df)
    print(f"{gene}: {len(df)} binding sites loaded")

if not all_rows:
    print("No data loaded. Exiting.")
    exit(1)

full = pd.concat(all_rows, ignore_index=True)

# Keep relevant columns
full = full[["tRF_id", "gene", "transcript", "target_position",
             "folding_energy_kcal", "p_value", "heteroduplex"]]

# Sort by p_value then folding energy
full = full.sort_values(["gene", "tRF_id", "p_value", "folding_energy_kcal"])

# Save full table
out_full = os.path.join(DATA_PROC, "rna22_binding_results.tsv")
full.to_csv(out_full, sep="\t", index=False)
print(f"\nFull binding table: {len(full)} rows → {out_full}")

# Significant: p < 0.05
sig = full[full["p_value"] < 0.05].copy()
out_sig = os.path.join(RESULTS, "rna22_binding_significant.tsv")
sig.to_csv(out_sig, sep="\t", index=False)
print(f"Significant (p<0.05): {len(sig)} rows → {out_sig}")

# Summary per tRF-gene pair: best hit (lowest p_value)
summary = (sig.sort_values("p_value")
             .groupby(["tRF_id", "gene"])
             .first()
             .reset_index()
             [["tRF_id", "gene", "target_position", "folding_energy_kcal", "p_value", "heteroduplex"]])
summary = summary.sort_values("p_value")

out_summary = os.path.join(RESULTS, "rna22_binding_summary.tsv")
summary.to_csv(out_summary, sep="\t", index=False)
print(f"Summary (best hit per tRF-gene pair): {len(summary)} rows → {out_summary}")

print("\nTop 20 binding pairs (by p-value):")
print(summary.head(20).to_string(index=False))
