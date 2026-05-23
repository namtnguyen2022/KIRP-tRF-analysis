"""
Phase 7a: Extract tRF sequences and fetch mRNA 3'UTR sequences
KIRP tRF Project

For each top anticorrelated tRF-gene pair:
  1. Extract tRF sequence from MINTbase expression files
  2. Fetch 3'UTR sequence from Ensembl REST API (canonical transcript)

Outputs:
  data_processed/trf_sequences.tsv
  data_processed/mrna_3utr_sequences.tsv
"""

import os
import re
import sys
import time
import requests
import pandas as pd

# Edit scripts/config.py to set ROOT for your machine, or change the fallback below:
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    from config import ROOT
except ImportError:
    ROOT = "/Users/nam.tnguyen2022/KIRP_tRF_project"
MINTBASE_DIR = os.path.join(ROOT, "data_raw", "mintbase_trf", "KIRP")
ANTI_FILE    = os.path.join(ROOT, "results", "tables", "anticorrelation_prioritized.tsv")
OUT_DIR      = os.path.join(ROOT, "data_processed")

ENSEMBL_API  = "https://rest.ensembl.org"

# ── 1. Load top prioritized pairs ────────────────────────────────────────────
anti = pd.read_csv(ANTI_FILE, sep="\t")
# Take top 20 pairs (covers top tRFs and all 3 genes)
top_pairs = anti.head(20).copy()

top_trfs   = top_pairs["trf"].unique().tolist()
top_genes  = top_pairs["gene"].unique().tolist()

print(f"Top tRFs to get sequences for: {len(top_trfs)}")
print(f"Genes to fetch 3'UTRs for:     {top_genes}")

# ── 2. Extract tRF sequences from MINTbase ───────────────────────────────────
print("\nExtracting tRF sequences from MINTbase files ...")
trf_seq_map = {}   # trf_id -> sequence

# Only need to find each tRF once — scan sample files until all found
needed = set(top_trfs)
subdirs = sorted([
    d for d in os.listdir(MINTBASE_DIR)
    if os.path.isdir(os.path.join(MINTBASE_DIR, d)) and d.startswith("TCGA")
])

for sample_dir in subdirs:
    if not needed:
        break
    sample_path = os.path.join(MINTBASE_DIR, sample_dir)
    for fname in os.listdir(sample_path):
        if "exclusive-tRFs.expression.txt" not in fname:
            continue
        fpath = os.path.join(sample_path, fname)
        try:
            df = pd.read_csv(fpath, sep="\t", low_memory=False)
            df.columns = df.columns.str.strip()
            id_col  = df.columns[0]
            seq_col = df.columns[1]   # tRF sequence
            for _, row in df.iterrows():
                trf_id = str(row[id_col]).strip()
                if trf_id in needed:
                    trf_seq_map[trf_id] = str(row[seq_col]).strip().upper().replace("T", "U")
                    needed.discard(trf_id)
        except Exception:
            pass
        break  # one file per sample dir is enough

print(f"  Sequences found: {len(trf_seq_map)} / {len(top_trfs)}")
if needed:
    print(f"  ⚠️  Not found: {needed}")

trf_df = pd.DataFrame([
    {"trf_id": t, "sequence": trf_seq_map.get(t, ""), "length": len(trf_seq_map.get(t, ""))}
    for t in top_trfs
])
trf_df.to_csv(os.path.join(OUT_DIR, "trf_sequences.tsv"), sep="\t", index=False)
print(f"  Saved: data_processed/trf_sequences.tsv")

# ── 3. Fetch 3'UTR sequences from Ensembl REST API ───────────────────────────
print("\nFetching 3'UTR sequences from Ensembl ...")

GENE_SYMBOLS = {
    "MET":    "ENSG00000105976",
    "CDKN2A": "ENSG00000147889",
    "TERT":   "ENSG00000164362",
}

def get_canonical_transcript(gene_id):
    """Get canonical transcript ID for an Ensembl gene."""
    url = f"{ENSEMBL_API}/lookup/id/{gene_id}?expand=1&content-type=application/json"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()
    # canonical_transcript field
    canonical = data.get("canonical_transcript")
    if canonical:
        return canonical.split(".")[0]
    # fallback: longest transcript
    transcripts = data.get("Transcript", [])
    if transcripts:
        longest = max(transcripts, key=lambda t: t.get("length", 0))
        return longest["id"]
    return None

def get_3utr_sequence(transcript_id):
    """Fetch the 3'UTR sequence for a transcript."""
    url = f"{ENSEMBL_API}/sequence/id/{transcript_id}?content-type=application/json&type=three_prime_utr"
    r = requests.get(url, timeout=30)
    if r.status_code == 400:
        # No 3'UTR annotated
        return None, None
    r.raise_for_status()
    data = r.json()
    seq = data.get("seq", "")
    desc = data.get("desc", "")
    return seq.upper(), desc

utr_records = []
for gene_name, gene_ensembl_id in GENE_SYMBOLS.items():
    if gene_name not in top_genes:
        continue
    print(f"  Fetching {gene_name} ({gene_ensembl_id}) ...")
    try:
        tx_id = get_canonical_transcript(gene_ensembl_id)
        print(f"    Canonical transcript: {tx_id}")
        time.sleep(0.3)

        utr_seq, desc = get_3utr_sequence(tx_id)
        if utr_seq:
            print(f"    3'UTR length: {len(utr_seq)} nt")
            utr_records.append({
                "gene":          gene_name,
                "ensembl_gene":  gene_ensembl_id,
                "transcript_id": tx_id,
                "utr3_length":   len(utr_seq),
                "utr3_sequence": utr_seq,
            })
        else:
            print(f"    ⚠️  No 3'UTR annotated, fetching full CDS+UTR ...")
            # fallback: fetch full mRNA and trim
            url2 = f"{ENSEMBL_API}/sequence/id/{tx_id}?content-type=application/json&type=cdna"
            r2 = requests.get(url2, timeout=30)
            r2.raise_for_status()
            cdna = r2.json().get("seq", "").upper()
            utr_records.append({
                "gene":          gene_name,
                "ensembl_gene":  gene_ensembl_id,
                "transcript_id": tx_id,
                "utr3_length":   len(cdna),
                "utr3_sequence": cdna,  # full cDNA fallback
            })
        time.sleep(0.5)
    except Exception as e:
        print(f"    ⚠️  Error: {e}")

utr_df = pd.DataFrame(utr_records)
utr_df.to_csv(os.path.join(OUT_DIR, "mrna_3utr_sequences.tsv"), sep="\t", index=False)
print(f"\nSaved: data_processed/mrna_3utr_sequences.tsv")
print(utr_df[["gene", "transcript_id", "utr3_length"]])
