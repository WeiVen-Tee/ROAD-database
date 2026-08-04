"""
Map PDB IDs (from mmCIF filenames) to NextGen SIFTS annotations (CATH, SCOP, UniProt).

Reads the NextGen annotation CSVs produced by extract_nextgen_annotations.py and
filters them to the PDB IDs present in a given CIF directory, then splits the unified
xref_segments CSV into database-specific files.

Normally the two steps share the same CIF directory, so the PDB ID sets already match.
The filter here acts as a safety net for cases where the CSVs were produced from a
different or larger CIF directory, or where some downloads failed in the first step.

Usage:
    python map_pdb_nextgen_annotations.py <mmcif_dir> [-v]

Outputs (result/mapping/):
    mapped_uniprot.csv, mapped_cath.csv, mapped_scop.csv,
    mapped_xref_all.csv, mapped_summary.csv
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import RESULT_MAPPING_DIR, NEXTGEN_OUT_DIR, LOG_DIR
from utils import pdb_ids_from_cif_dir, setup_logging

log = None  # set in main()


def main(mmcif_dir: str, verbose: bool = False) -> None:
    """Restrict the NextGen SIFTS CSVs to the PDB IDs present in mmcif_dir and split them.

    Outputs mapped_uniprot.csv (chain-level segments, the SIFTS arithmetic
    input), mapped_cath.csv, mapped_scop.csv, mapped_xref_all.csv, and
    mapped_summary.csv (per-PDB row counts). A PDB missing from the CSVs shows
    zero counts rather than raising.
    """
    global log
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'mapping', 'map_pdb_nextgen_annotations.log'))
    os.makedirs(RESULT_MAPPING_DIR, exist_ok=True)

    pdb_ids = set(pdb_ids_from_cif_dir(mmcif_dir))
    if not pdb_ids:
        raise ValueError(f"No *-assembly1.cif files found in: {mmcif_dir}")
    log.info(f"Loaded {len(pdb_ids):,} PDB IDs from {mmcif_dir}")
    print(f"Mapping {len(pdb_ids):,} PDB IDs to NextGen SIFTS annotations...")

    log.info("Loading annotation CSVs...")
    # Read all columns as strings to avoid pandas inferring wrong dtypes (e.g. PDB IDs
    # like '1234' being misread as integers).
    df_unp  = pd.read_csv(os.path.join(NEXTGEN_OUT_DIR, "pdbx_sifts_unp_segments.csv"),     dtype=str)
    df_xref = pd.read_csv(os.path.join(NEXTGEN_OUT_DIR, "pdbx_sifts_xref_db_segments.csv"), dtype=str)
    log.info(f"  unp_segments: {len(df_unp):,}  xref_segments: {len(df_xref):,}")

    def _filter(df: pd.DataFrame) -> pd.DataFrame:
        """Keep only rows whose pdb_id is in the target set (case-insensitive)."""
        return df[df["pdb_id"].str.upper().isin(pdb_ids)].reset_index(drop=True)

    df_uniprot   = _filter(df_unp)
    df_xref_all  = _filter(df_xref)
    # Split the xref_all DataFrame by database type into dedicated files.
    df_cath = df_xref_all[df_xref_all["xref_db"] == "CATH"].reset_index(drop=True)
    df_scop = df_xref_all[df_xref_all["xref_db"].isin(["SCOP2", "SCOP2B"])].reset_index(drop=True)

    for label, df in [("UniProt", df_uniprot), ("CATH", df_cath), ("SCOP", df_scop),
                      ("xref_all", df_xref_all)]:
        n = df["pdb_id"].str.upper().nunique()
        log.info(f"  {label}: {n:,}/{len(pdb_ids):,} PDB IDs covered, {len(df):,} rows")

    outputs = {
        "mapped_uniprot.csv":   df_uniprot,
        "mapped_cath.csv":      df_cath,
        "mapped_scop.csv":      df_scop,
        "mapped_xref_all.csv":  df_xref_all,
    }
    for fname, df in outputs.items():
        path = os.path.join(RESULT_MAPPING_DIR, fname)
        df.to_csv(path, index=False)
        log.info(f"Saved {fname}: {len(df):,} rows")

    # Per-PDB summary — one row per PDB ID, counts per annotation type.
    summary = pd.DataFrame({"pdb_id": sorted(pdb_ids)})
    for col, df in [("n_uniprot_segs", df_uniprot), ("n_cath_domains", df_cath),
                    ("n_scop_domains", df_scop)]:
        # Group by uppercase pdb_id and count rows; map back via the sorted pdb_id column.
        counts = df.groupby(df["pdb_id"].str.upper()).size()
        summary[col] = summary["pdb_id"].map(counts).fillna(0).astype(int)
    summary_path = os.path.join(RESULT_MAPPING_DIR, "mapped_summary.csv")
    summary.to_csv(summary_path, index=False)
    log.info(f"Saved mapped_summary.csv: {len(summary):,} rows")

    print(f"Done. {len(outputs) + 1} CSV files → {RESULT_MAPPING_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mmcif_dir", help="Directory containing *-assembly1.cif files")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.mmcif_dir, args.verbose)
