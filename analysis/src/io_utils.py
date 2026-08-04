"""
io_utils.py
===========
Shared CSV-loading and DataFrame-prep helpers.

Takes config values as explicit parameters rather than importing config, so
there is no import-order coupling to config.py.
"""

import logging
import os

import pandas as pd


def load_enriched_csv(enriched_csv_path: str,
                       dropna_subset: list | None = None) -> pd.DataFrame:
    """Load nonredundant_sites_annotated_enriched.csv from input_dataset/."""
    df = pd.read_csv(enriched_csv_path)
    if dropna_subset:
        df = df.dropna(subset=dropna_subset)
    return df


def assign_site_category_col(df: pd.DataFrame, assign_fn) -> pd.DataFrame:
    """Populate site_category from min_dist, keeping any value already present.

    collection's enriched CSV already carries site_category; recomputing it here
    only fills gaps, so the two projects cannot disagree on the classification.
    """
    derived = df['min_dist'].apply(assign_fn)
    if 'site_category' in df.columns:
        df['site_category'] = df['site_category'].where(df['site_category'].notna(), derived)
    else:
        df['site_category'] = derived
    return df


def extract_ccd_code_rsplit(df: pd.DataFrame, col: str = 'ligand_code') -> pd.DataFrame:
    df['ccd_code'] = df[col].str.rsplit(':', n=1).str[-1]
    return df


def output_is_cached(path: str, force: bool) -> bool:
    return os.path.exists(path) and not force


def save_csv(df: pd.DataFrame, out_path: str, logger: logging.Logger) -> None:
    """to_csv plus the standard 'Saved ->' print and row-count log line."""
    df.to_csv(out_path, index=False)
    print(f'  Saved → {out_path}')
    logger.info('Wrote %d rows to %s', len(df), out_path)
