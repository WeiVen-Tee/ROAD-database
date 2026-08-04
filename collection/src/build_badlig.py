"""
Build excluded-ligand lists from the PDB Chemical Component Dictionary (CCD) and BioLiP.

Sources:
  1. CCD (Components-rel-alt.cif): sugars, MW outliers
  2. BioLiP crystal-additive list
  3. Ad-hoc hardcoded set

Amino acids (standard or modified, CCD type PEPTIDE LINKING) are intentionally
not blacklisted here — see _get_bad() for why.

Outputs (result/ligands/):
  excluded_lig.p    — set[str]  all excluded-ligand CCD codes
  excluded_lig.txt  — same set, sorted, one code per line (human-readable)
  ions.p            — set[str]  ion/inorganic CCD codes (pdbx_type == HETAI)
  ions.txt          — same set, sorted, one code per line (human-readable)

Usage:
    python build_badlig.py [-v]
"""

import argparse
import os
import re
import sys
from itertools import islice

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import CCD_PATH, BIOLIP_PATH, RESULT_LIGANDS_DIR, BADLIG_ADHOC, MOL_WT_MIN, MOL_WT_MAX, LOG_DIR
from utils import save_pickle, setup_logging

log = None  # set in main()

_QUOTED = re.compile(r'^"(.*)"$')


def _unquote(val: str) -> str:
    """Strip surrounding double-quotes from a CIF string value if present.

    CIF files sometimes wrap values in double-quotes (e.g. "Amino acid").
    This function strips those quotes so the caller gets the raw string.
    """
    m = _QUOTED.match(val)
    return m.group(1) if m else val


def _parse_ccd(cif_path: str) -> list[dict]:
    """Parse Components-rel-alt.cif into [{id, type, pdbx_type, formula_weight}, ...].

    Hand-written rather than gemmi's CIF reader: the CCD is ~500 MB, and only
    these four fields are needed.
    """
    WANTED = {
        "_chem_comp.id", "_chem_comp.type",
        "_chem_comp.pdbx_type", "_chem_comp.formula_weight",
    }
    components: list[dict] = []
    current: dict = {}
    with open(cif_path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("data_"):
                # Each 'data_XXX' block is a new component; save the previous one.
                if current.get("id"):
                    components.append(current)
                current = {}
                continue
            if not line or line.startswith("#") or line.startswith("loop_"):
                continue
            parts = line.split(None, 1)
            if len(parts) == 2 and parts[0] in WANTED:
                # Strip the '_chem_comp.' prefix to get a short key like 'id' or 'type'.
                current[parts[0][len("_chem_comp."):]] = _unquote(parts[1].strip())
    if current.get("id"):
        components.append(current)
    return components


def _get_ions(components: list[dict]) -> set[str]:
    """CCD codes with pdbx_type == 'HETAI': metal ions and simple inorganics.

    Tracked separately from excluded ligands so they can annotate nearby ions.
    """
    return {c["id"] for c in components if c.get("pdbx_type", "").upper() == "HETAI"}


def _get_bad(components: list[dict]) -> set[str]:
    """CCD codes to exclude as non-drug-like: saccharides and MW outliers.

    Components outside [MOL_WT_MIN, MOL_WT_MAX] are buffer fragments or
    polymers rather than small-molecule ligands.

    Peptide-linking residues are skipped rather than filtered: whether such a
    code is a protein residue or a free ligand is per-structure context this
    function cannot see, and most standard amino acids fall below MOL_WT_MIN.
    That call is made at parse time from each residue's polymer membership.
    """
    bad: set[str] = set()
    for comp in components:
        cid       = comp.get("id", "")
        ctype     = comp.get("type", "").upper()
        pdbx_type = comp.get("pdbx_type", "").upper()
        if pdbx_type == "ATOMP" and "PEPTIDE LINKING" in ctype:
            continue
        elif pdbx_type == "ATOMS" and "SACCHARIDE" in ctype:
            bad.add(cid)
        else:
            try:
                mw = float(comp.get("formula_weight", "0"))
            except ValueError:
                mw = 0.0
            if mw < MOL_WT_MIN or mw > MOL_WT_MAX:
                bad.add(cid)
    return bad


def _get_biolip(biolip_path: str) -> set[str]:
    """CCD codes from the BioLiP crystal-additive list (PEGs, detergents, cryoprotectants).

    The file has two header lines before the data rows.
    """
    codes: set[str] = set()
    with open(biolip_path, encoding="utf-8") as fh:
        # islice tolerates a file with fewer than two header lines.
        for _ in islice(fh, 2):
            pass
        for line in fh:
            parts = line.split()
            if parts:
                codes.add(parts[0])
    return codes


def main(verbose: bool = False) -> None:
    """Write excluded_lig.p/.txt and ions.p/.txt.

    Excluded ligands merge three sources: CCD-derived (saccharides, MW
    outliers), the BioLiP crystal-additive list, and BADLIG_ADHOC. Ions stay a
    separate set so later steps can annotate them without treating them as
    binding-site ligands.
    """
    global log
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'build_badlig.log'))
    os.makedirs(RESULT_LIGANDS_DIR, exist_ok=True)

    log.info(f"Parsing CCD: {CCD_PATH}")
    print(f"Parsing CCD: {CCD_PATH} ...")
    components = _parse_ccd(CCD_PATH)
    log.info(f"{len(components):,} components loaded")

    ions    = _get_ions(components)
    bad_ccd = _get_bad(components)
    biolip  = _get_biolip(BIOLIP_PATH)
    # Union of all bad (i.e. exclude) sources; remove ions so they remain available for ion annotation.
    bad_lig = (bad_ccd | biolip | BADLIG_ADHOC) - ions

    log.info(f"Ions: {len(ions)}  bad_ccd: {len(bad_ccd)}  biolip: {len(biolip)}  adhoc: {len(BADLIG_ADHOC)}")
    log.info(f"Total excluded_lig: {len(bad_lig)}")
    print(f"Total excluded_lig: {len(bad_lig):,}  |  Ions: {len(ions)}")

    save_pickle(bad_lig, os.path.join(RESULT_LIGANDS_DIR, 'excluded_lig.p'))
    save_pickle(ions,    os.path.join(RESULT_LIGANDS_DIR, 'ions.p'))
    with open(os.path.join(RESULT_LIGANDS_DIR, 'excluded_lig.txt'), 'w') as f:
        f.write('\n'.join(sorted(bad_lig)) + '\n')
    with open(os.path.join(RESULT_LIGANDS_DIR, 'ions.txt'), 'w') as f:
        f.write('\n'.join(sorted(ions)) + '\n')
    log.info(f"Saved excluded_lig.p/.txt and ions.p/.txt → {RESULT_LIGANDS_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.verbose)
