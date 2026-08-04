# The ROAD (Repository for Orthosteric and Allosteric Design) database and analysis pipeline 

- Repository of allosteric and orthosteric (i.e. functional) sites and ligands. 

- Analysis of multiple characteristics of allosteric and orthosteric/functional sites and ligands.

## Site categories

Binding sites are classified by `min_dist`, the distance from the site to the
nearest UniProt-annotated Active/DNA-binding residue:

| Category    | Criterion               |
|--------------|--------------------------|
| Orthosteric  | `min_dist ≤ 5 Å`         |
| Allosteric   | `min_dist > 5 Å`         |

## Repository layout

```
collection/   Data collection pipeline: PDB/UniProt/SIFTS mapping, binding-site
              extraction, non-redundant clustering, metadata enrichment.
              Outputs are tabular data, and PDB files (.cif) 
analysis/     Analysis pipeline: characteristics of ligands, binding sites, and their interactions.
```

`collection` produces the non-redundant, annotated binding-site table, a copy of
which lives in `analysis/input_dataset/` for the analysis stages to read.
`analysis/src/config.py` also imports shared code and constants from
`collection/src/` via a relative path (`../collection`), so the two directories
must stay siblings under one parent as they are here.

## Setup

1. Create the conda environment:
   ```bash
   conda env create -f environment.yml   # collection + analysis
   conda activate ROAD-database
   ```
2. Point both projects at your local structure-file cache (files named
   `{PDB_ID}-assembly1.cif`). Neither repo hardcodes a machine-specific path —
   either set the `STRUCTURE_CIF_DIR` environment variable, or create a
   config from the shipped template and fill in `CIF_DIR`:
   ```bash
   cp collection/src/config_local.py.example collection/src/config_local.py
   # then edit collection/src/config_local.py and set CIF_DIR
   ```
   `analysis` reads `CIF_DIR` from `collection`'s config automatically, so you
   normally only need to set it once; `analysis/src/config_local.py.example`
   exists only as a fallback for the rare case where you run `analysis`
   standalone without a `collection/src/config_local.py` in place.
3. Hand `collection`'s output to `analysis` after running the collection
   pipeline — `analysis` reads it from a fixed location and does not reach into
   `collection/result/` for it:
   ```bash
   cp collection/result/nonredundant/nonredundant_sites_annotated_enriched.csv \
      analysis/input_dataset/
   ```
4. Obtain and place the third-party datasets as follows:
   - The PDB Chemical Component Dictionary (~500 MB), which both pipelines
     read independently from their own paths:
     `collection/data/Components-rel-alt.cif` and
     `analysis/ccd/Components-rel-alt.cif`. `collection/README.md` has the
     download link and the decompress/rename steps; the same file serves both.

## Pipeline stages

```
── collection ────────────────────────────────────────────────────────────────
nextgen_annotations, uniprot_mapping,
chain_unp_map, uniprot_names            SIFTS/UniProt mapping            → result/mapping/
ligand_blacklist                        build_badlig                     → result/ligands/
site_extraction, site_filter            binding-site extraction/filter   → result/binding_sites/
uniprot_features                        fetch_uniprot_features           → result/functional/
functional_site_mapping                 map_functional_sites             → result/functional/
distances                               calculate_distances              → result/distances/
nonredundant_build_*                    non-redundant clustering         → result/nonredundant/
protein_names                           Protein Identity (better_protein_names) → result/protein_names/
nonredundant_enrich_*                   distance/structural/functional enrich → result/nonredundant/
nonredundant_pdb_copy                   copy_nonredundant_pdb            → result/nonredundant/PDB{,_allosteric}.tar.gz

── analysis ──────────────────────────────────────────────────────────────────
Preparation             compute_exclusive_categories (inline)
Ligand Profiling        enrich_ligand_smiles.py → analyze_ligand_descriptors.py
Site Composition        analyze_aa_composition.py
Site Geometry           analyze_site.py            (size + SASA + curvature)
Interaction Contacts    analyze_interactions.py     (needs pdbe-arpeggio; see environment.yml)
Fragment Analysis       analyze_recap_fragments.py
Feature Build           build_site_ligand_features (inline; runs last, merging
                          every feature stage above into one table)
Protein Clustering      cluster_protein_similarity.py (needs foldseek; reads Feature Build)
```

Run collection first (`python collection/run_collection_pipeline.py`), then
analysis (`python analysis/run_analysis_pipeline.py`). See each subfolder's
README for stage-by-stage detail and CLI flags.

## Data

`ROAD_DATA/` holds the two result tables, one row per non-redundant
ligand-binding site:

| File | Shape | Contents |
|---|---|---|
| `nonredundant_sites_annotated_enriched.csv` | 27,164 × 20 | `collection`'s final output: the site residues, its orthosteric/allosteric category, distances to the functional site, and protein annotation. |
| `site_ligand_features.csv` | 27,164 × 91 | The same sites joined to `analysis`'s features: ligand descriptors, site geometry and solvent accessibility, amino-acid composition, and protein–ligand contact fractions. |

Every column is described in the companion `*_fields.txt` file. The two tables share
the same `ID` values, so they join one-to-one.

## Citation

See `CITATION.txt`.

## License

Copyright (C) 2026 Wei-Ven Tee

GNU General Public License v3.0 — see `LICENSE`.
