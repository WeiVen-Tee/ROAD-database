# analysis

Analyses the amino-acid composition, physicochemical properties, geometry,
and ligand-interaction profile of orthosteric vs allosteric binding sites,
using the non-redundant annotated site table produced by `../collection`.

## Setup

This project requires a sibling `../collection` checkout: it puts
`collection/src/` on `sys.path` and imports `config.py` and `utils.py` from it,
and reads `collection/result/` for `pdb_cif_cache.p` and the SCOP2/CATH mappings.

This is already set up and needs no action if you run the pipeline as-is; it only
matters if you modify collection. Both projects have a module named `config`, so
when collection's `utils.py` runs here it does `from config import ...` against
*this* `src/config.py`, not its own. So `src/config.py` must define every constant
collection's `utils.py` imports from config — if you add such an import to
`collection/src/utils.py`, add the constant here too, or the import breaks.

Copy the primary collection output into `analysis/input_dataset/` (create the
folder if it does not exist):

```bash
mkdir -p input_dataset
cp ../collection/result/nonredundant/nonredundant_sites_annotated_enriched.csv input_dataset/
```

CIF directory and conda environment: see the root `README.md`.

One third-party reference dataset is used but not committed to this repo —
download and place it at:
- `ccd/Components-rel-alt.cif` — PDB Chemical Component Dictionary

## Running

```bash
python run_analysis_pipeline.py                          # full run, all stages
python run_analysis_pipeline.py --verbose                 # full run, INFO logging
python run_analysis_pipeline.py --only smiles descriptors # run specific stages
python run_analysis_pipeline.py --only smiles descriptors aa_composition site recap_fragments build cluster
                                                           # every stage except interactions
                                                           # (e.g. no pdbe-arpeggio available)
python run_analysis_pipeline.py --only site --cif-dir PATH # Site Geometry only
                                                           # (for a remote run with real CIF files)
```

`--freesasa-bin`, `--foldseek-bin` and `--threads` are also accepted, for
environments where those executables are not on `PATH` under their usual names
or where the worker count needs capping.

`conda activate ROAD-database` — see the root `README.md` for environment setup.

Stages are cached — a stage is skipped only when **all** of its outputs exist,
so a run interrupted part-way through a stage re-runs it. Pass `--force` to
recompute regardless.

`--only` pulls in each named stage's inputs automatically, so `--only cluster`
also runs `build`. The one exception is `build` itself: all of its merges are
best-effort, so it never re-runs the feature stages — it merges whatever is on
disk and reports which expected columns were missing.

The default full run includes `interactions` (`src/analyze_interactions.py`),
which requires the `pdbe-arpeggio` binary (installed via the single `environment.yml`
env — see root `README.md`) and CIF files. If `pdbe-arpeggio` is unavailable,
use `--only` to list every stage except `interactions`.

**Caveats for `cluster`:**

- Cluster IDs are **file-local**. The same `cluster_struc` integer in
  `clusters_all.csv` and `clusters_allosteric.csv` does not mean the same cluster.
- `cluster_struc` is Foldseek-only (structure), and is
  all-singleton unless `--cif-dir` points at a populated directory of
  `{PDB_ID}-assembly1.cif` files.
- Foldseek parameters used here : `--cluster-mode 1`
  (connected-component) with `-c 0.8` / `--cov-mode 0` / `-e 0.01`.
  These are tunable only on the stage script itself, not through the orchestrator.
- `scop2_sf_ids` / `cath_topology_ids` are annotation-only and never gate
  clustering; name caches live in `data/scop2_sf_names.csv` and
  `data/cath-names.txt`.

## Pipeline stages

| Stage               | Script                             | Output                                   | Depends on                             |
|---------------------|-------------------------------------|-------------------------------------------|-----------------------------------------|
| `exclusive_categories` | `run_analysis_pipeline.py` (inline) | `data/ligand_exclusive_categories.csv`, `data/ligand_site_counts.csv` | —                    |
| `smiles`            | `src/enrich_ligand_smiles.py`      | `data/ligand_smiles.csv`                 | —                                       |
| `descriptors`       | `src/analyze_ligand_descriptors.py` | `data/ligand_descriptors.csv`, `data/plot_ligand_{descriptors,npr}.csv` | `smiles`        |
| `aa_composition`    | `src/analyze_aa_composition.py`    | `data/site_aa_composition.csv`, `data/aa_mean_fractions.csv` | `exclusive_categories`  |
| `site`              | `src/analyze_site.py`              | `data/sasa_results.csv`, `data/curvature_results.csv` | —                          |
| `interactions`      | `src/analyze_interactions.py`      | `data/interaction_contacts{,_shared}.csv` | `descriptors`, `exclusive_categories`  |
| `recap_fragments`   | `src/analyze_recap_fragments.py`   | `data/recap_fragments.csv`, `data/recap_fragment_associations.csv` | `smiles`, `exclusive_categories` |
| `build`             | `run_analysis_pipeline.py` (inline) | `data/site_ligand_features.csv`         | all of the above (best-effort)          |
| `cluster`           | `src/cluster_protein_similarity.py` | `data/clusters_{all,allosteric}.csv`, `data/largest_clusters.csv` | `build`         |

The `protein_name_pdb` column in the input table comes from `../collection`, not
this pipeline: collection's `protein_names` step writes it and folds it into the
enriched CSV. Analysis inherits the column and never generates it here.

## Scope

The stages above compute and persist per-site features as CSVs under `data/`.
Everything they consume is in `data/site_ligand_features.csv` (the `build` stage output) plus
the stage-specific outputs listed above.

Two of those outputs exist purely as convenient inputs for that step:
`plot_ligand_descriptors.csv` (long — `ccd_code, exclusive_category, metric,
value`) and `plot_ligand_npr.csv` (wide — `NPR1`/`NPR2` columns) are reshaped
views of `ligand_descriptors.csv`, in the form a plotting library consumes
directly. No other stage reads them, and nothing downstream depends on them, so
they can be ignored if you are driving your own analysis from
`site_ligand_features.csv`.

## Internal helpers

The following modules in `src/` are imported by the stage scripts and not
invoked directly:

- `src/cli_utils.py` — common argparse/logging/timing boilerplate
- `src/io_utils.py` — shared CSV loading and DataFrame helpers
- `src/site_pipeline_utils.py` — CIF-driven helpers for structure analysis
