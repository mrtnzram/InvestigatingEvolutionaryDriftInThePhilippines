# Investigating evolutionary drift in the Philippines Analysis

Quantifying and modeling Spanish colonial influence on Philippine languages through
computational-linguistic analysis. The project measures the structural similarity of
Philippine languages to a set of baseline languages (Spanish, plus English/Japanese and
unrelated controls) along two independent dimensions, then models and visualizes the
resulting structure geographically.

- **Phonemes** — phoneme inventories from PHOIBLE.
- **Grammar** — typological features from GRAMBANK (Spanish features cross-referenced from WALS).
- **Cognates** — lexical/Swadesh-list similarity via LingPy *(exploratory; excluded from the paper for now)*.

For each dimension the workflow is: build a binary feature dataframe → compute a weighted
cosine-similarity / distance matrix → fit distance-decay regressions and Mantel tests →
estimate an effective-migration surface with **FEEMS** (Fast EEMS, run in Python) →
overlay a minimum spanning tree (MST) and historical "waypoint" migration routes on a map.

> **Note:** the migration surface was migrated from MATLAB **EEMS** to Python **FEEMS**
> (`python/phoneme_feems.ipynb`). The phoneme pipeline is fully on FEEMS; grammar's FEEMS
> port is still pending.

This is a collection of analysis scripts run interactively, not a software package — there
is no build system or test suite.

## Repository layout

```
R/
  phoneme_analysis/   Numbered Ruhlen/PHOIBLE pipeline ([0]→[1]→…→[8]); FEEMS plotting in [6]/[7]
  grammar_analysis/   Numbered GRAMBANK pipeline ([0]→[1]); FEEMS port still pending
  cognate_analysis/   LingPy cognate analysis (gitignored; not used for the paper)
data/                 All R-pipeline data: inputs + intermediates (read/written via here("data", ...))
figures/              Generated plots, grouped by section/type (gitignored)
python/               FEEMS + waypoint/optimal-path notebooks (phoneme_feems.ipynb is the FEEMS engine)
archived_code/        Superseded exploratory scripts, kept for reference (gitignored)
swadeshlist_jsons/    Per-language Swadesh wordlists (cognate input)
```

## Paths convention

All R scripts resolve files with [`here`](https://here.r-lib.org/) anchored to the RStudio
project (`Indp Research Phillipine Languages.Rproj`), e.g. `read_csv(here("data", "PHOIBLEdf_PH.csv"))`.
Open the `.Rproj` (or set the working directory to the project root) before sourcing any script.

## Pipeline order

Scripts pass data through files in `data/`, so order matters. Within each analysis folder the
bracketed prefix is the execution order:

1. **`[0]_*database.R`** — fetch from PHOIBLE/GRAMBANK (via `lingtypology`), filter to the study
   languages, pivot to a binary feature matrix, write the feature matrix + IDF frequency table.
   - Phonemes use the Creanza/Ruhlen source: run **PART A** of `[0]_CREANZA_RUHLENdatabase.R`,
     then `[0]_Phylogenetic_Tree.R` (prunes the MCC tree to the study languages and returns
     `Ph_Languages_pruned`), then **PART B** of the database script. The Grambank unrelated-control
     set is built locally, so no grammar script needs to run first.
2. **Core similarity + downstream analyses.** For phonemes this is split into single-purpose files
   (run in this order): `[1]_PHONEME_cosine_similarity.R` (writes the cosine/cossim matrices) →
   `[2]_PHONEME_cosine_distribution_analysis.R` (ridge/density plots + Friedman/LMM) →
   `[3]_PHONEME_network_distance.R` (waypoint network + per-language land-penalized distance,
   writes `PHONEME_cossim_dist.csv`) → `[4]_PHONEME_regression.R` (distance-decay models) →
   `[5]_PHONEME_MMRR.R` (Dijkstra pairwise distances + multiple matrix regression of
   dissimilarity on geographic and phylogenetic distance) →
   `[6]_PHONEME_PGLS.R` (needs the phylogenetic tree). The grammar side still uses the combined
   `[1]_GRAMMARanalysisweighted_span.R`. (The old monolithic `[1]_PHONEMEanalysisweighted_span.R`
   and `[2]_PHONEMEanalysisweighted_otherlang.R` are retained for reference.)
3. *(Python — FEEMS)* run `python/phoneme_feems.ipynb` to estimate the effective-migration
   surface. It reads the feature matrix + frequency table from `data/` and writes the surface
   raster (`data/phoneme_surface_raster.csv`), node positions (`data/nodepos_phoneme.csv`), and
   run metadata (`data/phoneme_feems_meta.json`).
4. **`[6]_feems_plot_PA_span.R`** — render the FEEMS surface raster + per-language cosine points
   into a ggplot base map, saved as `data/base_plot_phoneme_FEEMS.rds`.
5. **`[7]_PA_weight_mst_feems_span.R`** — overlay MST edges / waypoint routes on the base map
   and mark the colonial capital; saves `figures/phoneme/mst_waypoints/PA_weight_mst_feems.png`.

The grammar analysis has not yet been ported to FEEMS; its migration-surface steps will be
duplicated from the phoneme `[6]`/`[7]` pattern.

## Python ↔ R coupling (important)

The waypoint network is produced in Python and consumed in R. The two notebooks in `python/`
are kept **as-is** (their internal file paths are still Windows paths and must be updated before
re-running):

- `python/waypointsystem.ipynb` writes `nodes.csv`, `edges.csv`.
- `python/optimal_path.ipynb` reads `GRAMFEATURE_match_df.csv`, `PHOIBLE_z_score_df.csv`; writes
  `mst_edges_{GA,PA}.csv`, `smooth_path_{GA,PA}.csv`.

These files live in `data/`. If you re-run the notebooks, point them at `data/` (or move their
outputs there) so the R scripts pick them up.

## Prerequisites

- **R 4.3+** with: `tidyverse`/`dplyr`/`readr`, `here`, `lingtypology`, `geosphere`, `rethinking`,
  `infotheo`, `proxy`, `reshape2`, `ggplot2`, `patchwork`, `igraph`, `sf`/`sfheaders`,
  `rnaturalearth`(`data`), `maps`, `scales`.
- **Python** with `feems` (+ its deps) for the effective-migration-surface step
  (`python/phoneme_feems.ipynb`).

## Not tracked in git

Figures, the cognate analysis, credentials, archived legacy scripts (`archived_code/`), large
source databases (`values.csv`, `languages.csv`, `logicalTLI_*`), and generated `.rds` objects
are gitignored (see `.gitignore`). Regenerate them by running the pipeline.
