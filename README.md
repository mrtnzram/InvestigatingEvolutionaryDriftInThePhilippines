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
feed distances into **EEMS** (Estimated Effective Migration Surfaces, run in MATLAB) →
overlay a minimum spanning tree (MST) and historical "waypoint" migration routes on a map.

This is a collection of analysis scripts run interactively, not a software package — there
is no build system or test suite.

## Repository layout

```
R/
  phoneme_analysis/   Numbered Ruhlen/PHOIBLE pipeline ([0]→[1]→…→[8]) + EEMS prep + shuffled null
  grammar_analysis/   Numbered GRAMBANK pipeline ([0]→[1]→[2]→[4]) + EEMS prep + shuffled null
  cognate_analysis/   LingPy cognate analysis (gitignored; not used for the paper)
  helpers/            default.eems.plots.R (vendored EEMS plotting helper)
data/                 All R-pipeline data: inputs + intermediates (read/written via here("data", ...))
figures/              Generated plots, grouped by section/type (gitignored)
python/               Waypoint/optimal-path notebooks + LingPy helper (run as-is; not maintained here)
matlab/               EEMS (MATLAB): one folder per analysis; MCMC run outputs are gitignored
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
3. **`*_EEMS.R`** — turn the similarity scores into the EEMS dissimilarity (`datapath.diffs`) and
   coordinate (`datapath.coord`) inputs, plus a shuffled null-model variant.
4. *(MATLAB)* run the EEMS MCMC in `matlab/matlab_eems_*/eems*.m`, consuming the `datapath.*` files.
5. **`eems_plot_*_span.R`** — render the EEMS MCMC output (`reemsplots2`) into a base map RDS.
   For phonemes this is `[7]_eems_plot_PA_span.R`; for grammar `[2]_eems_plot_GA_span.R`.
6. **`*_weight_mst_eems_span.R`** — overlay MST edges / waypoint routes on the base map
   (`[8]_PA_weight_mst_eems_span.R` for phonemes, `[4]_GA_weight_mst_eems_span.R` for grammar).

`*_shuffled` / `eems_plot_*_shuffled.R` are permutation null-model runs used to assess
significance against the real analysis.

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
  `rnaturalearth`(`data`), `maps`, `reemsplots2`.
- **MATLAB/Octave** for the EEMS MCMC step (third-party EEMS code vendored under `matlab/**/mscripts/`).

## Not tracked in git

Figures, the cognate analysis, credentials, archived legacy scripts, large source databases
(`values.csv`, `languages.csv`, `logicalTLI_*`), generated `.rds` objects, and EEMS MCMC run
outputs are gitignored (see `.gitignore`). Regenerate them by running the pipeline.
