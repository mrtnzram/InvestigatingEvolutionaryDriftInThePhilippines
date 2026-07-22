# AGENTS.md

This file provides guidance to Codex (Codex.ai/code) when working with code in this repository.

## What this project is

A computational-linguistics research project on Philippine languages. It quantifies
linguistic similarity between Philippine languages (plus "languages of interest" —
Spanish, English, Japanese — and a set of unrelated control languages) along three
independent dimensions, then visualizes the resulting structure geographically:

- **Phonemes** — phoneme inventories from the PHOIBLE database.
- **Grammar** — typological features from the GRAMBANK database.
- **Cognates** — lexical/Swadesh-list similarity computed with LingPy.

For each dimension the workflow is: build a binary feature dataframe → compute a
(weighted) similarity/distance matrix → estimate an effective-migration surface with
**FEEMS** (Fast EEMS, run in Python) → overlay a minimum spanning tree (MST)
and "waypoint" historical-route reconstruction on a map of the Philippines.

The migration surface was migrated from MATLAB **EEMS** to Python **FEEMS**
(`python/phoneme_feems.ipynb`). The phoneme pipeline is fully on FEEMS; the grammar
side has not been ported yet (its `[6]`/`[7]` FEEMS scripts are still to be built).

This is not a software package — it is a collection of analysis scripts run
interactively. There is no build system, no test suite, and no package manifest.

## Languages / environments

- **R** — primary analysis and plotting, under `R/` (`phoneme_analysis/`,
  `grammar_analysis/`, `cognate_analysis/`). Run interactively in RStudio;
  the project file is `Indp Research Phillipine Languages.Rproj`.
- **Python** — FEEMS migration-surface estimation (`python/phoneme_feems.ipynb`), plus
  LingPy cognate detection and the waypoint/optimal-path notebooks, under `python/`
  (`optimal_path.ipynb`, `waypointsystem.ipynb`, `enhancingdata_lingpy.py`) and
  `R/cognate_analysis/LingPY_analysis.py`.

There is no command to "run the project." Scripts are executed individually
(source an `.R` file in RStudio, or run a notebook cell).

## Directory layout

- `R/` — maintained analysis scripts (see Pipeline order).
- `data/` — all R-pipeline inputs + intermediates; every R script reads/writes here via
  `here("data", ...)`.
- `figures/` — generated plots, grouped by section (`phoneme/`, `grammar/`, `cognate/`,
  `shared/`) and type; gitignored.
- `python/` — FEEMS engine (`phoneme_feems.ipynb`) + waypoint/LingPy notebooks.
- `swadeshlist_jsons/` — Swadesh wordlist inputs for the cognate analysis.
- `archived_code/` — older/exploratory scripts, superseded by `R/`; gitignored.
- `credentials/` — gitignored.

## Pipeline order (important)

The cleaned, numbered pipelines live in `R/phoneme_analysis/` and `R/grammar_analysis/`.
The bracketed prefix is the execution order, and **scripts pass data through CSV/RDS
files in `data/`** (resolved via `here("data", ...)`), so order matters:

1. `[0]_*database.R` — fetch from PHOIBLE/GRAMBANK (via the `lingtypology` package),
   filter to the study languages, pivot to a binary feature matrix, write
   `data/PHOIBLEdf_PH.csv` / `data/GRAMBANKdf_PH.csv` and `data/*_freq.csv`.
2. `[1]_*analysisweighted_span.R` — compute weighted cosine / similarity matrices
   (IDF-style weighting from `*_freq.csv`), build the waypoint network, run Mantel tests
   and distance-decay models; write cosine/dissimilarity matrices and the waypoint plot RDS.
3. `python/phoneme_feems.ipynb` (Python — FEEMS) — estimate the effective-migration
   surface; reads the feature matrix + frequency table from `data/` and writes
   `data/phoneme_surface_raster.csv`, `data/nodepos_phoneme.csv`, and
   `data/phoneme_feems_meta.json`.
4. `[6]_feems_plot_PA_span.R` — render the FEEMS surface raster + per-language cosine
   points into a ggplot base map; save it as `data/base_plot_phoneme_FEEMS.rds`.
5. `[7]_PA_weight_mst_feems_span.R` — load the base map + waypoint plot RDS, overlay MST
   edges / historical routes, mark the capital, and `ggsave` the final figure to
   `figures/phoneme/mst_waypoints/PA_weight_mst_feems.png`.

The grammar side is not yet on FEEMS; its `[6]`/`[7]` FEEMS scripts will be duplicated
from the phoneme pattern (a `phoneme_feems.ipynb`-equivalent for grammar plus the two
R plotting scripts).

Older / exploratory versions of these analyses are kept in `archived_code/` (e.g.
`PHONEMEanalysis.R`, `GRAMMARanalysisweighted.R`); prefer the numbered scripts in `R/`.

Cognate analysis lives in `R/cognate_analysis/` (gitignored, not used for the paper):
`WORDLIST.R` parses the per-language Swadesh JSONs in `swadeshlist_jsons/` into wordlists,
the `LingPY_analysis` scripts convert orthography→IPA and run LexStat cognate detection,
and `Cognate Similarity Analysis.R` / `Heatmap.R` produce the figures.

## Critical gotchas

- **R paths use `here()`, not the working directory.** All R scripts under `R/` resolve
  files via `here("data", ...)` anchored to the `.Rproj`, so they work regardless of which
  folder they live in — but you must open the `.Rproj` (or set the working directory to the
  project root) first. The old Windows absolute paths have been removed from the R scripts.
  **The Python notebooks still contain Windows paths** (left intentionally unedited) — fix
  those before re-running.
- **`data/` is the data bus.** R scripts and the FEEMS notebook read and write intermediate
  CSVs/RDS in `data/`. Don't rename or move these (e.g. `PHOIBLEdf_PH.csv`,
  `GRAMBANKdf_PH.csv`, `*_cossim.csv`, `*_freq.csv`, `phoneme_surface_raster.csv`,
  `nodepos_phoneme.csv`, `base_plot_phoneme_FEEMS.rds`, `*_waypoint_plot.rds`) without
  tracing consumers — several are shared between the Python FEEMS notebook and the R scripts.
- **`GRAMBANKdf_PH.csv` is hand-curated.** `[0]_GRAMBANKdatabase.R` intentionally does **not**
  overwrite it (Spanish is added from WALS by hand); it reads the curated file back instead.
- **`*_cossim_shuffled.csv`** are permutation null-model matrices. They were produced by the
  now-removed EEMS scripts, but `R/cognate_analysis/Heatmap.R` still reads
  `PHONEME_cossim_shuffled.csv`, so the CSVs are kept in `data/` even though they can no
  longer be regenerated.
- **`PA` = Phoneme Analysis, `GA` = Grammar Analysis** throughout filenames.
- `credentials/default_credential_file.json` exists; do not commit secrets or echo
  its contents.

## Key R package dependencies

`lingtypology` (PHOIBLE/GRAMBANK access + `lat.lang`/`long.lang`/`map.feature`),
`tidyverse`/`dplyr`/`readr`, `geosphere`, `rethinking`, `infotheo`, `proxy`,
`reshape2`, `ggplot2`, `patchwork`, `igraph`, `sf`/`sfheaders`, `rnaturalearth`,
`scales`, and `maps` (for the FEEMS base map). Python side uses `feems` (migration
surface), plus `lingpy`, `segments`, `pandas`.
