# Genetic sPCA — server run

`[4.2]_GENETIC_sPCA_remote.R` is the genetic analogue of `[4.2]_PHONEME_sPCA.R` /
`[4.2]_GRAMMAR_sPCA.R`: residualizes a genotype matrix against the same
population-structure PCs `[4]_GENETIC_PVR.R` uses, then runs real
`adegenet::spca()` on the result. It's a headless R script meant to run on the
server where `phil_only_pruned` (the pruned, Philippines-only PLINK fileset
`pca_results_phil_only.eigenvec` was generated from) actually lives — the same
"scp up, run manually, scp back" pattern as `python/genetic_feems.py`
(`README_genetic_feems.md`), and for the same reason: the genotype data is
never loaded into this repo.

`phil_only_pruned` is assumed already pruned/QC'd; this script never shells out
to a PLINK binary.

## Files to scp to the server

Put these next to each other in one directory (they don't need to be a git
checkout — just the five files below):

| File | From |
|---|---|
| `[4.2]_GENETIC_sPCA_remote.R` | `R/genetic_analysis/` |
| `select_moran_eigenvectors.R` | `R/shared/` |
| `pca_results_phil_only.eigenvec` | `data/pvr/` |
| `GENETIC_final.csv` | `data/network_distance/` |
| `GENETIC_dist_matrix.csv` | `data/network_distance/` |

```bash
scp R/genetic_analysis/"[4.2]_GENETIC_sPCA_remote.R" R/shared/select_moran_eigenvectors.R \
    data/pvr/pca_results_phil_only.eigenvec data/network_distance/GENETIC_final.csv \
    data/network_distance/GENETIC_dist_matrix.csv \
    user@server:/path/to/genetic_spca/
```

`select_moran_eigenvectors.R` must sit in the same directory as the main
script — it's sourced relative to the script's own path, not the working
directory.

## Server setup

```bash
Rscript -e 'install.packages(c("bigsnpr","bigstatsr","spdep","adespatial","adegenet","dplyr","tibble"))'
```

`bigsnpr`/`bigstatsr` compile C++ (RcppArmadillo, OpenMP) — if the install
fails on a missing OpenMP linker symbol, that's a toolchain issue on that
machine, not this script; point `R_MAKEVARS_USER` at a Makevars with a working
`SHLIB_OPENMP_CFLAGS`/OpenMP-linked compiler, same fix this repo's own
`~/.R/Makevars` needed locally.

### If `install.packages()` fails on missing system libraries (sf/igraph/etc.)

`spdep` **Depends** on `sf` (needs GDAL/GEOS/PROJ/libxml2 at the system level);
`adespatial` **Imports** `shiny` (→ `bslib`/`sass`/`fs`) and `adephylo` (→
`phylobase` → `RNeXML` → `XML`/`httr`/`xml2`/`curl`/`openssl`). Both are hard
dependencies of packages this script actually calls (`mem()`, `multispati()`),
so there's no `dependencies = c("Imports")` trim around it — on a machine
without those system libraries (and without root to install them), source
compilation fails, sometimes as a runtime `dyn.load` error on an already-built
`.so` (e.g. `libxml2.so.16: cannot open shared object file`) rather than a
clean compile failure.

Fix: skip source compilation and pull precompiled binaries (system libs
bundled in the env) from conda-forge — no root needed. `bigsnpr`/`bigstatsr`
are both there too, already built with working OpenMP, sidestepping the
Makevars issue above as well.

```bash
# Only if neither is already on PATH:
curl -L -o ~/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash ~/miniforge.sh -b -p ~/miniforge3
source ~/miniforge3/etc/profile.d/conda.sh

mamba create -n genetic-spca -c conda-forge \
  r-base r-bigsnpr r-bigstatsr r-spdep r-adespatial r-adegenet r-dplyr r-tibble \
  -y
conda activate genetic-spca
```

Then run the script from within that environment (below) instead of the
system `Rscript`.

## Running

```bash
Rscript "[4.2]_GENETIC_sPCA_remote.R" \
    --bfile /path/to/phil_only_pruned \
    --eigenvec pca_results_phil_only.eigenvec \
    --genetic-final GENETIC_final.csv \
    --dist-matrix GENETIC_dist_matrix.csv \
    --out ./out
```

`--bfile` takes a PLINK prefix (no `.bed`/`.bim`/`.fam` extension). First run
converts `phil_only_pruned.bed` to a `.rds`/`.bk` file-backed pair next to it
(`bigsnpr::snp_readBed()`) — subsequent runs reuse that `.rds` instead of
re-converting.

Useful flags:
- `--id-prefix` — stripped from `.fam` FID before matching against
  `GENETIC_final.csv`'s `population` column (default: `""`, i.e. FID is
  assumed to already equal the population name, matching the convention
  `pca_results_phil_only.eigenvec` already uses — **not**
  `genetic_feems.py`'s `PHIL_<population>_<code>` convention, which is a
  different, larger fileset). Set this if `phil_only_pruned`'s `.fam` turns
  out to use a prefixed ID instead.
- `--alias` — extra `Old=New` population rename for `.fam` tokens that don't
  match `GENETIC_final.csv` verbatim (e.g. `ManoboRajahKabunsuwan=ManoboRK`,
  already built in, same as `[4]_GENETIC_PVR.R`/`genetic_feems.py`);
  repeatable.
- `--max-snps` — random subsample before residualizing (default: `0` = use
  all SNPs in the pruned panel). SNP count barely affects runtime here (the
  expensive step is an SVD linear in SNP count, not the cubic-in-SNP-count
  step the naive approach would hit — see the script's §6 comment), so this
  is a memory safety valve more than a speed one.

Check the console output after the run: it reports how many individuals had
no population match (should be 0, or a small number you can chase down with
`--id-prefix`/`--alias`) and the chosen spatial-weights threshold.

## Files to bring back

```bash
scp 'user@server:/path/to/genetic_spca/out/*' data/pvr/
```

| File | Consumed by |
|---|---|
| `GENETIC_sPCA_results.csv` | run summary (N, PCs selected, threshold, eigenvalue, variance explained, permutation p) |
| `GENETIC_sPCA_scores.csv` | per-population sPC1 score + coordinates — same shape as `PHONEME_sPCA_scores.csv`/`GRAMMAR_sPCA_scores.csv`, for the same point-symbol plot |
| `GENETIC_sPCA_loadings.csv` | per-SNP loading on sPC1, sorted by \|loading\| |
