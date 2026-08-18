# Genetic FEEMS — server run

`genetic_feems.py` fits a FEEMS effective-migration surface on the Philippine genetic
dataset (`Phil_1000g_SGDP_1.66M`, PLINK bed/bim/fam). It's the genetic analogue of
`phoneme_feems.ipynb` / `grammar_feems.ipynb`, but as a headless CLI script meant to
run on an external server where the genotype data actually lives — it is never loaded
into this repo.

The PLINK set is assumed already QC'd and pruned. This script does no filtering beyond
selecting samples with known coordinates; it never shells out to a PLINK binary.

## Files to scp to the server

Put these next to each other in one directory on the server (they don't need to be a
git checkout — just the four files below):

| File | From | Size |
|---|---|---|
| `genetic_feems.py` | `python/` | ~15 KB |
| `genetic_feems_env.yml` | `python/` | ~1 KB |
| `philippine_ethnolinguistic_coords.csv` | `data/genetic/initial_datasets/` | ~15 KB |
| `countries.geojson` | `data/` | 14.6 MB |

```bash
scp python/genetic_feems.py python/genetic_feems_env.yml \
    data/genetic/initial_datasets/philippine_ethnolinguistic_coords.csv data/countries.geojson \
    user@server:/path/to/genetic_feems/
```

**Not needed** — `grid_100.shp`/`.shx` (the FEEMS habitat grid) ship inside the `feems`
package itself and are resolved automatically once `genetic_feems_env.yml` is
installed. Only pass `--grid /path/to/grid_100.shp` if the server's `feems` install is
missing its package data.

## Server setup

Two steps, not one — bioconda's `feems` recipe pins `numpy<2` and hasn't been rebuilt
to track `tskit`/`msprime`'s move to NumPy 2.x, so no single conda solve can satisfy
both at once. Let conda solve everything under `feems`'s own `numpy<2` pin first (so
`tskit`/`msprime` install as builds actually compiled against that numpy), then upgrade
`numpy` alone with pip, outside the solver. That direction is safe — NumPy guarantees
code compiled against an older C-API keeps working under a newer NumPy at runtime; it's
only the reverse (newer-compiled code under an older NumPy) that breaks:

```bash
conda env create -f genetic_feems_env.yml   # or mamba, same file
conda activate genetic_feems
pip install --upgrade "numpy>=2,<3"
python -c "import numpy, feems; print(numpy.__version__, 'feems OK')"
```

Don't `pip install --force-reinstall tskit msprime` after the numpy upgrade — that
pulls numpy2-ABI wheels for them and reintroduces
`ImportError: numpy._core.multiarray failed to import`.

## Running

Two stages, run separately so a sample-ID mismatch is caught before spending time on
the fit:

```bash
# 1. Match .fam sample IDs (Phil_<population>) to coordinates — cheap, no .bed read.
python genetic_feems.py --stage match \
    --bfile /path/to/Phil_1000g_SGDP_1.66M \
    --coords philippine_ethnolinguistic_coords.csv \
    --out ./out

# Inspect ./out/GENETIC_feems_sample_match.csv — check matched-population count and
# the "populations with zero matched samples" line before continuing.

# 2. Fit the surface.
python genetic_feems.py --stage fit \
    --bfile /path/to/Phil_1000g_SGDP_1.66M \
    --coords philippine_ethnolinguistic_coords.csv \
    --out ./out --plot
```

Or run both in one call (default `--stage all`):

```bash
python genetic_feems.py --bfile /path/to/Phil_1000g_SGDP_1.66M \
    --coords philippine_ethnolinguistic_coords.csv --out ./out --plot
```

Useful flags:
- `--id-prefix` — defaults to `Phil_`; change if the `.fam` uses a different prefix.
- `--coast-buf` — habitat buffer in degrees, default `0.8` (0.7 strands the Sulu
  Archipelago; the genetic sample set includes Sulu/Tawi-Tawi/Basilan populations).
- `--n-folds` — CV folds, default `0` = leave-one-deme-out. The dominant runtime cost
  is folds × lambda-grid points, not SNP count (FEEMS collapses genotypes to per-deme
  allele frequencies once, up front) — lower this for a quick test run.
- `--max-snps N` — random SNP subsample for a fast smoke test; omit for the real run.
- `--lamb-min` / `--lamb-max` / `--lamb-n` — CV lambda grid; widen if the script warns
  the CV minimum landed on a grid boundary.

## Files to bring back

Copy `./out/*` back into local `data/genetic/feems/`:

```bash
scp 'user@server:/path/to/genetic_feems/out/*' data/genetic/feems/
```

| File | Consumed by |
|---|---|
| `genetic_surface_raster.csv` | the migration-surface raster (`lon,lat,log_w_ratio`) — feeds a future `R/genetic_analysis/[6]_feems_plot_GENETIC_span.R`, same way `R/phoneme_analysis/[6]_feems_plot_PA_span.R` consumes `phoneme_surface_raster.csv` |
| `nodepos_genetic.csv` | diagnostic (grid node positions + sample occupancy) |
| `edgew_genetic.csv` | diagnostic (fitted edge weights) |
| `genetic_feems_meta.json` | run config audit trail (habitat, lambda, sample/SNP counts) |
| `GENETIC_feems_sample_match.csv` | sample-matching audit trail from stage 1 |

Note: there is currently no `R/genetic_analysis/[6]` script that reads
`genetic_surface_raster.csv` — the genetic pipeline's `[7]_GENETIC_geoplots.R` still
plots on a plain grey map. Adding that `[6]` script is a separate follow-up.
