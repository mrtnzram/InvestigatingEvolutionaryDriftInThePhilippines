#!/usr/bin/env python3
"""
genetic_feems.py — FEEMS effective-migration surface estimation on Philippine
genetic PLINK data, run headless on a remote server over SSH.

Inputs:  --bfile PLINK prefix (bed/bim/fam, already QC'd/pruned — this script never
         shells out to a PLINK binary), --coords population coordinates CSV.
Outputs: genetic_surface_raster.csv, nodepos_genetic.csv, edgew_genetic.csv,
         genetic_feems_meta.json, GENETIC_feems_sample_match.csv. With --permute,
         also: genetic_surface_raster_null_mean.csv, genetic_null_density.png,
         genetic_feems_null_meta.json (and genetic_null_surface.png if --plot too).
         With --plot: genetic_feems_cv_curve.png, genetic_feems_surface.png.

Usage:
    python genetic_feems.py --stage match --bfile /path/to/Phil_1000g_SGDP_1.66M \
        --coords philippine_ethnolinguistic_coords.csv --out ./out

    python genetic_feems.py --stage fit --bfile /path/to/Phil_1000g_SGDP_1.66M \
        --coords philippine_ethnolinguistic_coords.csv --out ./out

    python genetic_feems.py --stage all ...   # match, then fit, in one call (default)

    # Geographic-shuffle null model (100 permutations, fixed lambda, no CV per
    # permutation) run immediately after a real fit:
    python genetic_feems.py --stage all --permute --bfile ... --coords ... --out ./out

    # Re-run just the null model later against an already-fit --out directory,
    # skipping CV by passing the lamb_cv it already found:
    python genetic_feems.py --stage fit --permute --lamb 0.0123 \\
        --bfile ... --coords ... --out ./out

"""

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# feems (bioconda, still pinned to numpy<2) references numpy aliases removed in
# NumPy 2.0 (np.Inf, np.NaN, np.float, ...); restore them here, before feems
# imports below, rather than patching the installed package.
_NP2_REMOVED_ALIASES = {
    "Inf": np.inf,
    "Infinity": np.inf,
    "NINF": -np.inf,
    "PINF": np.inf,
    "NAN": np.nan,
    "NaN": np.nan,
    "float": float,
    "int": int,
    "bool": bool,
    "complex": complex,
    "long": int,
    "unicode": str,
}
# np.object / np.str: still present in numpy>=2 (deprecated, warns on use), not shimmed.
for _name, _val in _NP2_REMOVED_ALIASES.items():
    if not hasattr(np, _name):
        setattr(np, _name, _val)

FAM_COLS = ["fid", "iid", "father", "mother", "sex", "phenotype"]


# --------------------------------------------------------------------------- #
# Stage: match — resolve .fam sample IDs to population coordinates; text-only,
# never touches .bed.
# --------------------------------------------------------------------------- #


def read_fam(bfile: str) -> pd.DataFrame:
    fam_path = Path(f"{bfile}.fam")
    if not fam_path.exists():
        sys.exit(f"error: .fam file not found at {fam_path}")
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None, names=FAM_COLS, dtype=str)
    fam.insert(0, "row", np.arange(len(fam)))  # original .fam order == G's sample axis
    return fam


def strip_prefix(s: str, prefix: str) -> str:
    if prefix and s.lower().startswith(prefix.lower()):
        return s[len(prefix) :]
    return s


# .fam population tokens don't always match the coords table's `population` verbatim
# — same renames as downstream R (e.g. R/genetic_analysis/[0]_GENETIC_ADMX_MALDER.R);
# keyed lowercase, extend via --alias rather than editing this table for one-offs.
DEFAULT_POPULATION_ALIASES = {
    "manoborajahkabunsuwan": "ManoboRK",
}


def population_candidates(sample_id: str, prefix: str) -> list:
    """.fam IDs are PHIL_<Population>_<SampleCode>; population names have no
    underscores, so try the first underscore token before the full remainder."""
    stripped = strip_prefix(sample_id, prefix)
    candidates = [stripped]
    if "_" in stripped:
        candidates.insert(0, stripped.split("_", 1)[0])
    return candidates


def match_samples(
    fam: pd.DataFrame, coords: pd.DataFrame, id_prefix: str, aliases: dict
) -> pd.DataFrame:
    pop_lookup = {str(p).lower(): p for p in coords["population"]}
    coords_by_pop = coords.set_index("population")

    matched_pop, matched_lat, matched_lon, matched = [], [], [], []
    for fid, iid in zip(fam["fid"], fam["iid"]):
        candidates = population_candidates(fid, id_prefix) + population_candidates(
            iid, id_prefix
        )
        hit = None
        for cand in candidates:
            cand_lower = cand.lower()
            key = pop_lookup.get(cand_lower)
            if key is None and cand_lower in aliases:
                key = pop_lookup.get(aliases[cand_lower].lower())
            if key is not None:
                hit = key
                break
        if hit is None:
            matched_pop.append(None)
            matched_lat.append(np.nan)
            matched_lon.append(np.nan)
            matched.append(False)
        else:
            row = coords_by_pop.loc[hit]
            matched_pop.append(hit)
            matched_lat.append(float(row["latitude"]))
            matched_lon.append(float(row["longitude"]))
            matched.append(True)

    out = fam.copy()
    out["population"] = matched_pop
    out["matched"] = matched
    out["latitude"] = matched_lat
    out["longitude"] = matched_lon
    return out


def parse_aliases(alias_args) -> dict:
    aliases = dict(DEFAULT_POPULATION_ALIASES)
    for spec in alias_args or []:
        if "=" not in spec:
            sys.exit(f"error: --alias '{spec}' must be in the form Old=New")
        old, new = spec.split("=", 1)
        aliases[old.strip().lower()] = new.strip()
    return aliases


def run_match(args) -> pd.DataFrame:
    fam = read_fam(args.bfile)
    coords = pd.read_csv(args.coords)
    if "population" not in coords.columns:
        sys.exit(f"error: {args.coords} has no 'population' column")

    aliases = parse_aliases(args.alias)
    match_df = match_samples(fam, coords, args.id_prefix, aliases)

    n_total = len(match_df)
    n_matched = int(match_df["matched"].sum())
    matched_pops = set(match_df.loc[match_df["matched"], "population"])
    zero_hit_pops = sorted(set(coords["population"]) - matched_pops)

    print(
        f"[match] .fam rows: {n_total} | matched: {n_matched} | dropped: {n_total - n_matched}"
    )
    print(
        f"[match] distinct populations matched: {len(matched_pops)} of {coords['population'].nunique()} in coords table"
    )
    if zero_hit_pops:
        print(
            f"[match] populations with zero matched samples ({len(zero_hit_pops)}): {zero_hit_pops}"
        )

    if n_matched < 10:
        sys.exit(
            f"error: only {n_matched} samples matched — this almost certainly means "
            f"the ID convention doesn't match --id-prefix '{args.id_prefix}' or the "
            f"coords table's 'population' column. Inspect the .fam and coords file "
            f"before proceeding; refusing to fit a surface on this few demes."
        )

    args.out.mkdir(parents=True, exist_ok=True)
    match_path = args.out / "GENETIC_feems_sample_match.csv"
    match_df.drop(columns=["row"]).to_csv(match_path, index=False)
    print(f"[match] wrote {match_path}")
    return match_df


# --------------------------------------------------------------------------- #
# Stage: fit — genotypes, habitat, graph, CV, fit, rasterize, export.
# --------------------------------------------------------------------------- #


def build_habitat(geojson_path: Path, country_name: str, coast_buf: float):
    from shapely.geometry import shape

    with open(geojson_path) as f:
        countries = json.load(f)
    feat = next(
        (
            f
            for f in countries["features"]
            if f["properties"].get("name") == country_name
        ),
        None,
    )
    if feat is None:
        sys.exit(f"error: no feature named '{country_name}' in {geojson_path}")
    geom = shape(feat["geometry"])
    habitat = geom.buffer(coast_buf).simplify(0.1)
    if habitat.geom_type != "Polygon":
        sys.exit(
            f"error: habitat is {habitat.geom_type} with {len(habitat.geoms)} pieces at "
            f"coast_buf={coast_buf} — raise --coast-buf until it merges into one polygon"
        )
    return np.array(habitat.exterior.coords)


def resolve_grid_path(grid_arg, grid_res: str) -> str:
    if grid_arg:
        if not os.path.exists(grid_arg):
            sys.exit(f"error: --grid path not found: {grid_arg}")
        return grid_arg
    import feems

    grid_path = os.path.join(os.path.dirname(feems.__file__), "data", grid_res)
    if not os.path.exists(grid_path):
        sys.exit(
            f"error: grid shapefile not found at {grid_path} (installed feems package) "
            f"and no --grid override given"
        )
    return grid_path


def resolve_geojson_path(geojson_arg) -> Path:
    if geojson_arg:
        p = Path(geojson_arg)
        if not p.exists():
            sys.exit(f"error: --geojson path not found: {p}")
        return p
    here = Path(__file__).resolve().parent
    candidates = [here / "countries.geojson"]
    try:
        import feems

        candidates.append(
            Path(os.path.dirname(feems.__file__)) / "data" / "countries.geojson"
        )
    except ImportError:
        pass
    for c in candidates:
        if c.exists():
            return c
    sys.exit(
        "error: countries.geojson not found next to the script or in the feems package; pass --geojson"
    )


def load_genotypes(bfile: str, match_df: pd.DataFrame, max_snps, seed: int):
    from pandas_plink import read_plink
    from sklearn.impute import SimpleImputer

    print(f"[fit] reading PLINK set: {bfile}")
    bim, fam, G = read_plink(bfile, verbose=False)
    n_variants, n_samples_total = G.shape

    if len(fam) != len(match_df):
        sys.exit(
            f"error: .fam row count from pandas_plink ({len(fam)}) != match table row "
            f"count ({len(match_df)}) — was the match table built from this same bfile?"
        )

    matched_rows = match_df.loc[match_df["matched"], "row"].to_numpy()
    n_matched = len(matched_rows)
    n_snps_total = n_variants

    if max_snps and max_snps < n_variants:
        rng = np.random.default_rng(seed)
        snp_idx = np.sort(rng.choice(n_variants, size=max_snps, replace=False))
        print(
            f"[fit] --max-snps subsample: {max_snps} of {n_variants} SNPs (seed={seed})"
        )
    else:
        snp_idx = slice(None)

    est_bytes = n_matched * (max_snps or n_variants) * 8
    print(
        f"[fit] genotype matrix: {n_matched} samples x {max_snps or n_variants} SNPs "
        f"(~{est_bytes / 1e9:.2f} GB as float64)"
    )

    # Subset both axes before compute() so only the needed genotypes are materialised.
    G_sub = G[snp_idx, :][:, matched_rows]
    genotypes = np.array(G_sub).T  # (n_matched, n_snps)

    n_missing = np.isnan(genotypes).sum()
    if n_missing:
        print(
            f"[fit] mean-imputing {n_missing} missing calls "
            f"({n_missing / genotypes.size * 100:.3f}% of matrix)"
        )
        n_before_impute = genotypes.shape[1]
        imp = SimpleImputer(missing_values=np.nan, strategy="mean")
        genotypes = imp.fit_transform(genotypes)  # drops SNPs with zero observed calls
        n_all_missing = n_before_impute - genotypes.shape[1]
        if n_all_missing:
            print(
                f"[fit] dropped {n_all_missing} SNPs with zero observed calls across matched samples"
            )

    # feems requires polymorphic SNPs and crashes on monomorphic ones (a broken debug
    # print, not a graceful error), so filter using feems' own invariance definition.
    variant_sum = genotypes.sum(axis=0)
    invariant = (variant_sum == 0) | (variant_sum == 2 * genotypes.shape[0])
    if invariant.any():
        print(
            f"[fit] dropping {int(invariant.sum())} of {genotypes.shape[1]} SNPs as invariant "
            f"(monomorphic across matched samples)"
        )
        genotypes = genotypes[:, ~invariant]

    return genotypes, n_snps_total, matched_rows


def rasterize(sp_graph, outer, raster_res: int):
    import matplotlib.tri as mtri
    from shapely.geometry import Point, Polygon

    n_nodes = sp_graph.node_pos.shape[0]
    edge_arr = np.array(sp_graph.edges)  # 0-indexed, (n_edges, 2)
    log_w = np.log10(sp_graph.w / sp_graph.w.mean())

    node_val = np.zeros(n_nodes)
    node_cnt = np.zeros(n_nodes)
    for k in range(len(edge_arr)):
        i, j = edge_arr[k]
        node_val[i] += log_w[k]
        node_cnt[i] += 1
        node_val[j] += log_w[k]
        node_cnt[j] += 1
    node_val = np.divide(
        node_val, node_cnt, out=np.full(n_nodes, np.nan), where=node_cnt > 0
    )

    triang = mtri.Triangulation(sp_graph.node_pos[:, 0], sp_graph.node_pos[:, 1])
    interp = mtri.LinearTriInterpolator(triang, node_val)

    lon_g = np.linspace(outer[:, 0].min(), outer[:, 0].max(), raster_res)
    lat_g = np.linspace(outer[:, 1].min(), outer[:, 1].max(), int(raster_res * 4 / 3))
    LON, LAT = np.meshgrid(lon_g, lat_g)
    Z = interp(LON, LAT)

    poly = Polygon(outer)
    mask = np.array(
        [
            [poly.contains(Point(LON[i, j], LAT[i, j])) for j in range(LON.shape[1])]
            for i in range(LON.shape[0])
        ]
    )
    Z = np.ma.array(Z, mask=~mask | np.ma.getmaskarray(Z))
    print(f"[fit] valid raster fraction: {(~np.ma.getmaskarray(Z)).mean():.3f}")
    return LON, LAT, Z


def _rasterize_fixed_grid(sp_graph, LON, LAT, valid_mask):
    """Same node-averaging + triangulation as rasterize(), but reuses a precomputed
    LON/LAT grid and point-in-polygon mask instead of rebuilding them. Valid because
    node_pos (and hence the mask) is a property of grid/edges/outer alone — it never
    changes across null permutations, only the fitted edge weights do. Used by
    run_permutations() to avoid repeating rasterize()'s point-in-polygon loop
    (the expensive part, at raster_res=300 a few hundred thousand Point.contains
    calls) n_perm times."""
    import matplotlib.tri as mtri

    n_nodes = sp_graph.node_pos.shape[0]
    edge_arr = np.array(sp_graph.edges)
    log_w = np.log10(sp_graph.w / sp_graph.w.mean())

    node_val = np.zeros(n_nodes)
    node_cnt = np.zeros(n_nodes)
    for k in range(len(edge_arr)):
        i, j = edge_arr[k]
        node_val[i] += log_w[k]
        node_cnt[i] += 1
        node_val[j] += log_w[k]
        node_cnt[j] += 1
    node_val = np.divide(
        node_val, node_cnt, out=np.full(n_nodes, np.nan), where=node_cnt > 0
    )

    triang = mtri.Triangulation(sp_graph.node_pos[:, 0], sp_graph.node_pos[:, 1])
    interp = mtri.LinearTriInterpolator(triang, node_val)
    Z = interp(LON, LAT)
    return np.ma.array(Z, mask=~valid_mask | np.ma.getmaskarray(Z))


def is_degenerate_fit(sp_graph, tol=1e-9):
    """True if every edge got (numerically) the same weight — the signature of a
    FEEMS fit that never actually fit anything.

    This happens when fit_null_model's constant-w MLE is unbounded (the likelihood
    keeps improving as w -> inf), so Nelder-Mead walks out to an arbitrary huge w0,
    L-BFGS starts there with non-finite gradients and returns its starting point
    unchanged. Because fit_null_model does not depend on lamb, EVERY lambda then
    returns the same uniform surface and the CV curve is perfectly flat. Observed
    on the 26-language GRAMBANK subset (16 demes, 10 singletons); guarded here
    because the failure is otherwise silent — it yields a plausible-looking flat
    map rather than an error.
    """
    w = np.asarray(sp_graph.w)
    return bool(np.ptp(w) <= tol * max(abs(float(np.mean(w))), 1.0))


def run_permutations(args, sp_graph, genotypes, coord, grid, edges, outer, LON, LAT, Z, lamb_cv):
    """Geographic-shuffle null model: reassign coord rows onto genotypes in random
    order (genotypes never reordered — only which sample owns which lon/lat is
    shuffled), refit, and summarise across args.n_perm permutations.

    Test statistic. Two are recorded:
      * T_het = sd of log10(w/wbar) over valid raster cells. Simple, but it is a
        function of the fitted w(lamb) and therefore inherits lamb's shrinkage:
        as lamb -> inf, T_het -> 0 regardless of signal. Holding lamb fixed across
        permutations embeds that nuisance parameter in the statistic, so a
        T_het-based p-value is only interpretable at a fixed, shared lamb.
      * T_cv  = min over lamb of the LODO CV error (--perm-cv). Profiling over lamb
        removes the nuisance parameter, and it tests the question the null model is
        actually asking: does the sample->location correspondence improve
        out-of-sample prediction? Costs one CV sweep per permutation.

    Note on the averaged null raster: independent permutations produce independently
    oriented surfaces, so their mean cancels toward 0 at roughly
    (per-permutation sd)/sqrt(n_perm) whether or not geography alone can generate
    structure. The averaged map is therefore NOT evidence of a flat null — the
    per-cell null sd (also exported) and the per-permutation statistic are the
    informative quantities.

    Resumable: progress is checkpointed to {args.out}/genetic_null_checkpoint.npz
    every 10 iterations, since this is the heaviest of the three FEEMS pipelines
    (1000+ samples) and may run unattended on a server. Pass --perm-restart to
    ignore an existing checkpoint and start over.
    """
    import time

    from feems import SpatialGraph

    n = coord.shape[0]
    valid_mask = ~np.ma.getmaskarray(Z)  # grid/outer/node_pos are permutation-invariant
    vmax_real = float(np.nanmax(np.abs(Z)))
    # Fixed histogram bins (not a growing pooled array) so density accumulates in
    # O(1) memory per iteration and survives a checkpoint/resume cleanly.
    hist_edges = np.linspace(-3 * vmax_real, 3 * vmax_real, 161)

    ckpt_path = args.out / "genetic_null_checkpoint.npz"
    start_i = 0
    raster_sum = np.zeros_like(LON, dtype=float)
    raster_sq = np.zeros_like(LON, dtype=float)  # for the per-cell null sd
    raster_count = np.zeros_like(LON, dtype=int)
    hist_counts = np.zeros(len(hist_edges) - 1, dtype=np.int64)
    perm_het = []  # per-permutation T_het
    perm_cv = []  # per-permutation T_cv (only when --perm-cv)
    n_completed = 0
    n_skipped = 0
    skipped_log = []

    if ckpt_path.exists() and not args.perm_restart:
        ckpt = np.load(ckpt_path)
        if int(ckpt["n_perm"]) == args.n_perm and int(ckpt["perm_seed"]) == args.perm_seed:
            start_i = int(ckpt["next_i"])
            raster_sum = ckpt["raster_sum"]
            raster_sq = ckpt["raster_sq"]
            raster_count = ckpt["raster_count"]
            hist_counts = ckpt["hist_counts"]
            perm_het = list(ckpt["perm_het"])
            perm_cv = list(ckpt["perm_cv"])
            n_completed = int(ckpt["n_completed"])
            n_skipped = int(ckpt["n_skipped"])
            print(f"[permute] resuming from checkpoint: {start_i}/{args.n_perm} already attempted")
        else:
            print(
                "[permute] checkpoint doesn't match --n-perm/--perm-seed — ignoring it "
                "(pass --perm-restart to silence this and always start fresh)"
            )

    if args.perm_cv:
        from feems.cross_validation import run_cv

        perm_lamb_grid = np.geomspace(args.lamb_min, args.lamb_max, args.lamb_n)[::-1]
        n_folds_perm = args.n_folds if args.n_folds > 0 else sp_graph.n_observed_nodes
        print(
            f"[permute] {args.n_perm} permutations, CV re-run per permutation "
            f"({args.lamb_n} lambdas x {n_folds_perm} folds each) — statistic: T_cv (lambda-free)"
        )
    else:
        print(
            f"[permute] {args.n_perm} permutations at FIXED lamb={lamb_cv:.6g} — statistic: T_het.\n"
            f"[permute] NOTE: T_het inherits lambda's shrinkage; a fixed-lambda p-value is only "
            f"comparable within this lambda. Use --perm-cv for the lambda-free T_cv statistic."
        )

    t0 = time.time()
    for i in range(start_i, args.n_perm):
        rng = np.random.default_rng(args.perm_seed + i)
        perm = rng.permutation(n)
        coord_null = coord[perm]
        try:
            sp_graph_null = SpatialGraph(genotypes, coord_null, grid, edges, scale_snps=True)
            # occupancy is a property of the fixed point set, not the label assignment —
            # this should never fire; it's a cheap sanity check against an implementation bug.
            assert sp_graph_null.n_observed_nodes == sp_graph.n_observed_nodes
            if args.perm_cv:
                cv_err_p = run_cv(
                    sp_graph_null, perm_lamb_grid, n_folds=n_folds_perm, factr=1e10
                )
                mean_err_p = np.mean(cv_err_p, axis=0).ravel()
                lamb_p = float(perm_lamb_grid[np.argmin(mean_err_p)])
                perm_cv.append(float(mean_err_p.min()))
            else:
                lamb_p = lamb_cv
            sp_graph_null.fit(lamb=lamb_p, optimize_q=None)
            if is_degenerate_fit(sp_graph_null):
                raise RuntimeError(
                    "degenerate fit (all edge weights equal) — the model did not "
                    "converge for this permutation"
                )
            Z_null = _rasterize_fixed_grid(sp_graph_null, LON, LAT, valid_mask)
        except Exception as e:
            n_skipped += 1
            skipped_log.append((i, str(e)))
            continue

        ok = ~np.ma.getmaskarray(Z_null)
        vals = Z_null[ok].data
        raster_sum[ok] += vals
        raster_sq[ok] += vals ** 2
        raster_count[ok] += 1
        counts, _ = np.histogram(vals, bins=hist_edges)
        hist_counts += counts
        perm_het.append(float(vals.std()))
        n_completed += 1

        if i == start_i:
            per_iter = time.time() - t0
            remaining = args.n_perm - start_i
            print(
                f"[permute] first iteration took {per_iter:.1f}s — "
                f"est. total ~{per_iter * remaining:.0f}s for the remaining {remaining}"
            )
        if (i + 1) % 10 == 0 or (i + 1) == args.n_perm:
            print(f"[permute] {i + 1}/{args.n_perm} done ({time.time() - t0:.0f}s elapsed)")
            np.savez(
                ckpt_path,
                raster_sum=raster_sum,
                raster_sq=raster_sq,
                raster_count=raster_count,
                hist_counts=hist_counts,
                hist_edges=hist_edges,
                perm_het=np.asarray(perm_het, dtype=float),
                perm_cv=np.asarray(perm_cv, dtype=float),
                next_i=i + 1,
                n_completed=n_completed,
                n_skipped=n_skipped,
                n_perm=args.n_perm,
                perm_seed=args.perm_seed,
            )

    print(f"[permute] completed {n_completed}/{args.n_perm} ({n_skipped} skipped): {skipped_log}")
    if n_completed == 0:
        sys.exit("error: no permutation completed — nothing to summarise")
    skip_frac = n_skipped / max(args.n_perm, 1)
    if skip_frac > 0.05:
        print(
            f"[permute] WARNING: {skip_frac:.0%} of permutations failed to converge and were "
            f"dropped. The p-values below are therefore CONDITIONAL ON CONVERGENCE, not a "
            f"clean permutation test: non-converged nulls are degenerate (uniform w, T_het=0), "
            f"so dropping them removes the low end of the null distribution and inflates the "
            f"T_het p-value. Report the skip rate alongside any p-value, and treat a skip rate "
            f"this high as evidence the model is poorly identified on this data."
        )

    null_mean = np.divide(
        raster_sum, raster_count, out=np.full_like(raster_sum, np.nan), where=raster_count > 0
    )
    # Per-cell null sd. This, not the mean, is what says how much structure geography
    # alone can generate at each location: independent permutations cancel in the mean
    # (~sd/sqrt(n_perm)) but not in the spread.
    null_var = np.divide(
        raster_sq, raster_count, out=np.full_like(raster_sq, np.nan), where=raster_count > 0
    ) - null_mean ** 2
    null_sd = np.sqrt(np.clip(null_var, 0.0, None))
    valid_mean = raster_count > 0

    # Per-cell z of the observed surface against the null, so the map can be read
    # cell-by-cell instead of collapsing to one scalar.
    with np.errstate(divide="ignore", invalid="ignore"):
        z_obs = np.where(null_sd > 0, (np.asarray(Z) - null_mean) / null_sd, np.nan)

    null_raster_df = pd.DataFrame(
        {
            "lon": LON[valid_mean],
            "lat": LAT[valid_mean],
            "log_w_ratio_null_mean": null_mean[valid_mean],
            "log_w_ratio_null_sd": null_sd[valid_mean],
            "log_w_ratio_obs": np.asarray(Z)[valid_mean],
            "z_obs_vs_null": z_obs[valid_mean],
        }
    )
    null_raster_path = args.out / "genetic_surface_raster_null_mean.csv"
    null_raster_df.to_csv(null_raster_path, index=False)
    print(f"[permute] wrote {null_raster_path} ({null_raster_df.shape})")

    # --- permutation test -----------------------------------------------------
    perm_het_arr = np.asarray(perm_het, dtype=float)
    obs_het = float(np.asarray(Z)[valid_mask].std())
    p_het = (1 + int(np.sum(perm_het_arr >= obs_het))) / (len(perm_het_arr) + 1)
    print(
        f"[permute] T_het observed {obs_het:.5f} | null mean {perm_het_arr.mean():.5f} "
        f"max {perm_het_arr.max():.5f} | p = {p_het:.4f}"
    )
    p_cv = None
    if args.perm_cv and perm_cv:
        perm_cv_arr = np.asarray(perm_cv, dtype=float)
        # obs_cv is the real run's profiled CV minimum, passed through args by run_fit.
        obs_cv = getattr(args, "_obs_cv_min", None)
        if obs_cv is not None:
            p_cv = (1 + int(np.sum(perm_cv_arr <= obs_cv))) / (len(perm_cv_arr) + 1)
            print(
                f"[permute] T_cv  observed {obs_cv:.6f} | null mean {perm_cv_arr.mean():.6f} "
                f"min {perm_cv_arr.min():.6f} | p = {p_cv:.4f}"
            )
        else:
            print("[permute] T_cv null collected but no observed CV minimum available "
                  "(run with --stage fit/all and without --lamb to enable the T_cv test)")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bin_centers = (hist_edges[:-1] + hist_edges[1:]) / 2
    bin_width = hist_edges[1] - hist_edges[0]
    total = hist_counts.sum()
    density = hist_counts / (total * bin_width) if total else hist_counts.astype(float)

    fig, ax = plt.subplots(figsize=(5, 3), dpi=120)
    ax.bar(bin_centers, density, width=bin_width, alpha=0.7)
    ax.set_xlabel(r"null $\log_{10}(w/\bar{w})$ (pooled raster cells, all permutations)")
    ax.set_ylabel("density")
    density_path = args.out / "genetic_null_density.png"
    fig.savefig(density_path, bbox_inches="tight")
    plt.close(fig)
    print(f"[permute] wrote {density_path}")

    if args.plot:
        from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

        feems_cmap = LinearSegmentedColormap.from_list(
            "feems_oc", ["#E8890C", "#FFFFFF", "#00C6D6"]
        )
        norm = TwoSlopeNorm(vmin=-vmax_real, vcenter=0.0, vmax=vmax_real)
        fig, ax = plt.subplots(figsize=(7, 9), dpi=150)
        ax.pcolormesh(
            LON, LAT, np.ma.array(null_mean, mask=~valid_mean),
            cmap=feems_cmap, norm=norm, shading="auto", alpha=0.85,
        )
        ax.scatter(coord[:, 0], coord[:, 1], s=12, c="black")
        ax.set_aspect("equal")
        surface_path = args.out / "genetic_null_surface.png"
        fig.savefig(surface_path, bbox_inches="tight")
        plt.close(fig)
        print(f"[permute] wrote {surface_path}")

    null_meta = {
        "n_perm": args.n_perm,
        "n_completed": n_completed,
        "n_skipped": n_skipped,
        "perm_seed": args.perm_seed,
        "lamb_handling": "re-selected per permutation (--perm-cv)" if args.perm_cv
                        else f"fixed at {lamb_cv}",
        "lamb_cv_reused": lamb_cv,
        "T_het_observed": obs_het,
        "T_het_null_mean": float(perm_het_arr.mean()),
        "T_het_null_max": float(perm_het_arr.max()),
        "T_het_p": p_het,
        "T_cv_p": p_cv,
    }
    null_meta_path = args.out / "genetic_feems_null_meta.json"
    with open(null_meta_path, "w") as f:
        json.dump(null_meta, f, indent=2)
    print(f"[permute] wrote {null_meta_path}")


def run_fit(args, match_df: pd.DataFrame):
    import networkx as nx
    from feems import SpatialGraph
    from feems.cross_validation import run_cv
    from feems.utils import prepare_graph_inputs
    from shapely.geometry import MultiPoint, Point

    genotypes, n_snps_total, matched_rows = load_genotypes(
        args.bfile, match_df, args.max_snps, args.seed
    )

    matched = match_df.loc[match_df["matched"]].set_index("row").loc[matched_rows]
    coord = matched[["longitude", "latitude"]].to_numpy()  # FEEMS wants (lon, lat)
    n_populations = matched["population"].nunique()
    print(f"[fit] {len(coord)} samples across {n_populations} populations")

    geojson_path = resolve_geojson_path(args.geojson)
    outer = build_habitat(geojson_path, args.country, args.coast_buf)
    print(
        f"[fit] habitat from {geojson_path.name} ({args.country}), coast_buf={args.coast_buf}, "
        f"{len(outer)} outer vertices"
    )

    grid_path = resolve_grid_path(args.grid, args.grid_res)
    outer, edges, grid, ipmap = prepare_graph_inputs(
        coord=coord,
        ggrid=grid_path,
        translated=False,
        buffer=0,
        outer=outer,
    )

    G_nx = nx.Graph()
    G_nx.add_nodes_from(range(grid.shape[0]))
    G_nx.add_edges_from([(a - 1, b - 1) for a, b in edges])
    n_components = nx.number_connected_components(G_nx)
    print(
        f"[fit] grid nodes: {grid.shape[0]} | edges: {edges.shape[0]} | components: {n_components}"
    )
    if n_components != 1:
        sys.exit(
            f"error: habitat fragmented into {n_components} disconnected pieces — raise --coast-buf"
        )

    node_hull = MultiPoint([tuple(g) for g in grid]).convex_hull
    stranded = [i for i, c in enumerate(coord) if not node_hull.contains(Point(c))]
    print(f"[fit] stranded samples: {len(stranded)}")
    if stranded:
        sys.exit(
            f"error: {len(stranded)} sample(s) fall outside the node mesh: {stranded}"
        )

    occ_idx = np.unique(ipmap)
    print(f"[fit] occupied demes: {len(occ_idx)} of {grid.shape[0]}")

    # scale_snps=True: FEEMS default for real genotypes (Patterson-scaled allele frequencies).
    sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)

    if args.lamb is not None:
        # Explicit lambda override — skip the expensive leave-one-deme-out CV sweep
        # entirely. This is what lets --permute be re-run later against an already-
        # fit --out directory without repeating CV (pass the lamb_cv a prior run
        # already found, e.g. read from that run's genetic_feems_meta.json).
        lamb_cv = float(args.lamb)
        print(f"[fit] --lamb given: skipping CV, lamb_cv = {lamb_cv:.6g}")
    else:
        n_folds = args.n_folds if args.n_folds > 0 else sp_graph.n_observed_nodes
        lamb_grid = np.geomspace(args.lamb_min, args.lamb_max, args.lamb_n)[::-1]
        print(f"[fit] running CV: {args.lamb_n} lambdas x {n_folds} folds")
        cv_err = run_cv(sp_graph, lamb_grid, n_folds=n_folds, factr=1e10)
        mean_cv_err = np.mean(cv_err, axis=0)
        argmin = int(np.argmin(mean_cv_err))
        lamb_cv = float(lamb_grid[argmin])
        # observed profiled CV minimum — the T_cv test statistic for --perm-cv
        args._obs_cv_min = float(mean_cv_err.min())
        print(f"[fit] lamb_cv = {lamb_cv:.6g}")
        if argmin in (0, len(lamb_grid) - 1):
            print(
                "[fit] WARNING: CV minimum at lambda-grid boundary — widen --lamb-min/--lamb-max"
            )

        if args.plot:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(figsize=(5, 3), dpi=120)
            ax.plot(np.log10(lamb_grid), mean_cv_err, "o-")
            ax.axvline(np.log10(lamb_cv), ls="--", c="r")
            ax.set_xlabel(r"$\log_{10}\lambda$")
            ax.set_ylabel("LOO CV error")
            fig.savefig(args.out / "genetic_feems_cv_curve.png", bbox_inches="tight")
            plt.close(fig)

    sp_graph.fit(lamb=lamb_cv, optimize_q=None)
    if is_degenerate_fit(sp_graph):
        sys.exit(
            "error: degenerate fit — every edge received the same weight, so no migration "
            "surface was actually estimated. This means the constant-w null model had no "
            "finite MLE (w -> inf) and every lambda returns the same uniform surface; check "
            "whether the CV curve is flat and how many demes are occupied. Refusing to "
            "export a surface that carries no information."
        )
    log_ratio = np.log10(sp_graph.w / sp_graph.w.mean())
    print(
        f"[fit] log10(w/wbar): min {log_ratio.min():.3f} max {log_ratio.max():.3f} sd {log_ratio.std():.3f}"
    )

    LON, LAT, Z = rasterize(sp_graph, outer, args.raster_res)

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm

        feems_cmap = LinearSegmentedColormap.from_list(
            "feems_oc", ["#E8890C", "#FFFFFF", "#00C6D6"]
        )
        vmax = float(np.nanmax(np.abs(Z)))
        norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)
        fig, ax = plt.subplots(figsize=(7, 9), dpi=150)
        ax.pcolormesh(
            LON, LAT, Z, cmap=feems_cmap, norm=norm, shading="auto", alpha=0.85
        )
        ax.scatter(coord[:, 0], coord[:, 1], s=12, c="black")
        ax.set_aspect("equal")
        fig.savefig(args.out / "genetic_feems_surface.png", bbox_inches="tight")
        plt.close(fig)

    # --- exports: schema matches the other FEEMS pipelines' outputs in data/ ---
    args.out.mkdir(parents=True, exist_ok=True)

    valid = ~np.ma.getmaskarray(Z)
    raster_df = pd.DataFrame(
        {"lon": LON[valid], "lat": LAT[valid], "log_w_ratio": Z[valid].data}
    )
    raster_path = args.out / "genetic_surface_raster.csv"
    raster_df.to_csv(raster_path, index=False)
    print(f"[fit] wrote {raster_path} ({raster_df.shape})")

    edge_arr = np.array(sp_graph.edges)
    edgew_df = pd.DataFrame(
        {
            "node_i": edge_arr[:, 0],
            "node_j": edge_arr[:, 1],
            "w": sp_graph.w,
            "log10_w_ratio": log_ratio,
        }
    )
    edgew_path = args.out / "edgew_genetic.csv"
    edgew_df.to_csv(edgew_path, index=False)
    print(f"[fit] wrote {edgew_path} ({edgew_df.shape})")

    nodepos_df = pd.DataFrame(
        {
            "lon": sp_graph.node_pos[:, 0],
            "lat": sp_graph.node_pos[:, 1],
            "n_samples": [
                sp_graph.nodes[n]["n_samples"] for n in range(len(sp_graph.nodes))
            ],
        }
    )
    nodepos_path = args.out / "nodepos_genetic.csv"
    nodepos_df.to_csv(nodepos_path, index=False)
    print(f"[fit] wrote {nodepos_path} ({nodepos_df.shape})")

    meta = {
        "grid_res": args.grid_res,
        "habitat_source": f"Natural Earth coastline ({geojson_path.name}, name='{args.country}')",
        "coast_buf": args.coast_buf,
        "n_nodes": int(grid.shape[0]),
        "n_edges": int(edges.shape[0]),
        "n_observed_nodes": int(sp_graph.n_observed_nodes),
        "scale_snps": True,
        "optimize_q": None,
        "lamb_cv": lamb_cv,
        "bfile": str(args.bfile),
        "n_samples": len(coord),
        "n_populations": int(n_populations),
        "n_snps": int(genotypes.shape[1]),
        "n_snps_total": int(n_snps_total),
    }
    meta_path = args.out / "genetic_feems_meta.json"
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[fit] wrote {meta_path}")
    print(json.dumps(meta, indent=2))

    if args.permute:
        run_permutations(args, sp_graph, genotypes, coord, grid, edges, outer, LON, LAT, Z, lamb_cv)


# --------------------------------------------------------------------------- #


def check_deps(stage: str):
    needed = ["numpy", "pandas"]
    if stage in ("fit", "all"):
        needed += [
            "pandas_plink",
            "sklearn",
            "shapely",
            "networkx",
            "feems",
            "matplotlib",
        ]
    missing = []
    for mod in needed:
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    if missing:
        sys.exit(
            f"error: missing required package(s): {', '.join(missing)}\n"
            f"Install with: conda env create -f genetic_feems_env.yml"
        )


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--bfile", required=True, help="PLINK prefix (no .bed/.bim/.fam extension)"
    )
    p.add_argument(
        "--coords",
        required=True,
        help="CSV with a 'population', 'latitude', 'longitude' column (e.g. philippine_ethnolinguistic_coords.csv)",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("./out"),
        help="output directory (default: ./out)",
    )
    p.add_argument("--stage", choices=["match", "fit", "all"], default="all")
    p.add_argument(
        "--id-prefix",
        default="PHIL_",
        help="group prefix stripped from FID/IID before matching against 'population'; the population is then read as the first underscore-delimited token of the remainder, e.g. PHIL_AgtaBulusan_BUL-02 -> AgtaBulusan (default: 'PHIL_')",
    )
    p.add_argument(
        "--alias",
        action="append",
        default=None,
        help="extra 'Old=New' population rename for .fam tokens that don't match the coords table verbatim (e.g. 'ManoboRajahKabunsuwan=ManoboRK', already built in); repeatable",
    )

    p.add_argument(
        "--geojson",
        default=None,
        help="path to countries.geojson (default: next to this script, else the installed feems package copy)",
    )
    p.add_argument(
        "--country",
        default="Philippines",
        help="properties.name to select from the geojson (default: Philippines)",
    )
    p.add_argument(
        "--coast-buf",
        type=float,
        default=0.8,
        help="habitat buffer in degrees (default 0.8 — 0.7 strands the Sulu Archipelago)",
    )
    p.add_argument(
        "--grid",
        default=None,
        help="path to a FEEMS grid shapefile (default: resolved from the installed feems package)",
    )
    p.add_argument(
        "--grid-res",
        default="grid_100.shp",
        help="grid filename within the feems package data dir (default: grid_100.shp)",
    )

    p.add_argument(
        "--n-folds",
        type=int,
        default=0,
        help="CV folds; 0 = leave-one-deme-out (default)",
    )
    p.add_argument("--lamb-min", type=float, default=1e-6)
    p.add_argument("--lamb-max", type=float, default=1e2)
    p.add_argument("--lamb-n", type=int, default=20)
    p.add_argument(
        "--raster-res",
        type=int,
        default=300,
        help="raster grid points along longitude (default 300)",
    )

    p.add_argument(
        "--max-snps",
        type=int,
        default=100000,
        help="randomly subsample to this many SNPs before fitting (default: 100000; "
        "pass 0 to use all). 100000 is what produced the recorded run in "
        "genetic_feems_meta.json (2246051 total -> 100000 subsampled -> 69833 after "
        "dropping invariant SNPs). FEEMS collapses genotypes to per-deme allele "
        "frequencies, so beyond ~1e5 SNPs the surface changes little while the "
        "per-permutation cost keeps scaling",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=42,
        help="RNG seed for --max-snps subsampling; with the same --bfile and --seed "
        "the subsample is reproducible across runs (default 42)",
    )
    p.add_argument(
        "--plot",
        action="store_true",
        help="also write feems_cv_curve.png / feems_surface.png (plain matplotlib, no cartopy)",
    )

    p.add_argument(
        "--lamb",
        type=float,
        default=None,
        help="skip CV and fit at this lambda directly — lets --permute be re-run "
        "later without repeating the expensive CV sweep (e.g. --lamb <lamb_cv from "
        "a prior genetic_feems_meta.json>)",
    )
    p.add_argument(
        "--permute",
        action="store_true",
        help="after the real fit, run a geographic-shuffle null model (see --n-perm)",
    )
    p.add_argument(
        "--n-perm", type=int, default=100, help="number of null permutations (default 100)"
    )
    p.add_argument(
        "--perm-seed",
        type=int,
        default=0,
        help="base RNG seed for the permutation loop — separate from --seed, which "
        "only controls --max-snps subsampling (default 0)",
    )
    p.add_argument(
        "--perm-restart",
        action="store_true",
        help="ignore any existing genetic_null_checkpoint.npz and restart the permutation loop from 0",
    )
    p.add_argument(
        "--perm-cv",
        action="store_true",
        help="re-run cross-validation inside every permutation so each gets its own "
        "lambda, and use the profiled CV minimum (T_cv) as the test statistic. "
        "Removes lambda from the statistic — without this, the fixed-lambda T_het "
        "statistic inherits lambda's shrinkage and is only comparable at that lambda. "
        "Costs one CV sweep per permutation.",
    )
    return p.parse_args()


def main():
    args = parse_args()
    check_deps(args.stage)

    match_df = None
    if args.stage in ("match", "all"):
        match_df = run_match(args)
    if args.stage in ("fit", "all"):
        if match_df is None:
            match_path = args.out / "GENETIC_feems_sample_match.csv"
            if not match_path.exists():
                sys.exit(
                    f"error: --stage fit needs {match_path} — run --stage match first"
                )
            match_df = pd.read_csv(match_path)
            match_df["row"] = np.arange(len(match_df))
        run_fit(args, match_df)


if __name__ == "__main__":
    main()
