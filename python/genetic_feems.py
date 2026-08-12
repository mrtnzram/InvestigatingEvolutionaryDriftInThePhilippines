#!/usr/bin/env python3
"""
genetic_feems.py — FEEMS effective-migration surface estimation on Philippine
genetic PLINK data, run headless on a remote server over SSH.

Inputs:  --bfile PLINK prefix (bed/bim/fam, already QC'd/pruned — this script never
         shells out to a PLINK binary), --coords population coordinates CSV.
Outputs: genetic_surface_raster.csv, nodepos_genetic.csv, edgew_genetic.csv,
         genetic_feems_meta.json, GENETIC_feems_sample_match.csv.

Usage:
    python genetic_feems.py --stage match --bfile /path/to/Phil_1000g_SGDP_1.66M \
        --coords philippine_ethnolinguistic_coords.csv --out ./out

    python genetic_feems.py --stage fit --bfile /path/to/Phil_1000g_SGDP_1.66M \
        --coords philippine_ethnolinguistic_coords.csv --out ./out

    python genetic_feems.py --stage all ...   # match, then fit, in one call (default)

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

    n_folds = args.n_folds if args.n_folds > 0 else sp_graph.n_observed_nodes
    lamb_grid = np.geomspace(args.lamb_min, args.lamb_max, args.lamb_n)[::-1]
    print(f"[fit] running CV: {args.lamb_n} lambdas x {n_folds} folds")
    cv_err = run_cv(sp_graph, lamb_grid, n_folds=n_folds, factr=1e10)
    mean_cv_err = np.mean(cv_err, axis=0)
    argmin = int(np.argmin(mean_cv_err))
    lamb_cv = float(lamb_grid[argmin])
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
        fig.savefig(args.out / "feems_cv_curve.png", bbox_inches="tight")
        plt.close(fig)

    sp_graph.fit(lamb=lamb_cv, optimize_q=None)
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
        fig.savefig(args.out / "feems_surface.png", bbox_inches="tight")
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
        default=None,
        help="randomly subsample to this many SNPs before fitting (default: use all)",
    )
    p.add_argument(
        "--seed", type=int, default=42, help="RNG seed for --max-snps subsampling"
    )
    p.add_argument(
        "--plot",
        action="store_true",
        help="also write feems_cv_curve.png / feems_surface.png (plain matplotlib, no cartopy)",
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
