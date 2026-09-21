# =============================================================================
# [1] Correlate Analysis — spatial KDE surfaces + pairwise linear regression
#
# Covers all six pairs among the four domains (phoneme, grammar, cognate,
# genetic). Where [0]_CORRELATE_admixture.R pairs by language identity and is
# therefore limited to each pair's intersection (n = 18-48), this pairs by
# GEOGRAPHIC LOCATION: every domain's scattered points are smoothed onto one
# shared grid over the Philippines, then each pair is regressed cell by cell.
# That uses every language in every domain and equalises point counts by
# construction, since all four are evaluated on the same grid.
#
# Surfaces are Nadaraya-Watson kernel means, NOT weighted KDEs:
#
#     m(s) = sum_i v_i K_h(s - s_i) / sum_i K_h(s - s_i)
#
# The denominator is the plain (unweighted) KDE — the local sampling density.
# Dividing by it turns a weighted sum into a weighted mean, so the surface
# estimates the value field rather than "where was this domain sampled". That
# matters here because sampling effort is wildly uneven: genetic contributes 75
# points and grammar 26, and without the normalisation grammar's surface would
# sit low everywhere purely for being sparse.
#
# INTERPRETATION CAVEAT: grid cells ~4 km apart smoothed at a 100 km bandwidth
# are near-perfectly autocorrelated, so grid_p_raw is ~0 for essentially every
# pair regardless of the truth — it is reported only to make that inflation
# visible. Read grid_p_dutilleul, which corrects for spatial autocorrelation via
# an effective sample size, and paired_p, which is the entity-level check.
#
# Inputs:  data/feems/phoneme_surface_raster.csv   (the shared grid + landmask)
#          data/network_distance/COGNATE_final.csv
#          data/cosine_distribution/PHONEME_cossim_marked.csv, data/phoneme/RUHLENdf_PH.csv
#          data/cosine_distribution/GRAMMAR_cossim_marked.csv, data/grammar/GRAMBANKdf_full.csv
#          data/network_distance/GENETIC_final.csv, data/genetic/GENETIC_subgroup_lookup.csv
# Outputs: data/correlate/CORRELATE_kde_surfaces.csv  (grid + 4 surfaces + masks)
#          data/correlate/CORRELATE_kde_pairwise.csv  (6 rows, the results table)
#          figures/correlate/kde_surfaces.png
#          figures/correlate/kde_<x>_vs_<y>.png  (6)
#          figures/correlate/kde_pairwise_grid.png
# =============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(geosphere)
library(lingtypology)
library(SpatialPack)

# `maps` is deliberately NOT attached: ggplot2::map_data() reaches into it
# without it being on the search path, and attaching it would mask purrr::map
# with maps::map (the same hazard [2]_*_cosine_distribution_analysis.R notes for
# mclust). purrr::map is namespaced below regardless.

dir.create(here("data", "correlate"), showWarnings = FALSE)
dir.create(here("figures", "correlate"), showWarnings = FALSE)

# Shared across domains so differences between surfaces reflect the data rather
# than differing amounts of smoothing.
BANDWIDTH_KM <- 100

# 2 bandwidths. Past this the Gaussian weight is exp(-2) = 0.135 and the local
# mean is extrapolation toward the domain's global mean, not estimation — it
# would still return a confident-looking number, so those cells are masked.
SUPPORT_MAX_KM <- 200

# Dutilleul's correction needs pairwise distances among cells; 54,804 cells is
# 3e9 pairs, so the corrected test runs on a systematically thinned subset.
DUTILLEUL_CELLS <- 2000

# Matches SPAN_SE_MAX in [0]_CORRELATE_admixture.R — a population whose span_se
# exceeds the largest observed span_admx cannot be measured at any true value.
SPAN_SE_MAX <- 0.083

# Plot-only thinning; every regression below uses the full supported set.
PLOT_CELLS <- 8000

DOMAIN_ORDER <- c("phoneme", "grammar", "cognate", "genetic")
DOMAIN_LABELS <- c(phoneme = "phoneme: cossim_span - cossim_unr",
                   grammar = "grammar: cossim_span - cossim_unr",
                   cognate = "cognate: Spanish loanwords",
                   genetic = "genetic: span_admx")


# ── 1. The shared grid ───────────────────────────────────────────────────────
# phoneme/genetic/cognate FEEMS rasters are bit-identical in lon/lat: a 300x400
# lattice already point-in-polygon masked to the buffered Philippine coastline
# by python/genetic_feems.py. Reused here purely as grid + landmask, so
# log_w_ratio (the FEEMS migration surface) is discarded.
grid <- read_csv(here("data", "feems", "phoneme_surface_raster.csv"),
                 show_col_types = FALSE) |>
  select(lon, lat)

grid_xy <- as.matrix(grid[, c("lon", "lat")])
message(nrow(grid), " grid cells (shared FEEMS habitat lattice).")


# ── 2. Domain point frames ───────────────────────────────────────────────────
# glottocode is carried for the paired half; bridging per domain follows
# [0]_CORRELATE_admixture.R §3-5.
COGNATE_final <- read_csv(here("data", "network_distance", "COGNATE_final.csv"),
                          show_col_types = FALSE)
cognate_pts <- COGNATE_final |>
  transmute(glottocode, longitude, latitude, value = number_of_loans)

GRAMMAR_cossim  <- read_csv(here("data", "cosine_distribution", "GRAMMAR_cossim_marked.csv"),
                            show_col_types = FALSE)
GRAMBANKdf_full <- read_csv(here("data", "grammar", "GRAMBANKdf_full.csv"),
                            show_col_types = FALSE)
grammar_pts <- GRAMMAR_cossim |>
  left_join(GRAMBANKdf_full |> select(language, glottocode = Language_ID) |>
              distinct(language, .keep_all = TRUE), by = "language") |>
  transmute(glottocode, longitude, latitude, value = cossim_span - cossim_unr)

PHONEME_cossim <- read_csv(here("data", "cosine_distribution", "PHONEME_cossim_marked.csv"),
                           show_col_types = FALSE)
RUHLENdf_PH    <- read_csv(here("data", "phoneme", "RUHLENdf_PH.csv"),
                           show_col_types = FALSE)
phoneme_pts <- PHONEME_cossim |>
  left_join(RUHLENdf_PH |> select(language, iso6393) |>
              distinct(language, .keep_all = TRUE), by = "language") |>
  transmute(glottocode = gltc.iso(iso6393), longitude, latitude,
            value = cossim_span - cossim_unr)

GENETIC_final  <- read_csv(here("data", "network_distance", "GENETIC_final.csv"),
                           show_col_types = FALSE)
GENETIC_lookup <- read_csv(here("data", "genetic", "GENETIC_subgroup_lookup.csv"),
                           show_col_types = FALSE)
genetic_pts <- GENETIC_final |>
  left_join(GENETIC_lookup |> select(population, glottocode), by = "population") |>
  filter(!is.na(span_admx), span_se < SPAN_SE_MAX) |>
  transmute(glottocode, longitude, latitude, value = span_admx)

domains <- list(phoneme = phoneme_pts, grammar = grammar_pts,
                cognate = cognate_pts, genetic = genetic_pts)

stopifnot("a domain frame has NA coordinates or values" =
            all(purrr::map_lgl(domains, \(d) !anyNA(d$longitude) && !anyNA(d$latitude) &&
                                      !anyNA(d$value))))

for (nm in DOMAIN_ORDER) {
  message(nm, ": ", nrow(domains[[nm]]), " points | value range [",
          round(min(domains[[nm]]$value), 4), ", ",
          round(max(domains[[nm]]$value), 4), "]")
}


# ── 3. Nadaraya-Watson surfaces ──────────────────────────────────────────────
# One domain at a time: the cell-to-point distance matrix is ~50 MB and is not
# retained past the two summaries taken from it.
nw_surface <- function(pts, grid_xy) {
  d <- distm(grid_xy, as.matrix(pts[, c("longitude", "latitude")]),
             fun = distHaversine) / 1000                     # km
  w <- exp(-0.5 * (d / BANDWIDTH_KM)^2)                      # gaussian kernel
  list(
    value     = as.vector((w %*% pts$value) / rowSums(w)),   # weighted mean
    supported = as.vector(apply(d, 1, min)) <= SUPPORT_MAX_KM
  )
}

surfaces <- grid
for (nm in DOMAIN_ORDER) {
  s <- nw_surface(domains[[nm]], grid_xy)
  surfaces[[paste0(nm, "_value")]]     <- s$value
  surfaces[[paste0(nm, "_supported")]] <- s$supported
  message(nm, ": ", sum(s$supported), " of ", nrow(grid), " cells supported (",
          round(100 * mean(s$supported)), "%).")
}

write.csv(surfaces, here("data", "correlate", "CORRELATE_kde_surfaces.csv"),
          row.names = FALSE)


# ── 4. Pairwise regressions ──────────────────────────────────────────────────
# combn over DOMAIN_ORDER gives the 6 unique pairs deterministically; the first
# domain is x, the second y. Same idiom as R/shared/pairwise_network_distance.R.
pairs <- t(combn(DOMAIN_ORDER, 2)) |>
  as.data.frame(stringsAsFactors = FALSE) |>
  setNames(c("domain_x", "domain_y")) |>
  as_tibble()

fit_pair <- function(dx, dy) {
  # --- grid regression, restricted to cells where BOTH domains are supported --
  sub <- surfaces |>
    filter(.data[[paste0(dx, "_supported")]], .data[[paste0(dy, "_supported")]]) |>
    transmute(lon, lat,
              x = .data[[paste0(dx, "_value")]],
              y = .data[[paste0(dy, "_value")]])

  m_grid <- lm(y ~ x, data = sub)
  cf     <- summary(m_grid)$coefficients

  # Dutilleul's modified t-test on a systematically thinned subset: corrects the
  # test for spatial autocorrelation by estimating an effective sample size.
  idx  <- unique(round(seq(1, nrow(sub), length.out = min(DUTILLEUL_CELLS, nrow(sub)))))
  thin <- sub[idx, ]
  p_dut <- tryCatch(
    modified.ttest(thin$x, thin$y, coords = as.matrix(thin[, c("lon", "lat")]))$p.value,
    error = function(e) NA_real_
  )

  # --- paired regression on shared glottocodes (the entity-level check) -------
  paired <- domains[[dx]] |>
    select(glottocode, x = value) |>
    inner_join(domains[[dy]] |> select(glottocode, y = value),
               by = "glottocode", relationship = "many-to-many")

  m_pair <- lm(y ~ x, data = paired)
  cfp    <- summary(m_pair)$coefficients

  tibble(
    domain_x         = dx,
    domain_y         = dy,
    grid_n           = nrow(sub),
    grid_slope       = cf["x", "Estimate"],
    grid_r2          = summary(m_grid)$r.squared,
    grid_p_raw       = cf["x", "Pr(>|t|)"],
    grid_p_dutilleul = p_dut,
    paired_n         = nrow(paired),
    paired_slope     = cfp["x", "Estimate"],
    paired_r2        = summary(m_pair)$r.squared,
    paired_p         = cfp["x", "Pr(>|t|)"]
  )
}

CORRELATE_kde_pairwise <- purrr::map2(pairs$domain_x, pairs$domain_y, fit_pair) |>
  list_rbind()

print(CORRELATE_kde_pairwise, width = Inf)

write.csv(CORRELATE_kde_pairwise,
          here("data", "correlate", "CORRELATE_kde_pairwise.csv"),
          row.names = FALSE)


# ── 5. Surface maps ──────────────────────────────────────────────────────────
# Map idiom follows the [6]_feems_plot_*.R house style: geom_tile raster, maps
# coastline drawn unfilled on top, coord_fixed, axes/grid/ticks dropped because
# the coastline is the only reference frame the reader needs.
map_subset <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

surface_map <- function(nm) {
  d <- surfaces |>
    filter(.data[[paste0(nm, "_supported")]]) |>
    transmute(lon, lat, value = .data[[paste0(nm, "_value")]])

  ggplot() +
    geom_tile(data = d, aes(x = lon, y = lat, fill = value)) +
    geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
                 fill = NA, colour = "black", linewidth = 0.25) +
    geom_point(data = domains[[nm]], aes(x = longitude, y = latitude),
               size = 0.8, colour = "grey15") +
    scale_fill_gradient(low = "white", high = "navy", name = NULL) +
    coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
    labs(title = DOMAIN_LABELS[[nm]],
         subtitle = paste0(nrow(domains[[nm]]), " points | ",
                           round(100 * mean(surfaces[[paste0(nm, "_supported")]])),
                           "% of grid supported")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8, colour = "grey35"),
          legend.key.width = unit(0.3, "cm"))
}

p_surfaces <- wrap_plots(purrr::map(DOMAIN_ORDER, surface_map), nrow = 2)
print(p_surfaces)
ggsave(here("figures", "correlate", "kde_surfaces.png"), p_surfaces,
       width = 10, height = 11, units = "in", dpi = 300, bg = "white")


# ── 6. Pair scatter panels ───────────────────────────────────────────────────
# Grid cells (grey) and the real paired points (accent) share one set of axes.
# Solid line = regression across all supported grid cells; dashed = regression
# on the paired points alone. No legend: the subtitle names both lines.
pair_panel <- function(row) {
  dx <- row$domain_x; dy <- row$domain_y

  sub <- surfaces |>
    filter(.data[[paste0(dx, "_supported")]], .data[[paste0(dy, "_supported")]]) |>
    transmute(x = .data[[paste0(dx, "_value")]], y = .data[[paste0(dy, "_value")]])

  paired <- domains[[dx]] |>
    select(glottocode, x = value) |>
    inner_join(domains[[dy]] |> select(glottocode, y = value),
               by = "glottocode", relationship = "many-to-many")

  m_grid <- lm(y ~ x, data = sub)
  m_pair <- lm(y ~ x, data = paired)

  plot_sub <- sub[unique(round(seq(1, nrow(sub),
                                   length.out = min(PLOT_CELLS, nrow(sub))))), ]

  fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)

  ggplot() +
    geom_point(data = plot_sub, aes(x, y), colour = "grey72", size = 0.35, alpha = 0.5) +
    geom_point(data = paired, aes(x, y), colour = "#2ca6a4", size = 1.6, alpha = 0.9) +
    geom_abline(intercept = coef(m_grid)[1], slope = coef(m_grid)[2],
                colour = "grey25", linewidth = 0.8) +
    geom_abline(intercept = coef(m_pair)[1], slope = coef(m_pair)[2],
                colour = "#2ca6a4", linewidth = 0.8, linetype = "dashed") +
    labs(title = paste0(dy, " ~ ", dx),
         subtitle = sprintf(
           "solid grid: r2=%.3f, Dutilleul p=%s  |  dashed paired (n=%d): r2=%.3f, p=%s",
           summary(m_grid)$r.squared, fmt_p(row$grid_p_dutilleul),
           nrow(paired), summary(m_pair)$r.squared, fmt_p(row$paired_p)),
         x = DOMAIN_LABELS[[dx]], y = DOMAIN_LABELS[[dy]]) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 7, colour = "grey35"),
          axis.title = element_text(size = 8))
}

panels <- purrr::map(seq_len(nrow(CORRELATE_kde_pairwise)),
              \(i) pair_panel(CORRELATE_kde_pairwise[i, ]))

for (i in seq_along(panels)) {
  r <- CORRELATE_kde_pairwise[i, ]
  ggsave(here("figures", "correlate",
              sprintf("kde_%s_vs_%s.png", r$domain_x, r$domain_y)),
         panels[[i]], width = 6, height = 4.5, units = "in", dpi = 300, bg = "white")
}

p_grid <- wrap_plots(panels, nrow = 2)
print(p_grid)
ggsave(here("figures", "correlate", "kde_pairwise_grid.png"), p_grid,
       width = 16, height = 9, units = "in", dpi = 300, bg = "white")
