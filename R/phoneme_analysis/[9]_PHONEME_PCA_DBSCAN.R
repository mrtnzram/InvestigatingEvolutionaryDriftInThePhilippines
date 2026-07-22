# =============================================================================
# [9] Phoneme Analysis — PCA + DBSCAN clustering
#
# Runs PCA directly on the IDF-weighted phoneme feature matrix (the matrix that
# feeds [1]'s cosine similarity, not the cosine matrix itself), then DBSCAN-
# clusters languages in PC1+PC2 space — a genealogy-agnostic grouping, useful
# to compare against the ABVD-tree-based family subgroups from [0].
#
#   Part A — Philippine languages only (58).
#   Part B — exploratory: full database (Philippine + interest + unrelated, 198).
#
# DEPENDENCY / RUN ORDER:
#   - Duplicates [1]'s IDF-weighting logic inline (recomputed here rather than
#     sourced, since [1] never exposes its weighted matrix — only the collapsed
#     cosine similarity escapes that function). If [1]'s weighting scheme ever
#     changes, this script needs revisiting too.
#   - Requires data/PHONEME_driver_table.csv to already exist (run
#     [1.5]_PHONEME_drivers.R at least once) — used to translate PCA loadings
#     from raw phoneme_XXX ids to real IPA symbols.
#   - Requires data/PHONEME_subgroup_lookup.csv (from [0]_Phylogenetic_Tree.R)
#     for the family-coloured plot and for deriving DBSCAN's minPts.
#   - source()s [3]_PHONEME_network_distance.R for the waypoint-map deliverable
#     (same pattern [8] uses).
#
# Input:   data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen_austronesian.csv,
#          data/PHONEME_driver_table.csv, data/PHONEME_subgroup_lookup.csv,
#          data/nodes.csv, data/edges.csv (via [3])
# Outputs: data/PHONEME_pca_dbscan_{ph,full}.csv, data/PHONEME_pca_loadings_{ph,full}.csv,
#          figures/phoneme/pca/PHONEME_pca_dbscan_{ph,full}.png,
#          figures/phoneme/pca/PHONEME_pca_subgroup_ph.png,
#          figures/phoneme/pca/PHONEME_pca_family_full.png,
#          figures/phoneme/pca/PHONEME_pca_loadings_{ph,full}.png,
#          figures/phoneme/pca/PHONEME_kNNdist_{ph,full}.png (diagnostic),
#          figures/phoneme/mst_waypoints/PHONEME_network_by_dbscan_{ph,full}.png
# =============================================================================

library(tidyverse)
library(here)
library(dbscan)
library(sf)
library(patchwork)

dir.create(here("figures", "phoneme", "pca"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Rebuild the IDF-weighted phoneme matrix (duplicated from [1] — see header)
# =============================================================================
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"), show_col_types = FALSE)

phoneme_cols <- RUHLENdf %>%
  dplyr::select(-language, -source, -iso6393, -Language_type, -latitude, -longitude) %>%
  colnames()

phoneme_freq <- read_csv(here("data", "phoneme_freq_ruhlen_austronesian.csv"), show_col_types = FALSE)

aligned_freq <- phoneme_freq %>%
  dplyr::filter(phoneme %in% phoneme_cols) %>%
  dplyr::arrange(match(phoneme, phoneme_cols))
idf_weights <- aligned_freq$IDF

binary_data <- RUHLENdf %>% dplyr::select(dplyr::all_of(phoneme_cols)) %>% as.matrix()
weighted_data <- sweep(binary_data, 2, idf_weights, FUN = "*")
rownames(weighted_data) <- RUHLENdf$language   # [1] leaves this NULL; PCA/plots need it here

# =============================================================================
# 2. Helper functions
# =============================================================================

# Drop columns with zero variance in the given row subset. Harmless numerically
# with scale.=FALSE, but keeps the matrix lean/auditable and guards against a
# future switch to scale.=TRUE (which WOULD error on a constant column).
drop_zero_variance <- function(mat, subset_name = "") {
  col_var <- apply(mat, 2, var)
  keep <- col_var > 0
  if (sum(!keep) > 0) {
    message(subset_name, ": dropping ", sum(!keep), " of ", ncol(mat),
            " zero-variance columns for this subset.")
  }
  mat[, keep, drop = FALSE]
}

# Simplified "kneedle" knee detection on a sorted, monotonically increasing
# curve: normalize to [0,1]x[0,1], return the index with max perpendicular
# distance from the chord connecting the first and last points.
find_knee <- function(y) {
  x  <- seq_along(y)
  xn <- (x - min(x)) / (max(x) - min(x))
  yn <- (y - min(y)) / (max(y) - min(y))
  x1 <- xn[1]; y1 <- yn[1]; x2 <- xn[length(xn)]; y2 <- yn[length(yn)]
  d <- abs((y2 - y1) * xn - (x2 - x1) * yn + x2 * y1 - y2 * x1) /
    sqrt((y2 - y1)^2 + (x2 - x1)^2)
  which.max(d)
}

# Choose eps via the k-NN distance elbow for a (PC1,PC2) score matrix, given
# minPts; optionally saves a diagnostic plot showing the chosen point.
choose_eps <- function(scores_mat, minPts, plot_path = NULL, title = "") {
  kdist  <- sort(dbscan::kNNdist(scores_mat, k = minPts))
  knee_i <- find_knee(kdist)
  eps    <- kdist[knee_i]
  if (!is.null(plot_path)) {
    df <- tibble(idx = seq_along(kdist), kdist = kdist)
    p <- ggplot(df, aes(idx, kdist)) +
      geom_line() +
      geom_point(data = df[knee_i, ], color = "firebrick", size = 2) +
      geom_hline(yintercept = eps, linetype = "dashed", color = "firebrick") +
      labs(title = title,
           subtitle = sprintf("minPts = %d | chosen eps = %.4f", minPts, eps),
           x = "Points sorted by k-NN distance", y = paste0(minPts, "-NN distance")) +
      theme_minimal()
    ggsave(plot_path, p, width = 5, height = 4, units = "in", dpi = 300)
  }
  eps
}

# minPts from the real linguistic subgroup sizes: smallest non-singleton size
# (a subgroup of 1 shouldn't define "how big a cluster must be" — DBSCAN's own
# noise category already handles isolates).
compute_minpts <- function(subgroup_lookup) {
  sizes <- table(subgroup_lookup$subgroup)
  as.integer(min(sizes[sizes > 1]))
}

# Cluster palette (same construction as [0]'s subgroup_pal), sized to the
# number of real (non-noise) clusters, plus a fixed grey for "Noise".
soften <- function(cols, s_mult = 0.55, v_mult = 0.95) {
  h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  grDevices::hsv(h["h", ], h["s", ] * s_mult, pmin(h["v", ] * v_mult, 1))
}
build_cluster_pal <- function(cluster_levels) {
  real_levels <- setdiff(cluster_levels, "Noise")
  .poly <- grDevices::palette.colors(NULL, "Polychrome 36")
  .lum  <- colSums(grDevices::col2rgb(.poly) * c(0.299, 0.587, 0.114))
  pal <- setNames(soften(.poly[.lum < 200])[seq_along(real_levels)], real_levels)
  c(pal, "Noise" = "grey60")
}

# Run PCA (center=TRUE, scale.=FALSE — see header rationale) + DBSCAN (PC1+PC2
# only) on a language x feature matrix. Returns per-language scores, cluster
# labels, variance-explained, and the fitted objects for downstream loadings.
run_pca_dbscan <- function(mat, minPts, subset_name, eps_plot_path = NULL) {
  mat <- drop_zero_variance(mat, subset_name)
  pca <- prcomp(mat, center = TRUE, scale. = FALSE)
  var_pct <- summary(pca)$importance[2, 1:2] * 100

  scores <- as_tibble(pca$x[, 1:2, drop = FALSE], rownames = "language")

  eps <- choose_eps(as.matrix(scores[, c("PC1", "PC2")]), minPts,
                    plot_path = eps_plot_path, title = subset_name)
  db  <- dbscan::dbscan(as.matrix(scores[, c("PC1", "PC2")]), eps = eps, minPts = minPts)

  cluster_chr    <- ifelse(db$cluster == 0, "Noise", as.character(db$cluster))
  numeric_levels <- sort(as.integer(setdiff(unique(cluster_chr), "Noise")))
  cluster_levels <- c(as.character(numeric_levels), "Noise")

  scores <- scores %>% mutate(cluster = factor(cluster_chr, levels = cluster_levels))

  message(subset_name, ": minPts=", minPts, ", eps=", round(eps, 4), " -> ",
          length(numeric_levels), " cluster(s), ",
          sum(cluster_chr == "Noise"), " noise point(s) of ", nrow(scores))

  list(pca = pca, scores = scores, var_pct = var_pct, minPts = minPts, eps = eps)
}

# Top-|loading| features per PC, joined to IPA/class where available.
get_loadings <- function(pca, driver_table) {
  rot <- pca$rotation[, 1:2, drop = FALSE]
  as_tibble(rot, rownames = "phoneme") %>%
    rename(PC1_loading = PC1, PC2_loading = PC2) %>%
    left_join(driver_table %>% dplyr::select(phoneme, ipa, class, IDF), by = "phoneme") %>%
    mutate(label = coalesce(ipa, phoneme)) %>%
    arrange(desc(abs(PC1_loading)))
}

top_loadings_plot <- function(loadings_df, pc_col, title) {
  d <- loadings_df %>%
    mutate(abs_loading = abs(.data[[pc_col]])) %>%
    slice_max(abs_loading, n = 15) %>%
    mutate(label = fct_reorder(label, .data[[pc_col]]))
  ggplot(d, aes(x = .data[[pc_col]], y = label, fill = .data[[pc_col]] > 0)) +
    geom_col(show.legend = FALSE) +
    scale_fill_manual(values = c(`TRUE` = "#2ca6a4", `FALSE` = "#e07a5f")) +
    labs(title = title, x = "Loading", y = NULL) +
    theme_minimal()
}

# =============================================================================
# 3. Network objects from [3] (for the map deliverable, both Parts)
# =============================================================================
source(here("R", "phoneme_analysis", "[3]_PHONEME_network_distance.R"), echo = FALSE)

subgroup_lookup <- read_csv(here("data", "PHONEME_subgroup_lookup.csv"), show_col_types = FALSE)
subgroup_pal    <- subgroup_lookup %>% distinct(subgroup, colour) %>% deframe()
driver_table    <- read_csv(here("data", "PHONEME_driver_table.csv"), show_col_types = FALSE)

minPts <- compute_minpts(subgroup_lookup)
message("Computed minPts = ", minPts, " (smallest non-singleton family-subgroup size).")

plot_network_by_cluster <- function(points_df, cluster_pal, refdf1 = MANILA,
                                    lon_range = c(116, 127), lat_range = c(4, 21),
                                    title = NULL) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70") +
    geom_sf(data = full_tree_sf %>% dplyr::filter(source == "main"),
            linewidth = 1, color = "grey40") +
    geom_sf(data = arrow_connectors,
            color = "grey60", linewidth = 0.5,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_point(data = points_df,
               aes(x = longitude, y = latitude, fill = cluster),
               size = 3, shape = 21, colour = "grey20") +
    geom_point(data = refdf1, aes(x = longitude, y = latitude),
               shape = 23, size = 3, fill = "white", stroke = 1) +
    scale_fill_manual(values = cluster_pal, guide = "none") +
    coord_sf(xlim = lon_range, ylim = lat_range) +
    theme_minimal() +
    labs(title = title, x = "Longitude", y = "Latitude")
}

build_cluster_legend_strip <- function(points_df, cluster_pal) {
  present <- as.character(unique(points_df$cluster))
  cluster_lat <- points_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(mean_lat = mean(latitude), .groups = "drop")
  legend_df <- tibble(cluster = names(cluster_pal)) %>%
    dplyr::filter(cluster %in% present) %>%
    dplyr::left_join(cluster_lat, by = "cluster") %>%
    dplyr::mutate(is_noise = cluster == "Noise") %>%
    dplyr::arrange(is_noise, dplyr::desc(mean_lat)) %>%
    dplyr::mutate(row = dplyr::row_number())

  ggplot(legend_df) +
    geom_tile(aes(x = 0, y = -row, fill = cluster), width = 0.9, height = 0.8) +
    geom_text(aes(x = 0.65, y = -row, label = cluster), hjust = 0, size = 2.8, colour = "grey15") +
    scale_fill_manual(values = cluster_pal, guide = "none") +
    scale_x_continuous(limits = c(-0.5, 6.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(6, 2, 6, 2))
}

# =============================================================================
# PART A — Philippine languages only
# =============================================================================
ph_mat <- weighted_data[RUHLENdf$Language_type == "Philippine Language", , drop = FALSE]

res_ph <- run_pca_dbscan(ph_mat, minPts, "Phoneme PH-only",
                         eps_plot_path = here("figures", "phoneme", "pca", "PHONEME_kNNdist_ph.png"))
scores_ph   <- res_ph$scores %>% mutate(Language_type = "Philippine Language", .after = language)
loadings_ph <- get_loadings(res_ph$pca, driver_table)

write_csv(scores_ph, here("data", "PHONEME_pca_dbscan_ph.csv"))
write_csv(loadings_ph, here("data", "PHONEME_pca_loadings_ph.csv"))

cluster_pal_ph <- build_cluster_pal(levels(scores_ph$cluster))

# (a) PCA scatter, coloured by DBSCAN cluster
p_pca_dbscan_ph <- ggplot(scores_ph, aes(PC1, PC2, colour = cluster)) +
  geom_point(size = 3) +
  scale_colour_manual(values = cluster_pal_ph, name = "DBSCAN cluster") +
  labs(title = "Phoneme PCA (Philippine languages) — coloured by DBSCAN cluster",
       x = sprintf("PC1 (%.1f%% var.)", res_ph$var_pct[1]),
       y = sprintf("PC2 (%.1f%% var.)", res_ph$var_pct[2])) +
  theme_minimal()
print(p_pca_dbscan_ph)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_dbscan_ph.png"), p_pca_dbscan_ph,
       width = 7, height = 5.5, units = "in", dpi = 300)

# (b) PCA scatter, coloured by language family (all 58 languages have a subgroup)
scores_ph_fam <- scores_ph %>%
  inner_join(subgroup_lookup %>% dplyr::select(language, subgroup), by = "language")

p_pca_subgroup_ph <- ggplot(scores_ph_fam, aes(PC1, PC2, colour = subgroup)) +
  geom_point(size = 3) +
  scale_colour_manual(values = subgroup_pal, name = "Family subgroup") +
  labs(title = "Phoneme PCA (Philippine languages) — coloured by family subgroup",
       x = sprintf("PC1 (%.1f%% var.)", res_ph$var_pct[1]),
       y = sprintf("PC2 (%.1f%% var.)", res_ph$var_pct[2])) +
  theme_minimal()
print(p_pca_subgroup_ph)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_subgroup_ph.png"), p_pca_subgroup_ph,
       width = 7.5, height = 5.5, units = "in", dpi = 300)

# (c) Map: waypoint network, points coloured by DBSCAN cluster
map_points_ph <- PHONEME_cossim %>%
  dplyr::left_join(scores_ph %>% dplyr::select(language, cluster), by = "language")
stopifnot("Some network languages have no PCA/DBSCAN cluster for Part A." =
            !anyNA(map_points_ph$cluster))

network_dbscan_ph <- plot_network_by_cluster(
  map_points_ph, cluster_pal_ph,
  title = "Waypoint network — coloured by DBSCAN cluster (Philippine-only PCA)")
legend_ph <- build_cluster_legend_strip(map_points_ph, cluster_pal_ph)
network_dbscan_ph_full <- (network_dbscan_ph | legend_ph) + plot_layout(widths = c(5, 1.2))
print(network_dbscan_ph_full)
ggsave(here("figures", "phoneme", "mst_waypoints", "PHONEME_network_by_dbscan_ph.png"),
       network_dbscan_ph_full, width = 11, height = 9, units = "in", dpi = 300)

# Loadings bar chart (top |loading| features per PC, IPA-labelled)
p_loadings_ph <- top_loadings_plot(loadings_ph, "PC1_loading", "Top PC1 loadings (Philippine-only)") +
  top_loadings_plot(loadings_ph, "PC2_loading", "Top PC2 loadings (Philippine-only)")
print(p_loadings_ph)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_loadings_ph.png"), p_loadings_ph,
       width = 10, height = 5, units = "in", dpi = 300)

# =============================================================================
# PART B — Exploratory: full database (Philippine + interest + unrelated)
# =============================================================================
res_full <- run_pca_dbscan(weighted_data, minPts, "Phoneme full database",
                           eps_plot_path = here("figures", "phoneme", "pca", "PHONEME_kNNdist_full.png"))
scores_full <- res_full$scores %>%
  dplyr::left_join(RUHLENdf %>% dplyr::select(language, Language_type), by = "language") %>%
  dplyr::relocate(Language_type, .after = language)
loadings_full <- get_loadings(res_full$pca, driver_table)

write_csv(scores_full, here("data", "PHONEME_pca_dbscan_full.csv"))
write_csv(loadings_full, here("data", "PHONEME_pca_loadings_full.csv"))

cluster_pal_full <- build_cluster_pal(levels(scores_full$cluster))
shape_vals <- c("Philippine Language" = 16, "Language of Interest" = 17, "Unrelated Language" = 4)

# (a) PCA scatter, coloured by DBSCAN cluster, shaped by baseline type
p_pca_dbscan_full <- ggplot(scores_full, aes(PC1, PC2, colour = cluster, shape = Language_type)) +
  geom_point(size = 3) +
  scale_colour_manual(values = cluster_pal_full, name = "DBSCAN cluster") +
  scale_shape_manual(values = shape_vals, name = "Baseline type") +
  labs(title = "Phoneme PCA (full database) — coloured by DBSCAN cluster",
       x = sprintf("PC1 (%.1f%% var.)", res_full$var_pct[1]),
       y = sprintf("PC2 (%.1f%% var.)", res_full$var_pct[2])) +
  theme_minimal()
print(p_pca_dbscan_full)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_dbscan_full.png"), p_pca_dbscan_full,
       width = 8, height = 6, units = "in", dpi = 300)

# (b) PCA scatter, coloured by family (PH: subgroup, dropping the (phoneme: none)
#     ungrouped languages; Interest: exact ridge-plot colours; Unrelated: grey70),
#     shaped by baseline type
ridge_colors <- c(Spanish = "#2ca6a4", Japanese = "#8fb339", English = "#e07a5f")

scores_full_family <- scores_full %>%
  dplyr::left_join(subgroup_lookup %>% dplyr::select(language, subgroup), by = "language") %>%
  dplyr::filter(!(Language_type == "Philippine Language" & is.na(subgroup))) %>%
  dplyr::mutate(
    family_key = dplyr::case_when(
      Language_type == "Philippine Language"  ~ subgroup,
      Language_type == "Language of Interest" ~ language,
      Language_type == "Unrelated Language"   ~ "Unrelated"
    )
  )
family_pal_full <- c(subgroup_pal, ridge_colors, "Unrelated" = "grey70")

p_pca_family_full <- ggplot(scores_full_family, aes(PC1, PC2, colour = family_key, shape = Language_type)) +
  geom_point(size = 3) +
  scale_colour_manual(values = family_pal_full, name = "Family / baseline") +
  scale_shape_manual(values = shape_vals, name = "Baseline type") +
  labs(title = "Phoneme PCA (full database) — coloured by language family",
       x = sprintf("PC1 (%.1f%% var.)", res_full$var_pct[1]),
       y = sprintf("PC2 (%.1f%% var.)", res_full$var_pct[2])) +
  theme_minimal()
print(p_pca_family_full)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_family_full.png"), p_pca_family_full,
       width = 9, height = 6.5, units = "in", dpi = 300)

# (c) Map: same waypoint network, points coloured by the FULL-database DBSCAN run
#     (map extent structurally can't show interest/unrelated languages, which sit
#     far outside lon 116-127 / lat 4-21 — only Philippine points are plotted)
map_points_full <- PHONEME_cossim %>%
  dplyr::left_join(scores_full %>% dplyr::select(language, cluster), by = "language")
stopifnot("Some network languages have no PCA/DBSCAN cluster for Part B." =
            !anyNA(map_points_full$cluster))

network_dbscan_full <- plot_network_by_cluster(
  map_points_full, cluster_pal_full,
  title = "Waypoint network — coloured by DBSCAN cluster (full-database PCA)")
legend_full <- build_cluster_legend_strip(map_points_full, cluster_pal_full)
network_dbscan_full_full <- (network_dbscan_full | legend_full) + plot_layout(widths = c(5, 1.2))
print(network_dbscan_full_full)
ggsave(here("figures", "phoneme", "mst_waypoints", "PHONEME_network_by_dbscan_full.png"),
       network_dbscan_full_full, width = 11, height = 9, units = "in", dpi = 300)

# Loadings bar chart (identical construction to Part A, full-database rotation)
p_loadings_full <- top_loadings_plot(loadings_full, "PC1_loading", "Top PC1 loadings (full database)") +
  top_loadings_plot(loadings_full, "PC2_loading", "Top PC2 loadings (full database)")
print(p_loadings_full)
ggsave(here("figures", "phoneme", "pca", "PHONEME_pca_loadings_full.png"), p_loadings_full,
       width = 10, height = 5, units = "in", dpi = 300)
