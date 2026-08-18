# =============================================================================
# [9] Grammar Analysis — Procrustes: PCA vs. geography
#
# Runs PCA (adegenet::dudi.pca) on the IDF-weighted GRAMBANK feature matrix (the
# matrix behind [1]'s cosine similarity, not the cosine matrix) for the 26
# Philippine languages, then fits a Procrustes rotation (ade4::procuste, via
# adegenet) of PC1/PC2 onto true longitude/latitude, reporting the fit's sum of
# squares and correlation.
#
# Run order: needs [0]_Phylogenetic_Tree.R to have been run
# (data/GRAMMAR_subgroup_lookup.csv). NOTE: [1]'s OHE + per-value IDF weighting
# is duplicated below because [1] never exposes its weighted matrix — a change
# to that scheme has to be mirrored here.
#
# Input:   data/GRAMBANKdf_full.csv, data/gramfeature_freq.csv,
#          data/GRAMMAR_subgroup_lookup.csv
# Outputs: data/GRAMMAR_procrustes_results.csv, data/GRAMMAR_procrustes_scores.csv,
#          figures/grammar/pca/GRAMMAR_pca.png,
#          figures/grammar/pca/GRAMMAR_procrustes.png
# =============================================================================

library(adegenet)   # loads ade4 (dudi.pca, procuste, procuste.randtest)
library(tidyverse)
library(here)
library(ggpubr)
library(maps)
library(fastDummies)

dir.create(here("figures", "grammar", "pca"), recursive = TRUE, showWarnings = FALSE)

N_PERM <- 999   # + observed arrangement = 1000 draws

# ── 1. Rebuild the IDF-weighted GRAMBANK matrix (duplicated from [1]) ───────
GRAMBANKdf_PH <- read_csv(here("data", "grammar", "initial_datasets", "GRAMBANKdf_full.csv"), show_col_types = FALSE)

feature_cols <- GRAMBANKdf_PH %>%
  dplyr::select(-Language_ID, -language, -Family_name, -Macroarea, -Longitude, -Latitude, -Language_Type) %>%
  colnames()

gramfeature_freq <- read_csv(here("data", "grammar", "initial_datasets", "gramfeature_freq.csv"), show_col_types = FALSE)

n_levels <- purrr::map_int(
  feature_cols,
  \(col) GRAMBANKdf_PH[[col]] |> na.omit() |> unique() |> length()
)
binary_cols      <- feature_cols[n_levels <= 2]
categorical_cols <- feature_cols[n_levels >  2]

GRAMBANKdf_PH_ohe <- GRAMBANKdf_PH |>
  dummy_cols(
    select_columns          = categorical_cols,
    remove_selected_columns = TRUE,
    remove_first_dummy      = FALSE,
    ignore_na                = TRUE
  )
ohe_pattern <- str_c("^(", str_c(categorical_cols, collapse = "|"), ")_")
ohe_cols    <- str_subset(names(GRAMBANKdf_PH_ohe), ohe_pattern)
feature_cols_ohe <- c(binary_cols, ohe_cols)
GRAMBANKdf_PH <- GRAMBANKdf_PH_ohe

# IDF-substitution weighting, matching [1]: OHE dummies take
# IDF(base_feature, level), binary columns take IDF(observed value).
long_data <- GRAMBANKdf_PH %>%
  dplyr::select(dplyr::all_of(feature_cols_ohe)) %>%
  dplyr::mutate(language = GRAMBANKdf_PH[["language"]]) %>%
  tidyr::pivot_longer(cols = -language, names_to = "feature", values_to = "value")

is_ohe <- long_data$feature %in% ohe_cols
long_data <- long_data %>%
  dplyr::mutate(
    lookup_feature = dplyr::if_else(is_ohe, stringr::str_remove(feature, "_[^_]+$"), feature),
    lookup_value   = dplyr::if_else(is_ohe, as.numeric(stringr::str_extract(feature, "(?<=_)[^_]+$")), value)
  )

weighted_long <- long_data %>%
  dplyr::left_join(gramfeature_freq, by = c("lookup_feature" = "feature", "lookup_value" = "value")) %>%
  dplyr::mutate(weighted_value = dplyr::if_else(feature %in% ohe_cols, value * IDF, IDF))

weighted_data <- weighted_long %>%
  dplyr::select(language, feature, weighted_value) %>%
  tidyr::pivot_wider(names_from = feature, values_from = weighted_value) %>%
  tibble::column_to_rownames("language") %>%
  as.matrix()
weighted_data <- weighted_data[GRAMBANKdf_PH$language, feature_cols_ohe, drop = FALSE]

# Spanish's GB415 has no WALS mapping and is the one expected NA; it is imputed
# to 0 (the OHE columns' native "absent" value) rather than dropping the row.
stopifnot("More NA cells than expected in the weighted matrix — investigate before imputing." =
            sum(is.na(weighted_data)) <= 1)
if (sum(is.na(weighted_data)) == 1) {
  na_idx <- which(is.na(weighted_data), arr.ind = TRUE)
  message("Imputing 1 NA cell to 0: language=", rownames(weighted_data)[na_idx[1, "row"]],
          ", feature=", colnames(weighted_data)[na_idx[1, "col"]])
  weighted_data[na_idx] <- 0
}

ph_mat <- weighted_data[GRAMBANKdf_PH$Language_Type == "Philippine Language", , drop = FALSE]
ph_mat <- ph_mat[, apply(ph_mat, 2, var) > 0, drop = FALSE]

ph_coords <- GRAMBANKdf_PH %>%
  dplyr::filter(Language_Type == "Philippine Language") %>%
  dplyr::select(language, longitude = Longitude, latitude = Latitude)

# ── 2. PCA (adegenet/ade4) ───────────────────────────────────────────────────
pca_fit <- dudi.pca(as.data.frame(ph_mat), center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
scores_df <- tibble(language = rownames(ph_mat), PC1 = pca_fit$li[, 1], PC2 = pca_fit$li[, 2])

# ── 3. Join coordinates + subgroup colour ────────────────────────────────────
subgroup_lookup <- read_csv(here("data", "grammar", "initial_datasets", "GRAMMAR_subgroup_lookup.csv"), show_col_types = FALSE)

analysis_df <- scores_df %>%
  left_join(ph_coords, by = "language") %>%
  left_join(subgroup_lookup %>% dplyr::select(language, subgroup, colour), by = "language") %>%
  arrange(language)

stopifnot(
  "Some PCA rows are missing a coordinate." = !anyNA(analysis_df$longitude),
  "Some PCA rows are missing a subgroup colour." = !anyNA(analysis_df$colour)
)
N <- nrow(analysis_df)

# ── 4. Procrustes fit ────────────────────────────────────────────────────────
pc_mat  <- as.data.frame(analysis_df[, c("PC1", "PC2")])
geo_mat <- as.data.frame(analysis_df[, c("longitude", "latitude")])
rownames(pc_mat) <- rownames(geo_mat) <- analysis_df$language

pr  <- procuste(pc_mat, geo_mat, scale = TRUE, nf = 2)
ss  <- sum((as.matrix(pr$rotX) - as.matrix(pr$tabY))^2)   # Gower's M^2, dimensionless
prt <- procuste.randtest(pc_mat, geo_mat, nrepet = N_PERM)
correlation <- prt$obs                                    # == sum(pr$d)
p_value     <- prt$pvalue

analysis_df$rotPC1 <- pr$rotX[[1]]
analysis_df$rotPC2 <- pr$rotX[[2]]

# De-normalize the Procrustes fit back to real degree units: rotX and tabY
# share the same ade4 unit-sum-of-squares frame, so inverting the same
# centering/scaling ade4 applied to geo_mat recovers a real-map position for
# each rotated PCA point (verified to round-trip exactly on tabY itself).
normY      <- sqrt(sum(scale(as.matrix(geo_mat), scale = FALSE)^2))
geo_center <- colMeans(as.matrix(geo_mat))
rot_geo    <- sweep(as.matrix(pr$rotX) * normY, 2, geo_center, FUN = "+")
analysis_df$rotLongitude <- rot_geo[, 1]
analysis_df$rotLatitude  <- rot_geo[, 2]

# Same transform, applied to arbitrary PC1/PC2 points (not just the observed
# rows): centers/scales by the fitted pc_mat's own mean/norm, applies the
# fitted rotation (loadX/loadY, i.e. u/v from ade4's internal SVD), then
# de-normalizes by geo_mat's mean/norm. Used below to draw the PC1/PC2 axes
# in geographic space (verified to reproduce pr$rotX exactly for the real data).
pc_center <- colMeans(as.matrix(pc_mat))
normX     <- sqrt(sum(scale(as.matrix(pc_mat), scale = FALSE)^2))
u_mat     <- as.matrix(pr$loadX)[c("PC1", "PC2"), ]
v_mat     <- as.matrix(pr$loadY)[c("longitude", "latitude"), ]

pc_to_geo <- function(pc1, pc2) {
  p <- sweep(cbind(PC1 = pc1, PC2 = pc2), 2, pc_center, "-") / normX
  sweep((p %*% u_mat %*% t(v_mat)) * normY, 2, geo_center, "+")
}

axis_lines <- bind_rows(
  as_tibble(pc_to_geo(range(analysis_df$PC1), rep(mean(analysis_df$PC2), 2))) %>%
    mutate(axis = "PC1", end = c("start", "end")),
  as_tibble(pc_to_geo(rep(mean(analysis_df$PC1), 2), range(analysis_df$PC2))) %>%
    mutate(axis = "PC2", end = c("start", "end"))
) %>%
  rename(x = longitude, y = latitude) %>%
  pivot_wider(names_from = end, values_from = c(x, y))

# ── 5. Results + scores tables ───────────────────────────────────────────────
results_df <- tibble(domain = "grammar", n = N, ss_gower = ss,
                      correlation = correlation, p_value = p_value, n_perm = N_PERM)
print(results_df)
write_csv(results_df, here("data", "grammar", "procrustes", "GRAMMAR_procrustes_results.csv"))
write_csv(analysis_df, here("data", "grammar", "procrustes", "GRAMMAR_procrustes_scores.csv"))

# ── 6. PCA figure (own file) ─────────────────────────────────────────────────
subgroup_pal <- analysis_df %>% distinct(subgroup, colour) %>% deframe()

p_pca <- ggplot(analysis_df, aes(PC1, PC2, colour = subgroup)) +
  geom_point(size = 3) +
  scale_colour_manual(values = subgroup_pal, name = "Subgroup") +
  guides(colour = guide_legend(ncol = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  labs(title = "PCA")
print(p_pca)
ggsave(here("figures", "grammar", "pca", "GRAMMAR_pca.png"), p_pca,
       width = 8, height = 6, units = "in", dpi = 300, bg = "white")

# ── 7. Procrustes-fit vs. actual-coordinates figure (own file) ──────────────
# Map styling matches the project's waypoint-map convention ([1]'s
# base_plot_phoneme_cosine.rds / the [7] geoplots): map_data("world") polygon,
# gray95/gray70 land, coord_fixed(115:130, 4:22), blank axes/grid. No MST
# edges or route arrows — those belong to the [3]/[7] waypoint figures, not
# this comparison.
map_theme <- theme_minimal() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7))

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

p_proc <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_segment(data = axis_lines, aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
               linetype = "dashed", colour = "grey30", linewidth = 0.6, inherit.aes = FALSE) +
  geom_point(data = analysis_df, aes(x = rotLongitude, y = rotLatitude, colour = subgroup),
             size = 3) +
  scale_colour_manual(values = subgroup_pal, name = "Subgroup") +
  guides(colour = guide_legend(ncol = 1)) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  map_theme

p_coords <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = analysis_df, aes(x = longitude, y = latitude, colour = subgroup),
             size = 3) +
  scale_colour_manual(values = subgroup_pal, name = "Subgroup") +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  map_theme

shared_legend <- get_legend(p_proc)
p_proc   <- p_proc   + theme(legend.position = "none")
p_coords <- p_coords + theme(legend.position = "none")

maps_row <- ggarrange(p_proc, p_coords, ncol = 2, nrow = 1, labels = c("A", "B"))
p_maps   <- ggarrange(maps_row, as_ggplot(shared_legend), ncol = 2, nrow = 1, widths = c(4, 1))
print(p_maps)
ggsave(here("figures", "grammar", "pca", "GRAMMAR_procrustes.png"), p_maps,
       width = 11, height = 6, units = "in", dpi = 300, bg = "white")
