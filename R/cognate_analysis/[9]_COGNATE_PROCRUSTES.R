# =============================================================================
# [9] Cognate Analysis — Procrustes: PCA vs. geography
#
# Runs PCA (adegenet::dudi.pca) on a one-hot (gloss, wordform) presence/absence
# matrix built from the raw ABVD wordlists, restricted to the 94-language
# analysis set, then fits a Procrustes rotation (ade4::procuste, via adegenet)
# of PC1/PC2 onto true longitude/latitude, reporting the fit's sum of squares
# and correlation.
#
# The one-hot matrix construction (split/trim/casefold each cell, one binary
# feature per distinct (gloss, wordform) pair, drop pairs attested in fewer
# than MIN_LANGS languages) matches [4.2]_COGNATE_sPCA.R's §2, minus that
# script's phylogenetic residualization step — this PCA is unconstrained.
#
# Input:   data/cognate/PH_df.csv (raw per-gloss wordforms),
#          data/cognate/COGNATE_final.csv (analysis-set glottocodes),
#          data/cognate/COGNATE_subgroup_lookup.csv
# Outputs: data/cognate/COGNATE_procrustes_results.csv,
#          data/cognate/COGNATE_procrustes_scores.csv,
#          figures/cognate/pca/COGNATE_procrustes.png
# =============================================================================

library(adegenet)   # loads ade4 (dudi.pca, procuste, procuste.randtest)
library(tidyverse)
library(here)
library(ggpubr)
library(maps)

dir.create(here("figures", "cognate", "pca"), recursive = TRUE, showWarnings = FALSE)

N_PERM    <- 999   # + observed arrangement = 1000 draws
MIN_LANGS <- 2      # drop (gloss, wordform) pairs attested in fewer languages

COGNATE_final <- read.csv(here("data", "cognate", "network_distance", "COGNATE_final.csv"))
analysis_glottocodes <- COGNATE_final$glottocode

# ── 1. One-hot (gloss, wordform) matrix at language grain ───────────────────
META_COLS <- c("language_id", "language", "latitude", "longitude",
               "glottocode", "author", "number_of_entries", "source_url")
PH_df <- read.csv(here("data", "cognate", "initial_datasets", "PH_df.csv"))
gloss_cols <- setdiff(names(PH_df), META_COLS)

ph_sub <- PH_df |>
  dplyr::select(glottocode, language, longitude, latitude, all_of(gloss_cols)) |>
  filter(glottocode %in% analysis_glottocodes)

long <- ph_sub |>
  dplyr::select(glottocode, all_of(gloss_cols)) |>
  pivot_longer(-glottocode, names_to = "gloss", values_to = "cell") |>
  filter(!is.na(cell), cell != "") |>
  separate_longer_delim(cell, "; ") |>
  mutate(form = str_trim(tolower(cell))) |>
  filter(form != "") |>
  distinct(glottocode, gloss, form)

pair_counts <- long |> count(gloss, form, name = "n_langs")
kept_pairs  <- pair_counts |> filter(n_langs >= MIN_LANGS)
long <- long |> semi_join(kept_pairs, by = c("gloss", "form"))

feature_wide <- long |>
  mutate(feature = paste(gloss, form, sep = "::"), present = 1L) |>
  distinct(glottocode, feature, present) |>
  pivot_wider(names_from = feature, values_from = present, values_fill = 0L)
feature_cols <- setdiff(names(feature_wide), "glottocode")

df <- ph_sub |>
  dplyr::select(glottocode, language, longitude, latitude) |>
  left_join(feature_wide, by = "glottocode") |>
  mutate(across(all_of(feature_cols), \(x) replace_na(x, 0L))) |>
  arrange(glottocode) |>
  as.data.frame()

X <- as.matrix(df[, feature_cols])
X <- X[, apply(X, 2, var) > 0, drop = FALSE]
# Distinct-presence-pattern duplicates add no new direction to the column
# space (see [4.2]'s note) — one representative per pattern is kept.
X <- X[, !duplicated(t(X)), drop = FALSE]
rownames(X) <- df$glottocode
message(nrow(df), " languages, ", ncol(X), " (gloss, wordform) features retained.")

# ── 2. PCA (adegenet/ade4) ───────────────────────────────────────────────────
pca_fit <- dudi.pca(as.data.frame(X), center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
scores_df <- tibble(glottocode = rownames(X), PC1 = pca_fit$li[, 1], PC2 = pca_fit$li[, 2])

# ── 3. Join coordinates + subgroup colour ────────────────────────────────────
subgroup_lookup <- read_csv(here("data", "cognate", "initial_datasets", "COGNATE_subgroup_lookup.csv"), show_col_types = FALSE)

analysis_df <- scores_df %>%
  left_join(df %>% dplyr::select(glottocode, language, longitude, latitude), by = "glottocode") %>%
  left_join(subgroup_lookup %>% dplyr::select(glottocode, subgroup, colour), by = "glottocode") %>%
  arrange(glottocode)

stopifnot(
  "Some PCA rows are missing a coordinate." = !anyNA(analysis_df$longitude),
  "Some PCA rows are missing a subgroup colour." = !anyNA(analysis_df$colour)
)
N <- nrow(analysis_df)

# ── 4. Procrustes fit ────────────────────────────────────────────────────────
pc_mat  <- as.data.frame(analysis_df[, c("PC1", "PC2")])
geo_mat <- as.data.frame(analysis_df[, c("longitude", "latitude")])
rownames(pc_mat) <- rownames(geo_mat) <- analysis_df$glottocode

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
results_df <- tibble(domain = "cognate", n = N, ss_gower = ss,
                      correlation = correlation, p_value = p_value, n_perm = N_PERM)
print(results_df)
write_csv(results_df, here("data", "cognate", "procrustes", "COGNATE_procrustes_results.csv"))
write_csv(analysis_df, here("data", "cognate", "procrustes", "COGNATE_procrustes_scores.csv"))

# ── 6. PCA figure (own file) ─────────────────────────────────────────────────
subgroup_pal <- analysis_df %>% distinct(subgroup, colour) %>% deframe()

p_pca <- ggplot(analysis_df, aes(PC1, PC2, colour = subgroup)) +
  geom_point(size = 3) +
  scale_colour_manual(values = subgroup_pal, name = "Subgroup") +
  guides(colour = guide_legend(ncol = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  labs(title = "PCA")
print(p_pca)
ggsave(here("figures", "cognate", "pca", "COGNATE_pca.png"), p_pca,
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
ggsave(here("figures", "cognate", "pca", "COGNATE_procrustes.png"), p_maps,
       width = 11, height = 6, units = "in", dpi = 300, bg = "white")
