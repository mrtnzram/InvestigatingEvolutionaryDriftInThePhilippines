# =============================================================================
# [9] Genetic Analysis — Procrustes: PCA vs. geography
#
# Reuses PLINK's population-structure PCA (data/Phil_2.24M_pca_results.eigenvec,
# the genetic analogue of the other domains' phylogenetic eigenvectors — see
# [4]_GENETIC_PVR.R's header) for the 115-population analysis set. The PCA
# figure plots all 1028 individuals; the Procrustes fit itself runs at
# population grain (mean PC1/PC2 per population) against true longitude/
# latitude, reporting the fit's sum of squares and correlation. Individuals are
# then mapped through that same population-level fit for the Procrustes-fit
# panel, shown alongside the population averages. No fresh PCA is run here —
# there is no local genotype matrix to run it on.
#
# Input:   data/Phil_2.24M_pca_results.eigenvec,
#          data/genetic/network_distance/GENETIC_final.csv,
#          data/genetic/initial_datasets/GENETIC_subgroup_lookup.csv
# Outputs: data/genetic/procrustes/GENETIC_procrustes_results.csv,
#          data/genetic/procrustes/GENETIC_procrustes_scores.csv,
#          data/genetic/procrustes/GENETIC_procrustes_individual_scores.csv,
#          figures/genetic/pca/GENETIC_pca.png,
#          figures/genetic/pca/GENETIC_procrustes.png
# =============================================================================

library(adegenet)   # loads ade4 (procuste, procuste.randtest)
library(tidyverse)
library(here)
library(ggpubr)
library(maps)

dir.create(here("figures", "genetic", "pca"), recursive = TRUE, showWarnings = FALSE)

N_PERM <- 999   # + observed arrangement = 1000 draws

GENETIC_final <- read.csv(here("data", "genetic", "network_distance", "GENETIC_final.csv"))

# ── 1. Population-structure PCA: individual grain + population means ────────
eigenvec_raw <- read.table(here("data", "Phil_2.24M_pca_results.eigenvec"),
                           header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eigenvec_raw) - 2
names(eigenvec_raw) <- c("population", "sample", paste0("PC", seq_len(n_pc)))

# Same population-label reconciliation [4]/[4.1]/[5] apply to the PLINK .mdist.id table.
eigenvec_raw$population <- ifelse(eigenvec_raw$population == "ManoboRajahKabunsuwan",
                                  "ManoboRK", eigenvec_raw$population)

eigenvec_pop <- eigenvec_raw |>
  group_by(population) |>
  summarise(across(starts_with("PC"), mean), .groups = "drop")

scores_df <- eigenvec_pop %>%
  dplyr::select(population, PC1, PC2)

# ── 2. Join coordinates + subgroup colour (population grain, for the fit) ───
subgroup_lookup <- read_csv(here("data", "genetic", "initial_datasets", "GENETIC_subgroup_lookup.csv"), show_col_types = FALSE)

analysis_df <- scores_df %>%
  left_join(GENETIC_final %>% dplyr::select(population, longitude, latitude), by = "population") %>%
  left_join(subgroup_lookup %>% dplyr::select(population, subgroup, colour), by = "population") %>%
  arrange(population)

stopifnot(
  "Some PCA rows are missing a coordinate." = !anyNA(analysis_df$longitude),
  "Some PCA rows are missing a subgroup colour." = !anyNA(analysis_df$colour)
)
N <- nrow(analysis_df)

# Individual-grain rows carry the same subgroup/colour as their population
# (no per-individual coordinate exists — genetic samples are geolocated at
# their population's site), used for the PCA figure and, after the fit below,
# for showing individual spread on the Procrustes-fit panel.
indiv_df <- eigenvec_raw %>%
  left_join(subgroup_lookup %>% dplyr::select(population, subgroup, colour), by = "population") %>%
  arrange(population, sample)

stopifnot("Some individual rows are missing a subgroup colour." = !anyNA(indiv_df$colour))
N_indiv <- nrow(indiv_df)

# ── 3. Procrustes fit (population grain) ─────────────────────────────────────
pc_mat  <- as.data.frame(analysis_df[, c("PC1", "PC2")])
geo_mat <- as.data.frame(analysis_df[, c("longitude", "latitude")])
rownames(pc_mat) <- rownames(geo_mat) <- analysis_df$population

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
# population means): centers/scales by the fitted pc_mat's own mean/norm,
# applies the fitted rotation (loadX/loadY, i.e. u/v from ade4's internal
# SVD), then de-normalizes by geo_mat's mean/norm. Used both to draw the
# PC1/PC2 axes in geographic space and to map each individual through the
# population-level fit below (verified to reproduce pr$rotX exactly for the
# real population-mean data).
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

# Individuals mapped through the population-level fit — not re-fit, just
# projected via the same pc_center/normX/u_mat/v_mat/normY/geo_center already
# derived from the 115-population fit above.
indiv_rot_geo <- pc_to_geo(indiv_df$PC1, indiv_df$PC2)
indiv_df$rotLongitude <- indiv_rot_geo[, 1]
indiv_df$rotLatitude  <- indiv_rot_geo[, 2]

# ── 4. Results + scores tables ───────────────────────────────────────────────
results_df <- tibble(domain = "genetic", n = N, n_individual = N_indiv, ss_gower = ss,
                      correlation = correlation, p_value = p_value, n_perm = N_PERM)
print(results_df)
write_csv(results_df, here("data", "genetic", "procrustes", "GENETIC_procrustes_results.csv"))
write_csv(analysis_df, here("data", "genetic", "procrustes", "GENETIC_procrustes_scores.csv"))
write_csv(indiv_df, here("data", "genetic", "procrustes", "GENETIC_procrustes_individual_scores.csv"))

# ── 5. PCA figure (own file, individual grain) ───────────────────────────────
subgroup_pal <- analysis_df %>% distinct(subgroup, colour) %>% deframe()

p_pca <- ggplot(indiv_df, aes(PC1, PC2, colour = subgroup)) +
  geom_point(size = 3,alpha = 0.8) +
  scale_colour_manual(values = subgroup_pal, name = "Subgroup") +
  guides(colour = guide_legend(ncol = 1)) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)) +
  theme_bw() +
  labs(title = "PCA")
print(p_pca)
ggsave(here("figures", "genetic", "pca", "GENETIC_pca.png"), p_pca,
       width = 8, height = 7, units = "in", dpi = 300, bg = "white")

# ── 6. Procrustes-fit vs. actual-coordinates figure (own file) ──────────────
# Map styling matches the project's waypoint-map convention ([1]'s
# base_plot_phoneme_cosine.rds / the [7] geoplots): map_data("world") polygon,
# gray95/gray70 land, coord_fixed(115:130, 4:22), blank axes/grid. No MST
# edges or route arrows — those belong to the [3]/[7] waypoint figures, not
# this comparison. The Procrustes-fit panel layers individuals (crosses, low
# alpha, mapped through the population fit) under population means (circles).
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
  geom_point(data = indiv_df, aes(x = rotLongitude, y = rotLatitude, colour = subgroup),
             shape = 4, size = 2, alpha = 0.35) +
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
p_maps   <- ggarrange(maps_row, as_ggplot(shared_legend), ncol = 2, nrow = 1, widths = c(4, 1.3))
print(p_maps)
ggsave(here("figures", "genetic", "pca", "GENETIC_procrustes.png"), p_maps,
       width = 12, height = 7, units = "in", dpi = 300, bg = "white")
