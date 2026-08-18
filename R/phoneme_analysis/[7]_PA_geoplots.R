# =============================================================================
# [7] Phoneme Analysis — Geo-plots: overlay waypoint routes on the base maps
# Combines the two [1]/[6] base maps (cosine-similarity points, FEEMS surface
# with subgroup points) with the historical migration-route arrows and marks
# the colonial capital (Manila), then saves each as its own figure. Also
# renders the [6] null-model base map (geographic-shuffle FEEMS null, same
# routes/capital/color scale as the real FEEMS figure) as its own figure.
#
# Input:   data/base_plot_phoneme_cosine.rds (from [1]),
#          data/base_plot_phoneme_FEEMS.rds, data/base_plot_phoneme_FEEMS_null.rds,
#          data/legend_strip_phoneme_subgroup.rds (from [6]),
#          data/phoneme_waypoint_plot.rds (arrow layers, from [3])
# Outputs: figures/phoneme/mst_waypoints/PA_weight_mst_cosine.png,
#          figures/phoneme/mst_waypoints/PA_weight_mst_feems.png,
#          figures/phoneme/mst_waypoints/PA_weight_mst_feems_null.png
# =============================================================================

library(ggplot2)
library(patchwork)
library(here)

# ---- Load the base maps and the waypoint-route arrow layers ----
arrows <- readRDS(here("data", "phoneme", "network_distance", "phoneme_waypoint_plot.rds"))

# ---- Shared overlay: crop to the base map, add the arrows + Manila marker ----
ref_coords1 <- c(121, 14.6)
capital_df  <- data.frame(x = ref_coords1[1], y = ref_coords1[2], label = "Capital")

add_routes_and_capital <- function(base) {
  p <- base + coord_sf(xlim = c(116, 127), ylim = c(4, 21))
  for (layer in arrows$layers) {
    p <- p + layer
  }
  p +
    geom_point(data = capital_df, aes(x = x, y = y, shape = label), color = "red", size = 4) +
    scale_shape_manual(values = c("Capital" = 18)) +
    guides(shape = guide_legend(title = NULL)) +
    theme(legend.position = "right")
}

# ---- Figure A: cosine similarity + routes -----------------------------------
base_plot_cosine <- readRDS(here("data", "phoneme", "cosine_similarity", "base_plot_phoneme_cosine.rds"))
final_plot_cosine <- add_routes_and_capital(base_plot_cosine)

print(final_plot_cosine)

ggsave(here("figures", "phoneme", "mst_waypoints", "PA_weight_mst_cosine.png"),
       final_plot_cosine, width = 7, height = 9, units = "in", dpi = 300)

# ---- Figure B: FEEMS surface + subgroup points + routes ---------------------
base_plot_feems  <- readRDS(here("data", "phoneme", "feems", "base_plot_phoneme_FEEMS.rds"))
legend_strip     <- readRDS(here("data", "phoneme", "feems", "legend_strip_phoneme_subgroup.rds"))

final_plot_feems <- (add_routes_and_capital(base_plot_feems) | legend_strip) +
  plot_layout(widths = c(5, 1.4))

print(final_plot_feems)

ggsave(here("figures", "phoneme", "mst_waypoints", "PA_weight_mst_feems.png"),
       final_plot_feems, width = 9, height = 9, units = "in", dpi = 300)

# ---- Figure C: FEEMS null-model surface + subgroup points + routes ----------
# Geographic-shuffle null (100 permutations, from phoneme_feems.ipynb), on the
# same color scale as Figure B (set jointly in [6]) so the two are comparable.
base_plot_feems_null  <- readRDS(here("data", "phoneme", "feems", "base_plot_phoneme_FEEMS_null.rds"))

final_plot_feems_null <- (add_routes_and_capital(base_plot_feems_null) | legend_strip) +
  plot_layout(widths = c(5, 1.4))

print(final_plot_feems_null)

ggsave(here("figures", "phoneme", "mst_waypoints", "PA_weight_mst_feems_null.png"),
       final_plot_feems_null, width = 9, height = 9, units = "in", dpi = 300)
