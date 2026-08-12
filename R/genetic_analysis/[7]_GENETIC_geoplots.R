# =============================================================================
# [7] Genetic Analysis — Geo-plots: overlay waypoint routes on the base maps
# Combines two independent base maps — Spanish admixture points on a plain
# outline, and the [6] FEEMS migration surface — with the [3] migration-route
# arrows and marks the colonial capital (Manila), then saves each as its own
# figure.
#
# Input:   data/GENETIC_final.csv (per-population admixture + coordinates, from [0]),
#          data/base_plot_genetic_FEEMS.rds (from [6]),
#          data/genetic_waypoint_plot.rds (arrow layers, from [3])
# Outputs: figures/genetic/mst_waypoints/GENETIC_span_admx_weighted.png,
#          figures/genetic/mst_waypoints/GENETIC_weight_mst_feems.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)

PH_genetic <- read.csv(here("data", "GENETIC_final.csv"))
arrows     <- readRDS(here("data", "genetic_waypoint_plot.rds"))

points_df <- PH_genetic |> filter(!is.na(span_admx), !is.na(span_w))

# SE caps recovered from the capped column so the legend endpoints track [3].
q_lo <- min(points_df$span_se_capped)
q_hi <- max(points_df$span_se_capped)

# Legend breaks are set in SE units and inverted to weights, so size reads as
# precision.
se_breaks   <- c(q_hi, 0.15, 0.0175, q_lo)
size_labels <- c(sprintf(">= %.2f", q_hi), "0.15", "0.018",
                 sprintf("<= %.3f", q_lo))

# ---- Shared overlay: crop to the base map, add the [3] arrows + Manila marker
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

dir.create(here("figures", "genetic", "mst_waypoints"),
           recursive = TRUE, showWarnings = FALSE)

# ---- Figure A: Spanish admixture + routes -----------------------------------
# White -> navy ramp anchored at 0, drawn as a translucent disc under a black
# ring; both layers share the `size` aesthetic so the ring cannot drift off
# the disc.
global_lim <- c(0, max(points_df$span_admx, na.rm = TRUE))

map_subset <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

base_plot_admx <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude, colour = span_admx, size = span_w),
             alpha = 0.7) +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude, size = span_w),
             shape = 21, colour = "black") +
  scale_colour_gradient(low = "white", high = "navy", limits = global_lim) +
  scale_size(
    trans  = "log10",
    range  = c(1, 5.5),
    breaks = 1 / se_breaks^2,
    labels = size_labels
  ) +
  guides(
    colour = guide_colorbar(title = "Spanish ancestry",
                            title.position = "top", title.hjust = 0.5, order = 1),
    size   = guide_legend(title = "qpAdm SE", override.aes = list(colour = "grey40"), order = 2)
  ) +
  theme_minimal() +
  # Coordinates carry no finding; the coastline is the only reference frame the
  # reader needs, so the axes, ticks and graticule are all dropped.
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

final_plot_admx <- add_routes_and_capital(base_plot_admx)

print(final_plot_admx)

ggsave(here("figures", "genetic", "mst_waypoints", "GENETIC_span_admx_weighted.png"),
       final_plot_admx, width = 7, height = 9, units = "in", dpi = 300)

# ---- Figure B: FEEMS surface + routes ----------------------------------------
base_plot_feems  <- readRDS(here("data", "base_plot_genetic_FEEMS.rds"))
final_plot_feems <- add_routes_and_capital(base_plot_feems)

print(final_plot_feems)

ggsave(here("figures", "genetic", "mst_waypoints", "GENETIC_weight_mst_feems.png"),
       final_plot_feems, width = 7, height = 9, units = "in", dpi = 300)
