# =============================================================================
# [7] Genetic Analysis — Geo-plot: Spanish admixture over the waypoint network
# Maps each population's qpAdm Spanish ancestry proportion on the Philippine base
# map, with the [3] migration-route arrows underneath and point area encoding
# span_w.
#
# Input:   data/GENETIC_final.csv, data/genetic_waypoint_plot.rds (from [3])
# Outputs: figures/genetic/mst_waypoints/GENETIC_span_admx_weighted.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(maps)
library(sf)
library(here)

PH_genetic <- read.csv(here("data", "GENETIC_final.csv"))
arrows     <- readRDS(here("data", "genetic_waypoint_plot.rds"))

points_df <- PH_genetic |> filter(!is.na(span_admx), !is.na(span_w))

# SE caps recovered from the capped column so the legend endpoints track [3].
q_lo <- min(points_df$span_se_capped)
q_hi <- max(points_df$span_se_capped)

# Legend breaks are set in SE units and inverted to weights, so size reads as
# precision.
se_breaks  <- c(q_hi, 0.15, 0.0175, q_lo)
size_labels <- c(sprintf(">= %.2f", q_hi), "0.15", "0.018",
                 sprintf("<= %.3f", q_lo))

# ---- Base map ---------------------------------------------------------------
map_subset <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  theme_minimal() +
  # Coordinates carry no finding; the coastline is the only reference frame the
  # reader needs, so the axes, ticks and graticule are all dropped.
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

# ---- Route overlay: crop to the base map, replay the [3] arrow layers -------
ref_coords1 <- c(121, 14.6)
capital_df  <- data.frame(x = ref_coords1[1], y = ref_coords1[2], label = "Capital")

add_routes <- function(base) {
  p <- base + coord_sf(xlim = c(116, 127), ylim = c(4, 21))
  for (layer in arrows$layers) {
    p <- p + layer
  }
  p
}

# ---- Spanish admixture points ----------------------------------------------
# White -> navy ramp anchored at 0, drawn as a translucent disc under a black
# ring, matching the [1] cosine base maps; both layers share the `size` aesthetic
# so the ring cannot drift off the disc. Capital marker goes last to stay on top.
global_lim <- c(0, max(points_df$span_admx, na.rm = TRUE))

final_plot <- add_routes(base_plot) +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude, colour = span_admx, size = span_w),
             alpha = 0.7) +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude, size = span_w),
             shape = 21, colour = "black") +
  geom_point(data = capital_df, aes(x = x, y = y, shape = label),
             color = "red", size = 4) +
  scale_colour_gradient(low = "white", high = "navy", limits = global_lim) +
  scale_size(
    trans  = "log10",
    range  = c(1, 5.5),
    breaks = 1 / se_breaks^2,
    labels = size_labels
  ) +
  scale_shape_manual(values = c("Capital" = 18)) +
  guides(
    colour = guide_colorbar(title = "Spanish ancestry",
                            title.position = "top", title.hjust = 0.5, order = 1),
    size   = guide_legend(title = "qpAdm SE", override.aes = list(colour = "grey40"), order = 2),
    shape  = guide_legend(title = NULL, order = 3)
  ) +
  theme(legend.position = "right")

print(final_plot)

dir.create(here("figures", "genetic", "mst_waypoints"),
           recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "genetic", "mst_waypoints", "GENETIC_span_admx_weighted.png"),
       final_plot, width = 7, height = 9, units = "in", dpi = 300)
