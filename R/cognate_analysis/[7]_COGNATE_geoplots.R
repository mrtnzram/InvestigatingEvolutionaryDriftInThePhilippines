# =============================================================================
# [7] Cognate Analysis — Geo-plot: Spanish loanwords over the waypoint network
# Maps each language's normalized Spanish loanword count on the Philippine base
# map, with the [3] migration-route arrows underneath.
#
# Input:   data/cognate/COGNATE_final.csv, data/cognate/cognate_waypoint_plot.rds (from [3])
# Outputs: figures/cognate/mst_waypoints/COGNATE_loans.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(maps)
library(sf)
library(here)

PH_cognate <- read.csv(here("data", "cognate", "COGNATE_final.csv"))
arrows     <- readRDS(here("data", "cognate", "cognate_waypoint_plot.rds"))

points_df <- PH_cognate |> filter(!is.na(loans_norm))

# Point size is constant: with no standard error there is nothing for point area
# to encode, and 108 languages crowd the map.
POINT_SIZE <- 2

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

# ---- Loanword points --------------------------------------------------------
# White -> navy ramp over the full [0, 1] normalized range, drawn as a
# translucent disc under a black ring; both layers use POINT_SIZE so the ring
# cannot drift off the disc. Capital marker goes last to stay on top.
final_plot <- add_routes(base_plot) +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude, colour = loans_norm),
             size = POINT_SIZE, alpha = 0.7) +
  geom_point(data = points_df,
             aes(x = longitude, y = latitude),
             shape = 21, colour = "black", size = POINT_SIZE) +
  geom_point(data = capital_df, aes(x = x, y = y, shape = label),
             color = "red", size = 4) +
  scale_colour_gradient(low = "white", high = "navy", limits = c(0, 1)) +
  scale_shape_manual(values = c("Capital" = 18)) +
  guides(
    colour = guide_colorbar(title = "Spanish loanwords",
                            title.position = "top", title.hjust = 0.5, order = 1),
    shape  = guide_legend(title = NULL, order = 2)
  ) +
  theme(legend.position = "right")

print(final_plot)

dir.create(here("figures", "cognate", "mst_waypoints"),
           recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "cognate", "mst_waypoints", "COGNATE_loans.png"),
       final_plot, width = 7, height = 9, units = "in", dpi = 300)
