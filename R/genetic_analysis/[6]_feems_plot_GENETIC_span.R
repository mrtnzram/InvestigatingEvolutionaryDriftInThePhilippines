# =============================================================================
# [6] Genetic Analysis — FEEMS migration-surface base map
# Renders the FEEMS effective-migration surface (the log10(w/w-bar) raster from
# python/genetic_feems.py) as a ggplot base map, with plain points marking
# sampled population locations. No subgroup or admixture colouring here — those
# are separate variables, plotted on their own in [7]'s span_admx figure. Saved
# as an RDS that [7] loads to overlay the waypoint routes.
#
# Input:   data/genetic_surface_raster.csv (FEEMS surface, from genetic_feems.py),
#          data/GENETIC_final.csv (per-population coordinates, from [0])
# Outputs: data/base_plot_genetic_FEEMS.rds
# Next:    [7]_GENETIC_geoplots.R
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)
library(scales)

gen_surface <- read.csv(here("data", "genetic_surface_raster.csv"))
PH_genetic  <- read.csv(here("data", "GENETIC_final.csv"))

vmax <- max(abs(gen_surface$log_w_ratio), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_tile(data = gen_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = PH_genetic, aes(x = longitude, y = latitude),
             size = 1.8, colour = "black", alpha = 0.6) +
  guides(
    fill = guide_colorbar(title = expression(log[10](w/bar(w))), title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    values = rescale(c(-vmax, 0, vmax), from = lims),  # forces white -> exactly 0
    limits = lims,
    na.value = "transparent"
  ) +
  theme_minimal() +
  # Coordinates carry no finding; the coastline is the only reference frame the
  # reader needs, so the axes, ticks and graticule are all dropped.
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

base_plot
saveRDS(base_plot, file = here("data", "base_plot_genetic_FEEMS.rds"))
