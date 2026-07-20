# =============================================================================
# [6] Grammar Analysis — FEEMS migration-surface base map
# Renders the FEEMS effective-migration surface (the log10(w/w-bar) raster from
# python/grammar_feems.ipynb) as a ggplot base map of the Philippines, with each
# language's cosine-similarity score overlaid as points. Saved as an RDS that
# [7] loads to overlay the waypoint routes.
#
# Input:   data/grammar_surface_raster.csv (FEEMS surface, from grammar_feems.ipynb),
#          data/GRAMMAR_cossim.csv (per-language cosine scores, from [1])
# Outputs: data/base_plot_grammar_FEEMS.rds
# Next:    [7]_GA_weight_mst_feems_span.R
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)
library(scales)

gram_surface   <- read.csv(here("data", "grammar_surface_raster.csv"))
GRAMMAR_cossim <- read.csv(here("data", "GRAMMAR_cossim.csv"))

# Overlay only the 26 tree-pruned Philippine languages, matching [7]'s waypoint
# arrows (from [3]) and [8]'s subgroup network — the 13 languages without an ABVD
# tree tip are not shown as points. The FEEMS surface itself is still estimated
# from all 39 languages (grammar_feems.ipynb) — dropping the southern Bornean
# samples would leave that part of the map purely interpolated — so only the
# discrete cosine-similarity markers are restricted here, not the raster.
tree_langs     <- read.csv(here("data", "GRAMMAR_subgroup_lookup.csv"))$language
GRAMMAR_cossim <- GRAMMAR_cossim[GRAMMAR_cossim$language %in% tree_langs, ]

global_lim <- c(0, max(GRAMMAR_cossim$cossim_span, na.rm = TRUE))

vmax <- max(abs(gram_surface$log_w_ratio), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_tile(data = gram_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = GRAMMAR_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 6, alpha = 0.7) +
  geom_point(data = GRAMMAR_cossim, aes(x = longitude, y = latitude),
             size = 6, shape = 21, color = "black") +
  scale_color_gradient(low = "white", high = "navy", limits = global_lim) +
  guides(
    fill = guide_colorbar(title = expression(log[10](w/bar(w))), title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title = "Cosine Similarity", title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_x_continuous(breaks = seq(115, 130, by = 2)) +
  scale_y_continuous(breaks = seq(4, 22, by = 2)) +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    values = rescale(c(-vmax, 0, vmax), from = lims),  # forces white -> exactly 0
    limits = lims,
    na.value = "transparent"
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

base_plot
saveRDS(base_plot, file = here("data", "base_plot_grammar_FEEMS.rds"))
