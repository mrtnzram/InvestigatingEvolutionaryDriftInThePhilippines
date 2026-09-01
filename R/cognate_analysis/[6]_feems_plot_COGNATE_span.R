# =============================================================================
# [6] Cognate Analysis — FEEMS migration-surface base map
# Renders the FEEMS effective-migration surface (the log10(w/w-bar) raster from
# python/cognate_feems.ipynb) as a ggplot base map, with each language's
# Glottolog subgroup membership overlaid as points (same palette/legend
# construction as the phoneme/grammar/genetic FEEMS maps), plus a matching
# legend strip. Also renders the geographic-shuffle null model's averaged
# surface (cognate_feems.ipynb's permutation cells) as a second base map, on
# the SAME color scale as the real one, so the two are directly comparable.
#
# Unlike phoneme/genetic, cognate has no [7]-equivalent MST/waypoint compositor
# to hand the RDS objects to ([7]_COGNATE_geoplots.R is an unrelated Spanish-
# loanword map) — both RDS objects are still saved for consistency with the
# other domains, but this script also ggsaves the two finished maps directly.
#
# Input:   data/cognate/feems/cognate_surface_raster.csv,
#          data/cognate/feems/cognate_surface_raster_null_mean.csv
#          (FEEMS surface + null-model average, from python/cognate_feems.ipynb),
#          data/cognate/network_distance/COGNATE_final.csv (per-language coordinates),
#          data/cognate/initial_datasets/COGNATE_subgroup_lookup.csv (subgroup + colour)
# Outputs: data/cognate/feems/base_plot_cognate_FEEMS.rds,
#          data/cognate/feems/base_plot_cognate_FEEMS_null.rds,
#          data/cognate/feems/legend_strip_cognate_subgroup.rds,
#          figures/cognate/feems/cognate_feems_surface.png,
#          figures/cognate/feems/cognate_feems_surface_null.png
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)
library(scales)
library(maps)

dir.create(here("figures", "cognate", "feems"), recursive = TRUE, showWarnings = FALSE)

cog_surface      <- read.csv(here("data", "cognate", "feems", "cognate_surface_raster.csv"))
cog_surface_null <- read.csv(here("data", "cognate", "feems", "cognate_surface_raster_null_mean.csv"))
COGNATE_final    <- read.csv(here("data", "cognate", "network_distance", "COGNATE_final.csv"))

# Shared color scale across the real and null maps (not each normalized to its own
# max) — this is what makes a flatter null surface visibly wash out next to the real
# one.
vmax <- max(abs(cog_surface$log_w_ratio), abs(cog_surface_null$log_w_ratio_null_mean), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

# ---- Subgroup palette + membership -------------------------------------------
subgroup_lookup <- read.csv(here("data", "cognate", "initial_datasets", "COGNATE_subgroup_lookup.csv"))
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

points_coloured <- COGNATE_final |>
  dplyr::left_join(dplyr::select(subgroup_lookup, glottocode, subgroup), by = "glottocode")
stopifnot(
  "Some languages have no subgroup colour — rerun [0]_Phylogenetic_Tree.R to refresh the lookup." =
    !anyNA(points_coloured$subgroup)
)

base_plot <- ggplot() +
  geom_tile(data = cog_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  # Language points, coloured by subgroup. Three concentric layers: a thin dark
  # outer rim, a white halo, then the subgroup colour at size 4 — identical
  # construction to the phoneme/grammar/genetic FEEMS maps. The halo lifts the
  # dot off the orange/cyan raster; the rim keeps the edge visible where the
  # raster passes through white (log ratio ~ 0), which a white-only ring would
  # disappear into.
  geom_point(data = points_coloured, aes(x = longitude, y = latitude),
             size = 4.8, colour = "grey20") +
  geom_point(data = points_coloured, aes(x = longitude, y = latitude),
             size = 4.4, colour = "white") +
  geom_point(data = points_coloured, aes(x = longitude, y = latitude, colour = subgroup),
             size = 4.0) +
  scale_colour_manual(values = pal, guide = "none") +
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

print(base_plot)
saveRDS(base_plot, file = here("data", "cognate", "feems", "base_plot_cognate_FEEMS.rds"))
ggsave(here("figures", "cognate", "feems", "cognate_feems_surface.png"), base_plot,
       width = 7, height = 9, units = "in", dpi = 300, bg = "white")

# ---- Null-model base map (same points, same shared color scale) -------------
# Geographic-shuffle null: average FEEMS surface across 100 permutations that
# shuffle which language sits at which location, from cognate_feems.ipynb's
# permutation cells. Same point overlay and `lims` as the real map above.
base_plot_null <- ggplot() +
  geom_tile(data = cog_surface_null, aes(x = lon, y = lat, fill = log_w_ratio_null_mean), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = points_coloured, aes(x = longitude, y = latitude),
             size = 4.8, colour = "grey20") +
  geom_point(data = points_coloured, aes(x = longitude, y = latitude),
             size = 4.4, colour = "white") +
  geom_point(data = points_coloured, aes(x = longitude, y = latitude, colour = subgroup),
             size = 4.0) +
  scale_colour_manual(values = pal, guide = "none") +
  guides(
    fill = guide_colorbar(title = expression(log[10](w/bar(w))), title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    values = rescale(c(-vmax, 0, vmax), from = lims),  # same lims as the real plot
    limits = lims,
    na.value = "transparent"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())

print(base_plot_null)
saveRDS(base_plot_null, file = here("data", "cognate", "feems", "base_plot_cognate_FEEMS_null.rds"))
ggsave(here("figures", "cognate", "feems", "cognate_feems_surface_null.png"), base_plot_null,
       width = 7, height = 9, units = "in", dpi = 300, bg = "white")

# ---- Legend strip (identical construction to phoneme/grammar/genetic) -------
# Hand-built stacked colour tiles, ordered by mean latitude (north first) so
# the strip reads top-to-bottom the way the subgroups sit on the map.
subgroup_lat <- points_coloured %>%
  dplyr::group_by(subgroup) %>%
  dplyr::summarise(mean_lat = mean(latitude), .groups = "drop")

legend_df <- subgroup_lookup %>%
  dplyr::distinct(subgroup, colour) %>%
  dplyr::left_join(subgroup_lat, by = "subgroup") %>%
  dplyr::arrange(dplyr::desc(mean_lat)) %>%
  dplyr::mutate(row = dplyr::row_number())

legend_strip <- ggplot(legend_df) +
  geom_tile(aes(x = 0, y = -row, fill = subgroup), width = 0.9, height = 0.8) +
  geom_text(aes(x = 0.65, y = -row, label = subgroup),
            hjust = 0, size = 2.8, colour = "grey15") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_x_continuous(limits = c(-0.5, 6.5), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(6, 2, 6, 2))

print(legend_strip)
saveRDS(legend_strip, file = here("data", "cognate", "feems", "legend_strip_cognate_subgroup.rds"))
