# =============================================================================
# [6] Grammar Analysis — FEEMS migration-surface base map
# Renders the FEEMS effective-migration surface (the log10(w/w-bar) raster from
# python/grammar_feems.ipynb) as a ggplot base map of the Philippines, with each
# language's subgroup membership overlaid as points (same palette/legend as
# [8]), plus a matching legend strip. Also renders the geographic-shuffle null
# model's averaged surface (same notebook) as a second base map, on the SAME
# color scale as the real one, so the two are directly comparable. Both saved
# as RDS that [7] loads to overlay the waypoint routes.
#
# Input:   data/grammar_surface_raster.csv, data/grammar_surface_raster_null_mean.csv
#          (FEEMS surface + null-model average, from grammar_feems.ipynb),
#          data/GRAMMAR_cossim.csv (per-language coordinates, from [1]),
#          data/GRAMMAR_subgroup_lookup.csv (subgroup + colour, from [0]_Phylogenetic_Tree.R)
# Outputs: data/base_plot_grammar_FEEMS.rds, data/base_plot_grammar_FEEMS_null.rds,
#          data/legend_strip_grammar_subgroup.rds
# Next:    [7]_GA_weight_mst_feems_span.R
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)
library(scales)

gram_surface      <- read.csv(here("data", "grammar_surface_raster.csv"))
gram_surface_null <- read.csv(here("data", "grammar_surface_raster_null_mean.csv"))
GRAMMAR_cossim    <- read.csv(here("data", "GRAMMAR_cossim.csv"))

# Only the subgroup-coloured markers are restricted to the 26 tree-pruned
# languages, matching [7] and [8]; the raster still comes from all 39, since
# dropping the southern Bornean samples would leave that area interpolated.
subgroup_lookup <- read.csv(here("data", "GRAMMAR_subgroup_lookup.csv"))
GRAMMAR_cossim   <- GRAMMAR_cossim[GRAMMAR_cossim$language %in% subgroup_lookup$language, ]

# Shared color scale across the real and null maps (not each normalized to its own
# max) — this is what makes a flatter null surface visibly wash out next to the real
# one; it also means the real plot's scale below is now anchored to this combined
# vmax, not gram_surface's own max as in earlier versions of this script.
vmax <- max(abs(gram_surface$log_w_ratio), abs(gram_surface_null$log_w_ratio_null_mean), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

# ---- Subgroup palette + membership (shared with [8]) ------------------------
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

points_coloured <- GRAMMAR_cossim |>
  dplyr::left_join(dplyr::select(subgroup_lookup, language, subgroup), by = "language")
stopifnot(
  "Some network languages have no subgroup colour — rerun [0] to refresh the lookup." =
    !anyNA(points_coloured$subgroup)
)

base_plot <- ggplot() +
  geom_tile(data = gram_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  # Language points, coloured by subgroup. Three concentric layers: a thin dark
  # outer rim, a white halo, then the subgroup colour at size 4 to match [1]'s
  # cosine points. The halo lifts the dot off the orange/cyan raster; the rim
  # keeps the edge visible where the raster passes through white (log ratio ~ 0),
  # which a white-only ring would disappear into. shape defaults to 19 so
  # `colour` fills the disc — `fill` is spoken for by the raster.
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

base_plot
saveRDS(base_plot, file = here("data", "base_plot_grammar_FEEMS.rds"))

# ---- Null-model base map (same points, same shared color scale) -------------
# Geographic-shuffle null: average FEEMS surface across 100 permutations that
# shuffle which language sits at which location, from grammar_feems.ipynb. Same
# point overlay and `lims` as the real map above, so [7] can render it as a
# directly comparable figure.
base_plot_null <- ggplot() +
  geom_tile(data = gram_surface_null, aes(x = lon, y = lat, fill = log_w_ratio_null_mean), alpha = 0.6) +
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

base_plot_null
saveRDS(base_plot_null, file = here("data", "base_plot_grammar_FEEMS_null.rds"))

# ---- Legend strip (identical construction to [8]) ---------------------------
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

legend_strip
saveRDS(legend_strip, file = here("data", "legend_strip_grammar_subgroup.rds"))
