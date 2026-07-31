# =============================================================================
# [6] Phoneme Analysis — FEEMS migration-surface base map
# Renders the FEEMS effective-migration surface (the log10(w/w-bar) raster from
# python/phoneme_feems.ipynb) as a ggplot base map of the Philippines, with each
# language's subgroup membership overlaid as points (same palette/legend as
# [8]), plus a matching legend strip. Saved as RDS that [7] loads to overlay
# the waypoint routes.
#
# Input:   data/phoneme_surface_raster.csv (FEEMS surface, from phoneme_feems.ipynb),
#          data/PHONEME_cossim.csv (per-language coordinates, from [1]),
#          data/PHONEME_subgroup_lookup.csv (subgroup + colour, from [0]_Phylogenetic_Tree.R)
# Outputs: data/base_plot_phoneme_FEEMS.rds, data/legend_strip_phoneme_subgroup.rds
# Next:    [7]_PA_geoplots.R
# =============================================================================

library(ggplot2)
library(dplyr)
library(here)
library(scales)

phon_surface   <- read.csv(here("data", "phoneme_surface_raster.csv"))
PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))

vmax <- max(abs(phon_surface$log_w_ratio), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

# ---- Subgroup palette + membership (shared with [8]) ------------------------
subgroup_lookup <- read.csv(here("data", "PHONEME_subgroup_lookup.csv"))
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

points_coloured <- PHONEME_cossim |>
  dplyr::left_join(dplyr::select(subgroup_lookup, language, subgroup), by = "language")
stopifnot(
  "Some network languages have no subgroup colour — rerun [0] to refresh the lookup." =
    !anyNA(points_coloured$subgroup)
)

base_plot <- ggplot() +
  geom_tile(data = phon_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
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
saveRDS(base_plot, file = here("data", "base_plot_phoneme_FEEMS.rds"))

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
saveRDS(legend_strip, file = here("data", "legend_strip_phoneme_subgroup.rds"))
