# =============================================================================
# [8] Grammar Analysis — EXPERIMENTAL: waypoint network vs. phylogeny
#
# Diagnostic, NOT part of the core pipeline. Pairs the [0] phylogeny PNG with
# [3]'s waypoint-network map, both coloured by the same level-4 Glottolog
# subgroup palette, to eyeball whether the network recreates the phylogeny, and
# saves the network panel separately with a matching legend strip.
#
# Run order: requires [0]_Phylogenetic_Tree.R (for the tree PNG and the shared
# subgroup lookup) and source()s [3]_GRAMMAR_network_distance.R for the network
# objects.
#
# Inputs:  figures/phylogenetic_tree/grammar_phylogenetic_tree.png, data/GRAMMAR_subgroup_lookup.csv
# Outputs: figures/tree_vs_network/grammar_tree_vs_network.png,
#          figures/mst_waypoints/GRAMMAR_network_by_subgroup.png
# =============================================================================

library(tidyverse)
library(sf)
library(ggplot2)
library(patchwork)
library(png)
library(grid)
library(here)

# ---- 1. Network objects from [3] --------------------------------------------
# Sourcing [3] leaves full_tree_sf, arrow_connectors, GRAMMAR_cossim, world_map,
# MANILA and plot_network() in the environment. NOTE: [3] attaches library(maps),
# which masks purrr::map — harmless here, since only dplyr verbs are used.
source(here("R", "grammar_analysis", "[3]_GRAMMAR_network_distance.R"), echo = FALSE)

# ---- 2. Shared subgroup palette (identical to the tree) ----------------------
subgroup_lookup <- read_csv(here("data", "grammar", "GRAMMAR_subgroup_lookup.csv"),
                            show_col_types = FALSE)
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

# Attach subgroup to the language points. Philippine languages the ABVD tree
# cannot place (no Sabah/Sama-Bajau coverage) have no subgroup colour, so they
# are dropped from this diagnostic only — the core pipeline still keeps them.
points_coloured <- GRAMMAR_cossim %>%
  dplyr::left_join(dplyr::select(subgroup_lookup, language, subgroup),
                   by = "language")

n_no_subgroup <- sum(is.na(points_coloured$subgroup))
if (n_no_subgroup > 0) {
  message(
    n_no_subgroup, " of ", nrow(points_coloured),
    " network languages have no tree placement / subgroup colour and are ",
    "dropped from this diagnostic: ",
    paste(points_coloured$language[is.na(points_coloured$subgroup)], collapse = ", ")
  )
  points_coloured <- points_coloured %>% dplyr::filter(!is.na(subgroup))
}

# ---- 3. Coloured network panel ----------------------------------------------
# [3]'s plot_network() with the points filled by subgroup and the connectors a
# uniform grey instead of the "H2 route" colouring. Legends are suppressed: the
# tree panel's colour strip is the shared key for both panels.
plot_network_coloured <- function(full_tree_sf, arrow_connectors,
                                  points_df, refdf1 = MANILA,
                                  lon_range = c(116, 127), lat_range = c(4, 21),
                                  title = "Waypoint network: points coloured by family subgroup") {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70") +
    geom_sf(data = full_tree_sf %>% dplyr::filter(source == "main"),
            linewidth = 1, color = "grey40") +
    geom_sf(data = arrow_connectors,
            color = "grey60", linewidth = 0.5,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_point(data = points_df,
               aes(x = longitude, y = latitude, fill = subgroup),
               size = 3, shape = 21, colour = "grey20") +
    geom_point(data = refdf1, aes(x = longitude, y = latitude),
               shape = 23, size = 3, fill = "white", stroke = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    coord_sf(xlim = lon_range, ylim = lat_range) +
    theme_minimal() +
    # Coordinates carry no finding; the coastline is the only reference frame the
    # reader needs, so the axes, ticks and graticule are all dropped.
    theme(panel.grid = element_blank(),
          axis.text  = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()) +
    labs(title = title)
}

network_panel <- plot_network_coloured(full_tree_sf, arrow_connectors,
                                       points_coloured, refdf1 = MANILA)

# ---- 4. Left panel: the phylogeny PNG as a raster ---------------------------
tree_png  <- png::readPNG(here("figures", "phylogenetic_tree", "grammar_phylogenetic_tree.png"))
tree_grob <- grid::rasterGrob(tree_png, interpolate = TRUE)

# ---- 5. Combine + save ------------------------------------------------------
# The two output folders are method-scoped, so neither is guaranteed to exist yet:
# [0] only creates figures/phylogenetic_tree/.
dir.create(here("figures", "tree_vs_network"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures", "mst_waypoints"), recursive = TRUE, showWarnings = FALSE)
combined <- (patchwork::wrap_elements(full = tree_grob) | network_panel) +
  plot_layout(widths = c(0.92, 1)) +
  plot_annotation(
    title = "Does the waypoint network recreate the phylogeny?",
    subtitle = "Left: ABVD phylogeny. Right: relative-migration network. Points/tips share the family-subgroup palette.",
    theme = theme(plot.title = element_text(face = "bold"))
  )
print(combined)

ggsave(here("figures", "tree_vs_network", "grammar_tree_vs_network.png"),
       combined, width = 15, height = 9, units = "in", dpi = 300)

# ---- 6. Standalone network figure, with a tree-style legend strip -----------
# The network panel alone needs its own key, built as a stacked colour tile per
# subgroup to match the tree's strip rather than ggplot's default swatches.
# Rows are ordered by mean latitude (north first) so the legend reads
# top-to-bottom the way the subgroups sit on the map.
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

# title = NULL drops plot_network_coloured()'s own internal title so only the
# outer plot_annotation() title below is shown (avoids a double heading).
network_standalone <- plot_network_coloured(full_tree_sf, arrow_connectors,
                                            points_coloured,
                                            refdf1 = MANILA, title = NULL)

network_with_legend <- (network_standalone | legend_strip) +
  plot_layout(widths = c(5, 1.2)) +
  plot_annotation(
    title = "Geographical Distribution of Languages along the Manila Galleon Trade Route",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  )
print(network_with_legend)

# Filed under figures/mst_waypoints/ with the other waypoint figures, and named
# GRAMMAR_ rather than left generic: this network is specific to the grammar dataset.
ggsave(here("figures", "mst_waypoints", "GRAMMAR_network_by_subgroup.png"),
       network_with_legend, width = 11, height = 9, units = "in", dpi = 300)
