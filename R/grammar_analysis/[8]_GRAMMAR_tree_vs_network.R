# =============================================================================
# [8] Grammar Analysis — EXPERIMENTAL: waypoint network vs. phylogeny
#
# Diagnostic, NOT part of the core pipeline. Puts two panels side by side to
# eyeball a single hypothesis: has the waypoint / relative-migration network
# essentially recreated the language phylogeny?
#
#   Left  — figures/shared/grammar_phylogenetic_tree.png (from [0]); tips coloured by
#           level-4 Glottolog family subgroup.
#   Right — [3]'s plot_network() map, with the language points coloured by the
#           SAME subgroup palette (they are uncoloured shape-21 dots in [3]).
#
# If the network recreates the tree, same-coloured points should cluster
# geographically on the right the way same-coloured tips cluster on the left.
#
# DEPENDENCY / RUN ORDER:
#   - REQUIRES [0]_Phylogenetic_Tree.R to have been run: needs
#     figures/shared/grammar_phylogenetic_tree.png AND data/GRAMMAR_subgroup_lookup.csv
#     (the shared subgroup -> colour map, so both panels match exactly).
#   - source()s [3]_GRAMMAR_network_distance.R to rebuild the network objects.
# Output: figures/shared/grammar_tree_vs_network.png (paired comparison)
#         figures/grammar/mst_waypoints/GRAMMAR_network_by_subgroup.png (network
#         panel alone, with a tree-style stacked colour-tile legend, titled for
#         the Manila Galleon trade-route framing) — filed under grammar/, not
#         shared/, since it's specific to the grammar network (mirrors the
#         PHONEME_ version under figures/phoneme/mst_waypoints/).
# =============================================================================

library(tidyverse)
library(sf)
library(ggplot2)
library(patchwork)
library(png)
library(grid)
library(here)

# ---- 1. Network objects from [3] --------------------------------------------
# [3] is self-contained (reads data/) and, as a side effect of sourcing, leaves
# full_tree_sf, arrow_connectors, GRAMMAR_cossim, world_map, MANILA
# and plot_network() in the environment. NOTE: [3] attaches library(maps), which
# masks purrr::map — this script uses only dplyr verbs, so that is harmless here.
source(here("R", "grammar_analysis", "[3]_GRAMMAR_network_distance.R"), echo = FALSE)

# ---- 2. Shared subgroup palette (identical to the tree) ----------------------
subgroup_lookup <- read_csv(here("data", "GRAMMAR_subgroup_lookup.csv"),
                            show_col_types = FALSE)
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

# Attach subgroup to the language points. Unlike phoneme (where all 58 network
# languages have a tree tip), several grammar Philippine languages have NO tree
# placement (the ABVD tree lacks Sabah/Sama-Bajau coverage — see
# [0]_Phylogenetic_Tree.R's pruning report, 13 of 39 dropped) and so have no
# subgroup colour. Those points are dropped from this diagnostic (not from the
# core pipeline) rather than erroring, with a message reporting how many.
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
# plot_network()'s body from [3], simplified for this comparison: the language
# points are shape-21 discs filled by subgroup (grey outline for contrast on the
# land polygons), and the connector arrows are drawn a uniform light grey, a
# little thinner than the main network edges, instead of [3]'s black/firebrick
# "H2 route" colouring. No legends: the subgroup fill is suppressed (the left
# panel's labelled colour strip is the shared key) and the connectors no longer
# carry a colour scale, so the map stays uncluttered. Main network edges have no
# arrowheads — they're a bidirectional route between nodes, so a directed arrow
# on one edge conflicts visually with the arrow its reverse would need; only the
# (directional) language -> node connectors keep arrowheads.
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
    labs(title = title, x = "Longitude", y = "Latitude")
}

network_panel <- plot_network_coloured(full_tree_sf, arrow_connectors,
                                       points_coloured, refdf1 = MANILA)

# ---- 4. Left panel: the phylogeny PNG as a raster ---------------------------
tree_png  <- png::readPNG(here("figures", "shared", "phylogenetic_tree.png"))
tree_grob <- grid::rasterGrob(tree_png, interpolate = TRUE)

# ---- 5. Combine + save ------------------------------------------------------
combined <- (patchwork::wrap_elements(full = tree_grob) | network_panel) +
  plot_layout(widths = c(0.92, 1)) +
  plot_annotation(
    title = "Does the waypoint network recreate the phylogeny?",
    subtitle = "Left: ABVD phylogeny. Right: relative-migration network. Points/tips share the family-subgroup palette.",
    theme = theme(plot.title = element_text(face = "bold"))
  )
print(combined)

ggsave(here("figures", "shared", "grammar_tree_vs_network.png"),
       combined, width = 15, height = 9, units = "in", dpi = 300)

# ---- 6. Standalone network figure, with a tree-style legend strip -----------
# The right-hand panel on its own (no phylogeny alongside), for use outside the
# tree-vs-network comparison. It needs its own legend since the tree's colour
# strip isn't present to serve as the shared key here — built to match that
# strip's look (a stacked colour tile + name per subgroup) rather than
# ggplot's default small-swatch legend, so the two figures read as one family.
# Rows are ordered by each subgroup's mean latitude (north first), not by
# clade size, so the legend reads top-to-bottom the same way the subgroups sit
# on the map (e.g. the northern Cordillera cluster near the top of the legend,
# the Mindanao clusters near the bottom).
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

# figures/grammar/mst_waypoints/, not figures/shared/: this network is specific
# to the grammar cosine-similarity dataset ([3]_GRAMMAR_network_distance.R),
# matching the existing PHONEME_/GRAMMAR_ split for the other waypoint figures
# in these two folders.
ggsave(here("figures", "grammar", "mst_waypoints", "GRAMMAR_network_by_subgroup.png"),
       network_with_legend, width = 11, height = 9, units = "in", dpi = 300)
