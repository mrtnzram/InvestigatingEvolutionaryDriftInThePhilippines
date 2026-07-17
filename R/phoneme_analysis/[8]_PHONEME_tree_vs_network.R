# =============================================================================
# [8] Phoneme Analysis — EXPERIMENTAL: waypoint network vs. phylogeny
#
# Diagnostic, NOT part of the core pipeline. Puts two panels side by side to
# eyeball a single hypothesis: has the waypoint / relative-migration network
# essentially recreated the language phylogeny?
#
#   Left  — figures/shared/phylogenetic_tree.png (from [0]); tips coloured by
#           level-4 Glottolog family subgroup.
#   Right — [3]'s plot_network() map, with the language points coloured by the
#           SAME subgroup palette (they are uncoloured shape-21 dots in [3]).
#
# If the network recreates the tree, same-coloured points should cluster
# geographically on the right the way same-coloured tips cluster on the left.
#
# DEPENDENCY / RUN ORDER:
#   - REQUIRES [0]_Phylogenetic_Tree.R to have been run: needs
#     figures/shared/phylogenetic_tree.png AND data/PHONEME_subgroup_lookup.csv
#     (the shared subgroup -> colour map, so both panels match exactly).
#   - source()s [3]_PHONEME_network_distance.R to rebuild the network objects.
# Output: figures/shared/tree_vs_network.png
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
# full_tree_sf, arrow_main, arrow_connectors, PHONEME_cossim, world_map, MANILA
# and plot_network() in the environment. NOTE: [3] attaches library(maps), which
# masks purrr::map — this script uses only dplyr verbs, so that is harmless here.
source(here("R", "phoneme_analysis", "[3]_PHONEME_network_distance.R"), echo = FALSE)

# ---- 2. Shared subgroup palette (identical to the tree) ----------------------
subgroup_lookup <- read_csv(here("data", "PHONEME_subgroup_lookup.csv"),
                            show_col_types = FALSE)
pal <- setNames(subgroup_lookup$colour, subgroup_lookup$subgroup)

# Attach subgroup to the language points. All 58 network languages are in the
# lookup, so every point gets a colour; stop loudly if that ever changes.
points_coloured <- PHONEME_cossim %>%
  dplyr::left_join(dplyr::select(subgroup_lookup, language, subgroup),
                   by = "language")
stopifnot(
  "Some network languages have no subgroup colour — rerun [0] to refresh the lookup." =
    !anyNA(points_coloured$subgroup)
)

# ---- 3. Coloured network panel ----------------------------------------------
# plot_network()'s body from [3], simplified for this comparison: the language
# points are shape-21 discs filled by subgroup (grey outline for contrast on the
# land polygons), and the connector arrows are drawn a uniform light grey, a
# little thinner than the main network edges, instead of [3]'s black/firebrick
# "H2 route" colouring. No legends: the subgroup fill is suppressed (the left
# panel's labelled colour strip is the shared key) and the connectors no longer
# carry a colour scale, so the map stays uncluttered.
plot_network_coloured <- function(full_tree_sf, arrow_main, arrow_connectors,
                                  points_df, refdf1 = MANILA,
                                  lon_range = c(116, 127), lat_range = c(4, 21)) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70") +
    geom_sf(data = full_tree_sf %>% dplyr::filter(source == "main"),
            linewidth = 1, color = "grey40") +
    geom_sf(data = arrow_connectors,
            color = "grey60", linewidth = 0.5,
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_sf(data = arrow_main,
            arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
            color = "grey40") +
    geom_point(data = points_df,
               aes(x = longitude, y = latitude, fill = subgroup),
               size = 3, shape = 21, colour = "grey20") +
    geom_point(data = refdf1, aes(x = longitude, y = latitude),
               shape = 23, size = 3, fill = "white", stroke = 1) +
    scale_fill_manual(values = pal, guide = "none") +
    coord_sf(xlim = lon_range, ylim = lat_range) +
    theme_minimal() +
    labs(title = "Waypoint network: points coloured by family subgroup",
         x = "Longitude", y = "Latitude")
}

network_panel <- plot_network_coloured(full_tree_sf, arrow_main,
                                       arrow_connectors, points_coloured,
                                       refdf1 = MANILA)

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

ggsave(here("figures", "shared", "tree_vs_network.png"),
       combined, width = 15, height = 9, units = "in", dpi = 300)
