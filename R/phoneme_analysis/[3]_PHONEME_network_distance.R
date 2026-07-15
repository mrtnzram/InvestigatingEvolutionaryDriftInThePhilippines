# ── [3]_PHONEME_network_distance.R ───────────────────────────────────────────
# Graph network creation + per-language distance calculation.
# Input:   data/PHONEME_cossim.csv, data/nodes.csv, data/edges.csv
# Outputs: data/PHONEME_cossim_dist.csv (cossim + geodist_H1_span, for the
#          regression file), data/phoneme_waypoint_plot.rds (overview arrow plot)
# Note:    pairwise inter-language distances for the Mantel test are computed
#          separately in [5]_PHONEME_mantel.R (Dijkstra routing).
#
# compute_shortest_path_df(). Instead of routing each language to a
# single reference point (ref_coords1) via Dijkstra, this computes each
# language's terrain-penalized cost to ENTER the navigable network — i.e.
# the connector from the language's coordinate to its nearest graph node.
#
# Dropped vs. the old version:
#   - ref_coords1 / find-nearest-node-to-Manila
#   - shortest_path_trace() (Dijkstra) — no longer needed; "distance to the
#     network" terminates the moment a language reaches any node
#   - per-row plot_path() — replaced by one overview arrow plot
#
# Kept:
#   - land/sea terrain penalty (land_penalty), now factored into one helper
#     used by both the network edges AND the language connectors, removing
#     the 3x duplicated land/sea-split logic from the original script
#   - the final arrow plot (full_tree_sf, main vs. connector, directional
#     arrowheads), saved as an .rds exactly as before
# ──────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(geosphere)
library(sf)
library(maps)
library(sfheaders)
library(here)

# ── 0. Load data ─────────────────────────────────────────────────────────────
PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim_marked.csv"))

nodes <- read.csv(here("data", "nodes.csv")) |> mutate(id = as.character(id))
edges <- read.csv(here("data", "edges.csv")) |>
  mutate(from = as.character(from), to = as.character(to))


# ── 1. Land mask (Philippines + Malaysia) ────────────────────────────────────
world_map <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

land_sf <- sf_polygon(obj = world_map, polygon_id = "group", x = "long", y = "lat") |>
  st_union() |>
  st_sf(geometry = _) |>
  st_set_crs(4326)


# ── 2. Shared helper: split a linestring into land/sea + terrain-penalize ────
# Used for BOTH network edges and language→node connectors, so the land/sea
# logic only exists once instead of being copy-pasted per use case.
split_and_penalize <- function(geom, land_sf, land_penalty) {
  
  land_part <- st_intersection(geom, st_geometry(land_sf))
  sea_part  <- st_difference(geom, st_geometry(land_sf))
  
  # st_intersection/st_difference at the bare sfc level RETAIN empty results
  # (unlike the sf data.frame method, which drops them) — filter explicitly
  # rather than relying on length() to tell empty from non-empty.
  land_part <- land_part[!st_is_empty(land_part)]
  sea_part  <- sea_part[!st_is_empty(sea_part)]
  
  land_len <- if (length(land_part) > 0) as.numeric(sum(st_length(land_part))) else 0
  sea_len  <- if (length(sea_part)  > 0) as.numeric(sum(st_length(sea_part)))  else 0
  
  list(
    land_len      = land_len,
    sea_len       = sea_len,
    crosses_land  = land_len > 0,
    weighted_cost = land_len * land_penalty + sea_len,
    land_geom     = if (length(land_part) > 0) land_part else NULL,  # sfc (may hold >1 piece)
    sea_geom      = if (length(sea_part)  > 0) sea_part  else NULL   # sfc (may hold >1 piece)
  )
}


# ── 3. Build the network edges (the "main" path) ─────────────────────────────
# NOTE: geometry is built and converted to a real sfc column (via st_as_sf)
# in its own stage BEFORE the land/sea split. Under rowwise(), a column
# built inline as list(st_linestring(...)) is already unwrapped to the bare
# sfg when referenced later in the SAME mutate() call — so `geometry[[1]]`
# would index into the geometry's underlying matrix instead of extracting
# it, handing split_and_penalize() a raw numeric. Once geometry is a proper
# sfc (post st_as_sf), a second rowwise() stage extracts each row's geometry
# correctly.
build_network_edges <- function(nodes, edges, land_sf, land_penalty) {
  
  edges_bi <- bind_rows(edges, edges |> rename(from = to, to = from)) |> distinct()
  
  edges_sf <- edges_bi |>
    rowwise() |>
    mutate(geometry = list(st_linestring(matrix(c(
      nodes$longitude[nodes$id == from], nodes$latitude[nodes$id == from],
      nodes$longitude[nodes$id == to],   nodes$latitude[nodes$id == to]
    ), ncol = 2, byrow = TRUE)))) |>
    ungroup() |>
    st_as_sf(crs = 4326)
  
  edges_sf |>
    rowwise() |>
    mutate(
      .split        = list(split_and_penalize(geometry, land_sf, land_penalty)),
      land_len      = .split$land_len,
      sea_len       = .split$sea_len,
      crosses_land  = .split$crosses_land,
      weighted_cost = .split$weighted_cost,
      land_geom     = list(.split$land_geom),
      sea_geom      = list(.split$sea_geom)
    ) |>
    ungroup() |>
    select(-.split) |>
    st_as_sf()
}


# ── 4. Nearest-node finder ────────────────────────────────────────────────────
find_nearest_node <- function(coords, nodes) {
  d <- distHaversine(matrix(c(nodes$longitude, nodes$latitude), ncol = 2), coords)
  nodes$id[which.min(d)]
}


# ── Manila reference (hoisted so plots can use it too) ──────────────────────
MANILA <- data.frame(
  longitude = 121,
  latitude  = 14.6
)


# ── 5. Per-language: terrain-penalized distance to the network ───────────────
# geodist_H1_span : penalized cost from the language point to its nearest
#                   network node (the "get onto the network" leg).
# geodist_H2_span : penalized cost to Manila, taking whichever is cheaper —
#                   (connector + network traversal) or a direct line.
# using_network   : TRUE when the network route won; also selects which
#                   split geometry is retained for plotting.
compute_network_distance_df <- function(df, nodes, edges, land_sf,
                                        refdf1 = MANILA,
                                        land_penalty = 4.44) {
  
  stopifnot(land_penalty >= 1)
  
  nodes <- nodes |>
    dplyr::mutate(id = as.character(id))
  
  edges <- edges |>
    dplyr::mutate(
      from = as.character(from),
      to   = as.character(to)
    )
  
  #------------------------------------------------------------
  # Build weighted graph
  #------------------------------------------------------------
  g <- igraph::graph_from_data_frame(
    d = edges |>
      dplyr::select(from, to, weighted_cost),
    directed = FALSE,
    vertices = nodes
  )
  
  #------------------------------------------------------------
  # Compute shortest weighted distance from Manila node to all nodes
  #------------------------------------------------------------
  manila_coords <- c(refdf1$longitude, refdf1$latitude)
  
  manila_node <- find_nearest_node(manila_coords, nodes)
  
  net_to_manila <- igraph::distances(
    g,
    v = manila_node,
    to = igraph::V(g),
    weights = igraph::E(g)$weighted_cost
  )
  
  net_dist_km <- tibble(
    nearest_node = colnames(net_to_manila),
    net_km = as.numeric(net_to_manila[1, ]) / 1000
  )
  
  #------------------------------------------------------------
  # Build connector to nearest node
  #------------------------------------------------------------
  df_geom <- df |>
    rowwise() |>
    mutate(
      start_coords = list(c(longitude, latitude)),
      nearest_node = find_nearest_node(start_coords, nodes),
      
      connector_geom = list(
        st_linestring(rbind(
          start_coords,
          c(
            nodes$longitude[nodes$id == nearest_node],
            nodes$latitude[nodes$id == nearest_node]
          )
        ))
      )
    ) |>
    ungroup() |>
    mutate(
      connector_geom = st_sfc(connector_geom, crs = 4326)
    )
  
  #------------------------------------------------------------
  # Build direct connector to Manila
  #------------------------------------------------------------
  df_geom <- df_geom |>
    rowwise() |>
    mutate(
      manila_geom = list(
        st_linestring(rbind(
          start_coords,
          manila_coords
        ))
      )
    ) |>
    ungroup() |>
    mutate(
      manila_geom = st_sfc(manila_geom, crs = 4326)
    )
  
  #------------------------------------------------------------
  # Compute connector costs, pick the cheaper H2 route, and retain
  # the geometry belonging to whichever route won.
  #------------------------------------------------------------
  df_geom |>
    rowwise() |>
    mutate(
      # Point -> nearest node
      .split_node = list(
        split_and_penalize(connector_geom, land_sf, land_penalty)
      ),
      
      # Point -> Manila directly
      .split_manila = list(
        split_and_penalize(manila_geom, land_sf, land_penalty)
      ),
      
      geodist_H1_span  = .split_node$weighted_cost   / 1000,   # m -> km
      direct_to_manila = .split_manila$weighted_cost / 1000    # m -> km
    ) |>
    ungroup() |>
    left_join(net_dist_km, by = join_by(nearest_node)) |>
    mutate(
      network_to_manila = geodist_H1_span + net_km,
      
      # igraph::distances() yields Inf for unreachable nodes, so a
      # disconnected graph falls back to the direct line rather than NA.
      using_network = !is.na(network_to_manila) &
        network_to_manila < direct_to_manila,
      
      geodist_H2_span = pmin(network_to_manila, direct_to_manila, na.rm = TRUE)
    ) |>
    # Geometry follows the winning route: these feed connector_sf in §7.
    rowwise() |>
    mutate(
      land_geom    = list(if (using_network) .split_node$land_geom else .split_manila$land_geom),
      sea_geom     = list(if (using_network) .split_node$sea_geom  else .split_manila$sea_geom),
      land_len     = if (using_network) .split_node$land_len     else .split_manila$land_len,
      sea_len      = if (using_network) .split_node$sea_len      else .split_manila$sea_len,
      crosses_land = if (using_network) .split_node$crosses_land else .split_manila$crosses_land
    ) |>
    ungroup() |>
    dplyr::select(
      -.split_node, -.split_manila,
      -net_km, -network_to_manila, -direct_to_manila
    )
}


# ── 6. Run ─────────────────────────────────────────────────────────────────
network_edges <- build_network_edges(
  nodes,
  edges,
  land_sf,
  land_penalty = 4.44
)

PHONEME_cossim <- compute_network_distance_df(
  PHONEME_cossim,
  nodes,
  network_edges,
  land_sf,
  refdf1 = MANILA,
  land_penalty = 4.44
)

PHONEME_cossim |>
  select(any_of(c("language", "nearest_node", "geodist_H1_span",
                  "land_len", "sea_len", "crosses_land"))) |>
  arrange(geodist_H1_span) |>
  print(n = 20)

PHONEME_cossim |>
  select(any_of(c("language", "nearest_node", "geodist_H2_span",
                  "using_network", "land_len", "sea_len", "crosses_land"))) |>
  arrange(geodist_H2_span) |>
  print(n = 20)

# Diagnostic: how often does the network actually beat the direct line?
# At land_penalty = 4.44 expect more FALSE than under the old 44.18 setting.
PHONEME_cossim |>
  count(using_network) |>
  print()

# Write the geodist-augmented table (flat columns only; drops the geometry /
# list-columns) for [4]_PHONEME_regression.R to consume. The *_influenced /
# sig_* columns are carried forward when present (written by [2]); any_of() keeps
# this runnable if [2] hasn't been sourced yet.
PHONEME_cossim |>
  select(language, latitude, longitude, starts_with("cossim_"),
         any_of(c("span_influenced", "jap_influenced", "eng_influenced")),
         geodist_H1_span, geodist_H2_span, using_network) |>
  write.csv(file = here("data", "PHONEME_cossim_dist.csv"), row.names = FALSE)


# ── 7. Assemble full_tree_sf (main edges + per-language connectors) ─────────
# Mirrors the original land_segments_sf / sea_segments_sf / connector_sf
# assembly, but built once via map()/list_rbind() instead of manual loops,
# and with a single connector per language (no "end" connector, since there
# is no longer a ref_coords1 target). The connector geometry is now whichever
# route won in §5, so `using_network` is carried through for styling.
geom_rows <- function(geom_list, crosses_land_value, source_value,
                      row_attrs = NULL) {
  # NOTE: purrr::map is explicitly namespaced — library(maps) is loaded after
  # library(tidyverse), so the unqualified map() resolves to maps::map()
  # (cartographic plotting), not purrr::map() (list mapping), and silently
  # fails with a cryptic "database type not supported" error.
  #
  # to_sfg_list() normalizes whatever split_and_penalize() stored for a row
  # — a bare sfg, a (possibly multi-piece) sfc, or NULL — into a flat list of
  # bare sfg objects, which is the only thing st_sfc() accepts. as.list() on
  # an sfc is the one operation sf guarantees always drops correctly to sfg;
  # relying on `[[1]]` further upstream wasn't safe across this rowwise chain.
  to_sfg_list <- function(x) {
    if (is.null(x))         return(list())
    if (inherits(x, "sfc")) return(as.list(x))
    if (inherits(x, "sfg")) return(list(x))
    list()
  }
  
  unwrapped   <- purrr::map(geom_list, 1)                  # undo rowwise's list(x) wrapper
  sfg_per_row <- purrr::map(unwrapped, to_sfg_list)        # per-row list of sfg
  n_per_row   <- purrr::map_int(sfg_per_row, length)       # 0 when a row has no piece
  sfg_list    <- do.call(c, sfg_per_row)                   # flatten to one list of sfg
  
  if (length(sfg_list) == 0) return(NULL)
  
  out <- st_sf(
    geometry     = st_sfc(sfg_list, crs = 4326),
    crosses_land = crosses_land_value,
    source       = source_value
  )
  
  # rep(times = n_per_row) keeps per-language attrs aligned with the
  # flattened geometry: a row contributing zero pieces contributes zero attrs.
  if (!is.null(row_attrs)) {
    out$using_network <- rep(row_attrs, times = n_per_row)
  }
  
  out
}

main_sf <- bind_rows(
  geom_rows(edge_lines$land_geom, TRUE,  "main"),
  geom_rows(edge_lines$sea_geom,  FALSE, "main")
)

connector_sf <- bind_rows(
  geom_rows(PHONEME_cossim$land_geom, TRUE,  "connector", PHONEME_cossim$using_network),
  geom_rows(PHONEME_cossim$sea_geom,  FALSE, "connector", PHONEME_cossim$using_network)
)

# main_sf has no using_network; bind_rows fills NA, which is harmless because
# the plots filter on `source` before mapping colour.
full_tree_sf <- bind_rows(main_sf, connector_sf)


# ── 8. Simplify each main edge to a straight start→end segment, for arrows ──
# (Land/sea-split edges can be multi-vertex; arrows look cleaner on the
# overall edge direction rather than every sub-segment.)
simplify_to_straight <- function(sf_lines) {
  sf_lines |>
    mutate(.row = row_number()) |>
    group_by(.row) |>
    group_modify(~ {
      coords <- st_coordinates(.x)
      if (nrow(coords) < 2) return(.x)
      new_geom <- st_sfc(
        st_linestring(rbind(coords[1, c("X", "Y")], coords[nrow(coords), c("X", "Y")])),
        crs = st_crs(.x)
      )
      st_set_geometry(.x, new_geom)
    }) |>
    ungroup() |>
    select(-.row) |>
    st_as_sf()
}

arrow_main <- full_tree_sf |>
  filter(source == "main") |>
  simplify_to_straight()

arrow_connectors <- full_tree_sf |>
  filter(source == "connector")


# ── 9. Overview plot: network + language points + directional arrows ────────
# Connector arrows are coloured by using_network: black arrows stop at a
# network node (main edges carry them onward to Manila), red arrows run
# straight to Manila.
plot_network <- function(full_tree_sf, arrow_main, arrow_connectors,
                         points_df, refdf1 = MANILA,
                         lon_range = c(116, 127), lat_range = c(4, 21)) {
  ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                 fill = "gray95", color = "gray70") +
    geom_sf(data = full_tree_sf |> filter(source == "main"),
            linewidth = 1, color = "grey40") +
    geom_sf(data = arrow_connectors,
            aes(color = using_network),
            arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_sf(data = arrow_main,
            arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
            color = "grey40") +
    geom_point(data = points_df, aes(x = longitude, y = latitude),
               size = 3, shape = 21) +
    geom_point(data = refdf1, aes(x = longitude, y = latitude),
               shape = 23, size = 3, fill = "white", stroke = 1) +
    scale_color_manual(
      values = c(`TRUE` = "black", `FALSE` = "firebrick"),
      labels = c(`TRUE` = "via network node", `FALSE` = "direct to Manila"),
      name   = "H2 route"
    ) +
    coord_sf(xlim = lon_range, ylim = lat_range) +
    theme_minimal() +
    labs(title = "Phonemic Historical Routes Waypoint System",
         x = "Longitude", y = "Latitude")
}

print(plot_network(full_tree_sf, arrow_main, arrow_connectors,
                   PHONEME_cossim, refdf1 = MANILA))


# ── 10. Arrow-only plot (main tree + connector arrows, all black) ───────────
# Built as a transparent overlay layer for the EEMS raster + cossim dots.
arrow_plot <- ggplot() +
  geom_sf(data = full_tree_sf |> filter(source == "main"),
          color = "black", linewidth = 1) +
  geom_sf(data = arrow_connectors,
          color = "black", linewidth = 1,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background  = element_rect(fill = "transparent", color = NA)
  )

print(arrow_plot)
saveRDS(arrow_plot, file = here("data", "phoneme_waypoint_plot.rds"))
