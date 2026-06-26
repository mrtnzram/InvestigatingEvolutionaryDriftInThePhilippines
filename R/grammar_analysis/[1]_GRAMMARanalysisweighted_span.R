library(readr)
library(tidyverse)
library(dplyr)
library(geosphere)
library(rethinking)
library(lingtypology)
library(infotheo)
library(reshape2)
library(ggplot2)
library(purrr)
library(here)
library(ggridges)
library(fastDummies)


# ---- Prepping Globals --------------------------------------------------------------------------

GRAMBANKdf_PH <- read_csv(here("data", "GRAMBANKdf_PH.csv"))

GRAMBANKdf_PH <- GRAMBANKdf_PH %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

ph_lang <- c('Tiruray', 'Maranao', 'Central Sama', 'Hiligaynon', 'Aklanon', 'Filipino', 'Coastal-Naga Bikol','Tagalog','Pampanga','Pangasinan','Iloko','Dupaninan Agta',
             'Cebuano','Southern Sama')
int_lang <- c('English','Spanish','Japanese')
unr_lang <- c("Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Liv",
              "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")

GRAMBANK_query <- grambank.feature(c('gb020','gb021','gb022','gb023','gb028','gb030','gb031','gb035','gb036','gb037','gb042','gb043','gb044','gb051','gb052','gb053',
                                     'gb054','gb065','gb070','gb071','gb072','gb073','gb079','gb080','gb082','gb083','gb084','gb086','gb089','gb090','gb091','gb092',
                                     'gb093','gb094','gb107','gb121','gb130','gb131','gb137','gb138','gb171','gb172','gb186','gb192','gb196','gb197','gb316','gb318',
                                     'gb321','gb415'),na.rm = FALSE)

feature_cols <- intersect(colnames(GRAMBANKdf_PH), colnames(GRAMBANK_query)) %>%
  setdiff(c("longitude", "latitude"))

feature_cols


gramfeature_freq <- read_csv(here("data", "gramfeature_freq.csv"))

# ---- one-hot encode categorical columns ----------

# Classify feature_cols by cardinality
n_levels <- map_int(
  feature_cols,
  \(col) GRAMBANKdf_PH[[col]] |> na.omit() |> unique() |> length()
)

binary_cols      <- feature_cols[n_levels <= 2]   # 0/1  → already valid; keep
categorical_cols <- feature_cols[n_levels >  2]   # 0/1/2 → expand to 3 dummies

cat(sprintf(
  "Binary (kept):       %d cols\nCategorical (OHE'd): %d cols\n",
  length(binary_cols), length(categorical_cols)
))

# One-hot encode the 3-level columns
GRAMBANKdf_PH_ohe <- GRAMBANKdf_PH |>
  dummy_cols(
    select_columns          = categorical_cols,
    remove_selected_columns = TRUE,   # drop original integer column
    remove_first_dummy      = FALSE,  # full OHE — keep all 3 levels
    ignore_na               = TRUE    # NA row → NA across all 3 dummies
  )

# Reconstruct updated feature_cols vector
# fastDummies names dummies as <col>_<value> (e.g. CONJ_COORD_0, _1, _2)
ohe_pattern  <- str_c("^(", str_c(categorical_cols, collapse = "|"), ")_")
ohe_cols     <- str_subset(names(GRAMBANKdf_PH_ohe), ohe_pattern)

feature_cols_ohe <- c(binary_cols, ohe_cols)

cat(sprintf(
  "Feature matrix: %d binary + %d OHE cols = %d total features\n",
  length(binary_cols), length(ohe_cols), length(feature_cols_ohe)
))

GRAMBANKdf_PH <- GRAMBANKdf_PH_ohe

# ----- cosine similarity ------------------------

calculate_weighted_cosine_similarity <- function(GRAMBANKdf_PH, gramfeature_freq, feature_cols, id_col = "language", feature_col = "feature") {
  
  aligned_freq <- gramfeature_freq %>%
    filter(feature %in% feature_cols) %>%
    distinct(feature, .keep_all = TRUE) %>%
    arrange(match(feature, feature_cols))
  
  idf_weights <- aligned_freq$IDF
  binary_data <- GRAMBANKdf_PH %>%
    select(all_of(aligned_freq$feature)) %>%
    as.matrix()
  

  long_data <- GRAMBANKdf_PH %>%
    select(all_of(feature_cols)) %>%
    mutate(language = GRAMBANKdf_PH[[id_col]]) %>%
    pivot_longer(cols = -language, names_to = "feature", values_to = "value")
  
  weighted_long <- long_data %>%
    left_join(gramfeature_freq, by = c("feature", "value")) %>%
    mutate(weighted_value = IDF)
  
  weighted_data <- weighted_long %>%
    select(language, feature, weighted_value) %>%
    pivot_wider(names_from = feature, values_from = weighted_value) %>%
    column_to_rownames("language") %>%
    as.matrix()
  
  weighted_data["Spanish", ]
  
  
  
  language_ids <- GRAMBANKdf_PH[[id_col]]
  
  # Step 4: Compute cosine similarity
  n_languages <- nrow(weighted_data)
  cosine_matrix_grammar <- matrix(0, nrow = n_languages, ncol = n_languages,
                          dimnames = list(language_ids, language_ids))
  
  epsilon <- 1e-9
   
  for (i in 1:n_languages) {
    for (j in i:n_languages) {
      vec_a <- weighted_data[i, ]
      vec_b <- weighted_data[j, ]
      
      # Handle partial NAs
      valid_idx <- which(!is.na(vec_a) & !is.na(vec_b))
      
      if (length(valid_idx) == 0) {
        score <- NA
      } else {
        dot_product <- sum(vec_a[valid_idx] * vec_b[valid_idx])
        magnitude_a <- sqrt(sum(vec_a[valid_idx]^2))
        magnitude_b <- sqrt(sum(vec_b[valid_idx]^2))
        
        denominator <- magnitude_a * magnitude_b
        score <- ifelse(denominator == 0, NA, dot_product / denominator)
      }
      
      cosine_matrix_grammar[i, j] <- score
      cosine_matrix_grammar[j, i] <- score
    }
  }
  
  return(cosine_matrix_grammar)
}

feature

cosine_matrix_grammar <- calculate_weighted_cosine_similarity(
  GRAMBANKdf_PH, 
  gramfeature_freq, 
  feature_cols_ohe, 
  id_col = "language")



# INVESTIGATE --------------------------------------------------

ordered_languages <- GRAMBANKdf_PH %>%
  arrange(Language_Type) %>%
  pull(language)


cosine_matrix_grammar <- cosine_matrix_grammar[ordered_languages, ordered_languages]

setdiff(ordered_languages, rownames(cosine_matrix_grammar))


cosine_matrix_grammar['Filipino','Spanish'] 
cosine_matrix_grammar['Filipino','Japanese'] 
cosine_matrix_grammar['Filipino','English'] 

ph_lang <- GRAMBANKdf_PH %>% 
  filter(Language_Type == 'Philippine Language') %>% 
  pull(language)

ph_lang %in% rownames(cosine_matrix_grammar)


sub_matrixspan <- cosine_matrix_grammar[ph_lang, 'Spanish']
sub_matrixjap <- cosine_matrix_grammar[ph_lang, 'Japanese']
sub_matrixeng <- cosine_matrix_grammar[ph_lang, 'English']

colnames(cosine_matrix_grammar)



df_span <- as_tibble(sub_matrixspan, rownames = 'language')
colnames(df_span)[2] <- 'cossim_span'

df_jap <- as_tibble(sub_matrixjap, rownames = 'language')
colnames(df_jap)[2] <- 'cossim_jap'

df_eng <- as_tibble(sub_matrixeng, rownames = 'language')
colnames(df_eng)[2] <- 'cossim_eng'


sub_matrixunr <- cosine_matrix_grammar[ph_lang, unr_lang]
mean_scores_unr <- rowMeans(sub_matrixunr)
mean_scores_unr_matrix <- as.matrix(mean_scores_unr)

df_unr <- as_tibble(mean_scores_unr_matrix, rownames = 'language')
colnames(df_unr)[2] <- 'cossim_unr'

colnames(mean_scores_unr_matrix) <- "Unrelated"


melted_matrix <- melt(cosine_matrix_grammar)

# heatmap
ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are squar



combined_scores <- data.frame(
  Spanish = sub_matrixspan,
  Japanese = sub_matrixjap,
  English = sub_matrixeng,
  Unr = mean_scores_unr_matrix
) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Language",
    values_to = "Similarity_Score"
  )

combined_scores_summary <- combined_scores %>%
  group_by(Language) %>%
  summarize(mean_score = mean(Similarity_Score))


cossim_grammar_density_ridge <- ggplot(combined_scores, aes(x = Similarity_Score, y = Language, fill = Language)) +
  
  geom_density_ridges(alpha = 0.5, scale = 1.2, color = "black") + 
  
  geom_segment(
    data = combined_scores_summary,
    aes(
      x = mean_score, 
      xend = mean_score, 
      y = as.numeric(factor(Language)), 
      yend = as.numeric(factor(Language)) + 0.9, 
      color = Language
    ),
    linetype = "dashed",
    size = 1.2,
    inherit.aes = FALSE # Prevents conflicting with the main plot's y aesthetic mapping
  ) +
  labs(
    title = "Grammar Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Language"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.05)) +
  scale_y_discrete(expand = c(0.01,0)) +
  theme(legend.position = "none")

#---
cossim_grammar_density_ridge
#---

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_ridgeplot.png"),
  plot = cossim_grammar_density_ridge,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)


grammar_cos_s <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','Spanish')), aes(x = Similarity_Score, fill = Language)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = combined_scores_summary %>% filter(Language %in% c('Unrelated','Spanish')),
    aes(xintercept = mean_score, color = Language),
    linetype = "dashed",
    size = 1.2
  ) +
  labs(
    title = "Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Density"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.05))

grammar_cos_e <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','English')), aes(x = Similarity_Score, fill = Language)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = combined_scores_summary %>% filter(Language %in% c('Unrelated','English')),
    aes(xintercept = mean_score, color = Language),
    linetype = "dashed",
    size = 1.2
  ) +
  labs(
    title = "Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Density"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.05))

grammar_cos_j <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','Japanese')), aes(x = Similarity_Score, fill = Language)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = combined_scores_summary %>% filter(Language %in% c('Unrelated','Japanese')),
    aes(xintercept = mean_score, color = Language),
    linetype = "dashed",
    size = 1.2
  ) +
  labs(
    title = "Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Density"
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.05))

grammar_cos_s + grammar_cos_e + grammar_cos_j

# ----- grammar_cossim ver 1 -----------------

GRAMMAR_cossim <- df_span %>% 
  left_join(df_jap,by = 'language') %>% 
  left_join(df_eng,by = 'language') %>% 
  left_join(df_unr,by = 'language') %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

write.csv(GRAMMAR_cossim, file = here("data", "GRAMMAR_cossim.csv"), row.names = TRUE)


# ---- weighted geo_distance from capital using waypoints -------

library(dplyr)
library(geosphere)
library(purrr)
library(ggplot2)
library(sf)
library(maps)
library(sfheaders)

GRAMMAR_cossim <- read.csv(here("data", "GRAMMAR_cossim.csv"))
df <- GRAMMAR_cossim
# load MGT route tree
nodes <- read.csv(here("data", "nodes.csv"))
edges <- read.csv(here("data", "edges.csv"))

nodes$id <- as.character(nodes$id)
edges$from <- as.character(edges$from)
edges$to <- as.character(edges$to)
reverse_edges <- edges %>% rename(from = to, to = from)
edges <- bind_rows(edges, reverse_edges) %>% distinct()

ref_coords1 <- c(121,14.6)

compute_shortest_path_df <- function(df, ref_coords1, nodes, edges, land_penalty = 4.44) {
  
  land_penalty <- 4.44
  df <- GRAMMAR_cossim
  # Ensure IDs are character
  nodes <- nodes %>% mutate(id = as.character(id))
  edges <- edges %>% mutate(from = as.character(from), to = as.character(to))
  
  # Add reverse edges for bidirectional routing
  reverse_edges <- edges %>% rename(from = to, to = from)
  edges <- bind_rows(edges, reverse_edges) %>% distinct()
  
  # Load land polygons (Philippines + Malaysia)
  world_map <- map_data("world") %>% filter(region %in% c("Philippines", "Malaysia"))
  
  land_sf <- sf_polygon(
    obj = world_map,
    polygon_id = "group",
    x = "long",
    y = "lat"
  ) %>%
    st_union() %>%
    st_sf(geometry = .)
  
  # Compute edge weights and terrain-aware cost
  edges <- edges %>%
    rowwise() %>%
    mutate(weight = {
      from_coords <- nodes %>% filter(id == from)
      to_coords <- nodes %>% filter(id == to)
      if (nrow(from_coords) == 0 || nrow(to_coords) == 0) NA_real_
      else distHaversine(c(from_coords$longitude, from_coords$latitude),
                         c(to_coords$longitude, to_coords$latitude))
    }) %>%
    ungroup()
  
  land_sf <- land_sf %>% 
    st_set_crs(4326)
  
  edge_lines <- edges %>%
    rowwise() %>%
    mutate(
      geometry = list(st_linestring(matrix(c(
        nodes$longitude[nodes$id == from],
        nodes$latitude[nodes$id == from],
        nodes$longitude[nodes$id == to],
        nodes$latitude[nodes$id == to]
      ), ncol = 2, byrow = TRUE)))
    ) %>%
    ungroup() %>%
    st_as_sf(crs = 4326)
  
  
  
  edge_lines <- edge_lines %>%
    rowwise() %>%
    mutate(
      land_part = list(st_intersection(geometry, land_sf)),
      sea_part  = list(st_difference(geometry, land_sf)),
      
      land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
      sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),
      
      weighted_cost = land_len * land_penalty + sea_len,
      crosses_land  = land_len > 0
      
    ) %>%
    ungroup()
  
  # Build complete graph with weighted costs
  all_ids <- unique(c(edge_lines$from, edge_lines$to))
  graph <- lapply(all_ids, function(id) {
    neighbors <- edge_lines %>% filter(from == id) %>% select(to, weighted_cost)
    if (nrow(neighbors) == 0) tibble(to = character(), weight = numeric())
    else rename(neighbors, weight = weighted_cost)
  })
  names(graph) <- all_ids
  
  # Helper: find nearest node to a coordinate
  find_nearest_node <- function(coords) {
    distances <- distHaversine(matrix(c(nodes$longitude, nodes$latitude), ncol = 2),
                               coords)
    nodes$id[which.min(distances)]
  }
  
  # Shortest path function (returns distance and trace)
  shortest_path_trace <- function(start_id, end_id, graph) {
    visited <- setNames(rep(FALSE, length(graph)), names(graph))
    dist <- setNames(rep(Inf, length(graph)), names(graph))
    prev <- setNames(rep(NA_character_, length(graph)), names(graph))
    dist[start_id] <- 0
    queue <- data.frame(id = start_id, dist = 0)
    
    while (nrow(queue) > 0) {
      queue <- queue[order(queue$dist), ]
      current <- queue$id[1]
      current_dist <- queue$dist[1]
      queue <- queue[-1, ]
      
      if (visited[current]) next
      visited[current] <- TRUE
      
      neighbors <- graph[[current]]
      if (is.null(neighbors)) next
      
      for (i in seq_len(nrow(neighbors))) {
        neighbor <- neighbors$to[i]
        weight <- neighbors$weight[i]
        if (is.na(weight)) next
        if (dist[neighbor] > current_dist + weight) {
          dist[neighbor] <- current_dist + weight
          prev[neighbor] <- current
          queue <- rbind(queue, data.frame(id = neighbor, dist = dist[neighbor]))
        }
      }
    }
    
    if (!is.finite(dist[end_id])) return(list(distance = NA_real_, path = NULL))
    
    # Reconstruct path
    path <- end_id
    while (!is.na(prev[path[1]])) {
      path <- c(prev[path[1]], path)
    }
    return(list(distance = dist[end_id], path = path))
  }
  
  # Plot function with land overlay
  plot_path <- function(path_ids, nodes,
                        start_coords = NULL,
                        end_coords = NULL,
                        land_part_start = NULL,
                        sea_part_start  = NULL,
                        land_part_end   = NULL,
                        sea_part_end    = NULL) {
    
    path_df <- nodes %>%
      filter(id %in% path_ids) %>%
      arrange(factor(id, levels = path_ids))
    
    # Build connector segment sf objects with land/sea labels
    connector_sf <- list()
    
    if (!is.null(land_part_start)) {
      connector_sf <- append(connector_sf, list(
        st_sf(geometry = st_sfc(land_part_start, crs = 4326), crosses_land = TRUE)
      ))
    }
    if (!is.null(sea_part_start)) {
      connector_sf <- append(connector_sf, list(
        st_sf(geometry = st_sfc(sea_part_start, crs = 4326), crosses_land = FALSE)
      ))
    }
    if (!is.null(land_part_end)) {
      connector_sf <- append(connector_sf, list(
        st_sf(geometry = st_sfc(land_part_end, crs = 4326), crosses_land = TRUE)
      ))
    }
    if (!is.null(sea_part_end)) {
      connector_sf <- append(connector_sf, list(
        st_sf(geometry = st_sfc(sea_part_end, crs = 4326), crosses_land = FALSE)
      ))
    }
    
    connector_segments <- do.call(rbind, connector_sf)
    
    # Build plot
    ggplot() +
      geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
                   fill = "gray95", color = "gray70") +
      geom_sf(data = connector_segments, aes(color = crosses_land), size = 1.2) +
      geom_path(data = path_df, aes(x = longitude, y = latitude),
                color = "black", size = 1.2) +
      geom_point(data = path_df, aes(x = longitude, y = latitude),
                 color = "black", size = 2) +
      { if (!is.null(start_coords)) geom_point(aes(x = start_coords[1], y = start_coords[2]),
                                               color = "red", size = 3) } +
      { if (!is.null(end_coords)) geom_point(aes(x = end_coords[1], y = end_coords[2]),
                                             color = "green", size = 3) } +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
      coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
      theme_minimal() +
      labs(title = paste("Path:", paste(path_ids, collapse = " → ")),
           color = "Crosses Land")
  }
  
  
  # Precompute destination node
  ref_nearest <- find_nearest_node(ref_coords1)
  
  df <- df %>%
    rowwise() %>%
    mutate(
      start_coords = list(c(longitude, latitude)),
      start_nearest = find_nearest_node(start_coords),
      
      connector_start_geom = list(st_linestring(rbind(
        start_coords,
        c(nodes$longitude[nodes$id == start_nearest],
          nodes$latitude[nodes$id == start_nearest])
      ))),
      
      connector_end_geom = list(st_linestring(rbind(
        ref_coords1,
        c(nodes$longitude[nodes$id == ref_nearest],
          nodes$latitude[nodes$id == ref_nearest])
      )))
    ) %>%
    ungroup()
  
  
  
  connector_lines <- df %>%
    select(connector_start_geom, connector_end_geom) %>%
    mutate(
      connector_start_geom = st_sfc(connector_start_geom, crs = 4326),
      connector_end_geom   = st_sfc(connector_end_geom, crs = 4326)
    )
  
  
  connector_lines <- connector_lines %>%
    rowwise() %>%
    mutate(
      land_part_start = list(st_intersection(connector_start_geom, land_sf)),
      sea_part_start  = list(st_difference(connector_start_geom, land_sf)),
      
      land_len_start = as.numeric(if (!is.null(land_part_start) && length(land_part_start) > 0) st_length(land_part_start) else 0),
      sea_len_start  = as.numeric(if (!is.null(sea_part_start)  && length(sea_part_start)  > 0) st_length(sea_part_start)  else 0),
      
      land_part_end = list(st_intersection(connector_end_geom, land_sf)),
      sea_part_end  = list(st_difference(connector_end_geom, land_sf)),
      
      land_len_end = as.numeric(if (!is.null(land_part_end) && length(land_part_end) > 0) st_length(land_part_end) else 0),
      sea_len_end  = as.numeric(if (!is.null(sea_part_end)  && length(sea_part_end)  > 0) st_length(sea_part_end)  else 0),
      
      connector_start_penalty = land_len_start * land_penalty + sea_len_start,
      connector_end_penalty   = land_len_end   * land_penalty + sea_len_end
    ) %>%
    ungroup()
  
  
  df <- df %>%
    mutate(row_id = row_number()) %>%
    left_join(connector_lines %>% mutate(row_id = row_number()), by = "row_id")
  
  extract_sfg <- function(segment_list) {
    lapply(segment_list, function(x) {
      if (!is.null(x) && length(x) > 0 && inherits(x[[1]], "sfg")) {
        x[[1]]
      } else {
        NULL
      }
    })
  }
  
  
  df$land_geom_start <- extract_sfg(df$land_part_start)
  
  
  df <- df %>%
    mutate(
      land_geom_start = extract_sfg(land_part_start),
      sea_geom_start  = extract_sfg(sea_part_start),
      land_geom_end   = extract_sfg(land_part_end),
      sea_geom_end    = extract_sfg(sea_part_end)
    )
  
  df <- df %>%
    rowwise() %>%
    mutate(
      trace_result = list(shortest_path_trace(start_nearest, ref_nearest, graph)),
      tree_dist = trace_result$distance,
      path_nodes = list(trace_result$path),
      
      geodist_H1_span = if (is.na(tree_dist)) NA_real_ else
        (connector_start_penalty + tree_dist + connector_end_penalty) / 1000,
      plot = list(plot_path(
        path_ids = trace_result$path,
        nodes = nodes,
        start_coords = start_coords,
        end_coords = ref_coords1,
        land_part_start = land_geom_start,
        sea_part_start  = sea_geom_start,
        land_part_end   = land_geom_end,
        sea_part_end    = sea_geom_end
      ))
    ) %>%
    ungroup()
  
  
  return(df)
}

GRAMMAR_cossim <- compute_shortest_path_df(GRAMMAR_cossim, ref_coords1, nodes, edges, land_penalty = 44.18)

df <- GRAMMAR_cossim 

# ------- ROUTE PLOTTING ----------------
world_map <- map_data("world") %>%
  filter(region %in% c("Philippines", "Malaysia"))

ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_sf(data = edge_lines, aes(color = crosses_land), size = 1.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal() +
  labs(title = "Edge Segments: Land Crossing Detection",
       color = "Crosses Land")

# Extract valid land segments
land_segments <- edge_lines$land_part[
  sapply(edge_lines$land_part, function(x) {
    !is.null(x) && length(x) > 0 && inherits(x[[1]], "sfg")
  })
]

# Flatten to sfg list
land_geoms <- lapply(land_segments, `[[`, 1)

# Build sf object
land_segments_sf <- st_sf(geometry = st_sfc(land_geoms, crs = 4326))


# Extract valid sea segments
sea_segments <- edge_lines$sea_part[
  sapply(edge_lines$sea_part, function(x) {
    !is.null(x) && length(x) > 0 && inherits(x[[1]], "sfg")
  })
]

sea_geoms <- lapply(sea_segments, `[[`, 1)

sea_segments_sf <- st_sf(geometry = st_sfc(sea_geoms, crs = 4326))

# ----- plot paths -------------------------

print(GRAMMAR_cossim$plot[[12]])

land_segments_sf$crosses_land <- TRUE
sea_segments_sf$crosses_land  <- FALSE
main_path_sf <- rbind(land_segments_sf, sea_segments_sf)

# Build arrow_sf and connector_sf from a single loop
arrow_geom_list <- list()
arrow_land_flag <- logical()

for (i in seq_len(nrow(df))) {
  if (!is.null(df$land_geom_start[[i]])) {
    arrow_geom_list <- append(arrow_geom_list, list(df$land_geom_start[[i]]))
    arrow_land_flag <- append(arrow_land_flag, TRUE)
  }
  if (!is.null(df$sea_geom_start[[i]])) {
    arrow_geom_list <- append(arrow_geom_list, list(df$sea_geom_start[[i]]))
    arrow_land_flag <- append(arrow_land_flag, FALSE)
  }
  if (!is.null(df$land_geom_end[[i]])) {
    arrow_geom_list <- append(arrow_geom_list, list(df$land_geom_end[[i]]))
    arrow_land_flag <- append(arrow_land_flag, TRUE)
  }
  if (!is.null(df$sea_geom_end[[i]])) {
    arrow_geom_list <- append(arrow_geom_list, list(df$sea_geom_end[[i]]))
    arrow_land_flag <- append(arrow_land_flag, FALSE)
  }
}

arrow_sf <- st_sf(
  crosses_land = arrow_land_flag,
  geometry     = st_sfc(arrow_geom_list, crs = 4326)
)

connector_sf <- arrow_sf

# Tag and combine
main_path_sf$source <- "main"
connector_sf$source <- "connector"
main_path_sf        <- st_set_crs(main_path_sf, 4326)
full_tree_sf        <- rbind(main_path_sf, connector_sf)

# Diagnostic plot
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_sf(data = full_tree_sf, size = 1.5, color = "black") +
  geom_sf(data = arrow_sf,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "black") +
  geom_point(data = GRAMMAR_cossim, aes(x = longitude, y = latitude),
             size = 3, shape = 21) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal() +
  labs(title = "Grammatic Historical Routes Waypoint System",
       x = "Longitude", y = "Latitude")

# Arrow overlay plot (saved)
arrow_plot_GRAMMAR <- ggplot() +
  geom_sf(data = full_tree_sf,
          mapping = aes(geometry = geometry),
          color = "black", linewidth = 0.7) +
  geom_sf(data = arrow_sf,
          mapping = aes(geometry = geometry),
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "black", linewidth = 0.7) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal()

arrow_plot_GRAMMAR
saveRDS(arrow_plot_GRAMMAR, file = here("data", "grammar_waypoint_plot.rds"))

# ----- linear model -------------------------------------------

lm_span <- lm(cossim_span ~ geodist_H1_span, data = GRAMMAR_cossim)
summary_lm_span <- summary(lm_span)

slope <- round(coef(lm_span)[2],5)
coefficient <- round(coef(lm_span)[1],5)
r_squared <- round(summary_lm_span$r.squared,3)
linear_model_equation <- paste0("Y = ", slope, "x + ", coefficient)

ggplot(data = GRAMMAR_cossim, aes(x = geodist_H1_span, y = cossim_span)) +
  geom_point() +
  geom_smooth(method = 'lm',se = FALSE) +
  theme_bw() +
  labs(title = 'Linear Model',
         x = 'Relative Migration Distance (km)',
         y = 'Cosine Similarity')

# ---- exponential case -------------------------------------

library(rethinking)

# Standardize predictors
GRAMMAR_cossim$dist_std <- standardize(GRAMMAR_cossim$geodist_H1_span)

# Model: exponential decay
m_exp <- rethinking::map(
  alist(
    cossim_span ~ dnorm(mu, sigma),
    mu <- a * exp(-b * dist_std),
    a ~ dnorm(0.5, 0.2),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = GRAMMAR_cossim
)

# Summarize
precis(m_exp)

# Generate predictions
dist_seq <- seq(from = min(GRAMMAR_cossim$dist_std), to = max(GRAMMAR_cossim$dist_std), length.out = 100)
preds <- link(m_exp, data = data.frame(dist_std = dist_seq))
mu_mean <- apply(preds, 2, mean)
mu_PI <- apply(preds, 2, PI)
mu_PI_t <- t(mu_PI)
# Plot
scatter_df <- GRAMMAR_cossim

# Line + ribbon data
ribbon_df <- data.frame(
  dist = dist_seq,
  mean = mu_mean,
  lower = mu_PI_t[,1],
  upper = mu_PI_t[,2]
)



ggplot() +
  geom_point(data = scatter_df, aes(x = dist_std, y = cossim_span),
             color = "black", size = 2) +
  geom_line(data = ribbon_df, aes(x = dist, y = mean),
            color = "blue", linewidth = 1.2) +
  theme_minimal() +
  labs(
    x = "Relative Migration Distance (km)",
    y = "Cosine Similarity",
    title = "Exponential Model"
  ) + 
  theme_bw()


# Extract coefficients
coefs <- coef(m_exp)
a <- round(coefs["a"], 5)
b <- round(coefs["b"], 5)

# Write out the equation
exp_eq <- paste0("μ = ", a, " * exp(-", b, " * x)")
print(exp_eq)

# ------- loess regression ---------------



# ----- spline regression

library(splines)

spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = GRAMMAR_cossim)
x_vals <- GRAMMAR_cossim$geodist_H1_span
y_spline <- predict(spline_fit, newdata = data.frame(geodist_H1_span = x_vals))

ggplot(GRAMMAR_cossim, aes(x = geodist_H1_span, y = cossim_span)) +
  geom_point(color = "gray") +
  geom_line(aes(y = y_spline), color = "blue", linewidth = 1.2) +
  theme_bw() +
  labs(
    title = "Spline Regression",
    x = "Relative Migration Distance (km)",
    y = "Cosine Similarity"
  ) 

paste0("y(x) = ", paste0("β", 1:5, "·B", 1:5, "(x)", collapse = " + "))
coefs <- coef(spline_fit)

paste0(
  "y(x) = ", round(coefs[1], 5), " + ",
  paste0(round(coefs[-1], 5), "·B", 1:5, "(x)", collapse = " + ")
)

# 1. Extract the spline basis object from the model frame
spline_basis_object <- model.frame(spline_fit)$'bs(geodist_H1_span, df = 5)'

# 2. Extract the 'knots' attribute from that object
knots <- attr(spline_basis_object, "knots")

# Print the resulting knot vector
print(knots)

# ----- comparing errors ---------------------------


spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = GRAMMAR_cossim)

# Linear model predictions
GRAMMAR_cossim$pred_linear <- predict(lm_span)

# Exponential model predictions using link()
exp_preds <- link(m_exp)
GRAMMAR_cossim$pred_exp <- apply(exp_preds, 2, mean)

# RMSE calculation
rmse_linear <- sqrt(mean((GRAMMAR_cossim$cossim_span - GRAMMAR_cossim$pred_linear)^2))
rmse_exp <- sqrt(mean((GRAMMAR_cossim$cossim_span - GRAMMAR_cossim$pred_exp)^2))
rmse_spline <- sqrt(mean((GRAMMAR_cossim$cossim_span - y_spline)^2))

p_linear <- summary(lm_span)$coefficients[2, 4]  # Slope p-value
p_spline <- pf(summary(spline_fit)$fstatistic[1],
               summary(spline_fit)$fstatistic[2],
               summary(spline_fit)$fstatistic[3],
               lower.tail = FALSE)

exp_post <- extract.samples(m_exp)
pr_b_gt_0 <- mean(exp_post$b > 0)  # Replace 'b' with your actual slope parameter name

# Compare

cossim_span_range <- max(GRAMMAR_cossim$cossim_span) - min(GRAMMAR_cossim$cossim_span)

rmse_df <- data.frame(
  Model = c("Linear", "Exponential", "Spline"),
  RMSE = round(c(rmse_linear, rmse_exp, rmse_spline), 5),
  normalized_rmse = round(c(rmse_linear, rmse_exp, rmse_spline) / cossim_span_range, 5),
  p_value = c(round(p_linear, 5), NA, round(p_spline, 5)),
  posterior_prob = c(NA, round(pr_b_gt_0, 3), NA)
)

print(rmse_df)

# ---- dissimilarity matrix ----------------------------------

ph_lang <- GRAMBANKdf_PH %>% 
  filter(Language_Type == 'Philippine Language') %>% 
  pull(language)

cosine_matrix_phil_grammar <- cosine_matrix_grammar[ph_lang, ph_lang]
lang_order <- rownames(cosine_matrix_phil_grammar)

cosine_matrix_phil_grammar <- 1-cosine_matrix_phil_grammar
diag(cosine_matrix_phil_grammar) <- 0

GRAMMAR_diss_matrix <- cosine_matrix_phil_grammar

# ----- distance matrix --------------------------------------

phil_df <- GRAMBANKdf_PH %>%
  filter(Language_Type == "Philippine Language") %>%
  mutate(
    start_coords = map2(Longitude, Latitude, ~ c(.x, .y)),
    nearest_node = map_chr(start_coords, find_nearest_node)
  )

land_penalty <- 4.44

connector_df <- phil_df %>%
  mutate(
    connector_geom = map2(start_coords, nearest_node, ~ st_linestring(rbind(
      .x,
      c(nodes$longitude[nodes$id == .y], nodes$latitude[nodes$id == .y])
    )))
  )

connector_df <- connector_df %>%
  mutate(connector_geom_sfc = st_sfc(connector_geom, crs = 4326))

connector_df <- connector_df %>%
  rowwise() %>%
  mutate(
    land_part = list(st_intersection(connector_geom_sfc, land_sf)),
    sea_part  = list(st_difference(connector_geom_sfc, land_sf)),
    
    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),
    
    connector_penalty = land_len * land_penalty + sea_len
  ) %>%
  ungroup()


phil_pairs <- expand_grid(lang1 = phil_df$language, lang2 = phil_df$language) %>%
  filter(lang1 != lang2) %>%
  left_join(phil_df %>% select(language, node1 = nearest_node), by = c("lang1" = "language")) %>%
  left_join(phil_df %>% select(language, node2 = nearest_node), by = c("lang2" = "language")) %>%
  left_join(connector_df %>% select(language, penalty1 = connector_penalty), by = c("lang1" = "language")) %>%
  left_join(connector_df %>% select(language, penalty2 = connector_penalty), by = c("lang2" = "language"))

phil_pairs <- phil_pairs %>%
  rowwise() %>%
  mutate(
    trace = list(shortest_path_trace(node1, node2, graph)),
    tree_dist = trace$distance,
    geodist_H1_span = if (is.na(tree_dist)) NA_real_ else
      (penalty1 + tree_dist + penalty2) / 1000
  ) %>%
  ungroup()

dist_matrix <- phil_pairs %>%
  select(lang1, lang2, geodist_H1_span) %>%
  pivot_wider(names_from = lang2, values_from = geodist_H1_span) %>%
  column_to_rownames("lang1") %>%
  as.matrix()

GRAMMAR_dist_matrix <- dist_matrix[ph_lang, ph_lang]
GRAMMAR_dist_matrix[is.na(GRAMMAR_dist_matrix)] <- 0

# ---- plot matrices -----------------------------------

melt_grammar_dist_matrix <- melt(GRAMMAR_dist_matrix)
melt_grammar_diss_matrix <- melt(GRAMMAR_diss_matrix)

dist_matrix_p <- ggplot(melt_grammar_dist_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are square

diss_matrix_p  <- ggplot(melt_grammar_diss_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are square

dist_matrix_p + diss_matrix_p 


# ---- mantel test ----------------------------------------
library(vegan)

# Convert to distance objects
x_dist <- as.dist(GRAMMAR_dist_matrix)
y_dist <- as.dist(GRAMMAR_diss_matrix)

# Run Mantel test
mantel_result <- mantel(x_dist, y_dist, method = "spearman", permutations = 999)

print(mantel_result)

library(ggplot2)

# Convert matrices to vectors
x_vec <- as.vector(as.dist(GRAMMAR_dist_matrix))
y_vec <- as.vector(as.dist(GRAMMAR_diss_matrix))




# Plot
ggplot(data.frame(x = x_vec, y = y_vec), aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  theme_bw() +
  labs(
    title = "Mantel Test: Grammar Feature Distance vs. Dissimilarity",
    x = "Relative Migration Pairwise Distance (km)",
    y = "Grammatic Dissimilarity"
  )

# ----- saving dataframes -----------------------------------

write.csv(cosine_matrix_grammar, file = here("data", "GRAMMAR_cosine_matrix.csv"), row.names = TRUE)


GRAMMAR_cossim <- df_span %>% 
  left_join(df_jap,by = 'language') %>% 
  left_join(df_eng,by = 'language') %>% 
  left_join(df_unr,by = 'language') %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

write.csv(GRAMMAR_cossim, file = here("data", "GRAMMAR_cossim.csv"), row.names = TRUE)




