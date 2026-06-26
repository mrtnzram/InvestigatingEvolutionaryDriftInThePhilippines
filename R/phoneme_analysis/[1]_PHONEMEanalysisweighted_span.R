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
library(patchwork)
library(igraph)
library(sfheaders)
library(here)
library(ggridges)


# ---- Prepping Globals --------------------------------------------------------------------------

PHOIBLEdf <- read_csv(here("data", "PHOIBLEdf_PH.csv"))

PHOIBLEdf <- PHOIBLEdf %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

ph_lang <- c('Tiruray', 'Maranao', 'Central Sama', 'Hiligaynon', 'Aklanon', 'Filipino', 'Coastal-Naga Bikol','Tagalog','Pampanga','Pangasinan','Iloko','Dupaninan Agta',
             'Cebuano','Southern Sama')
int_lang <- c('English','Spanish','Japanese')
unr_lang <- c("Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Liv",
              "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")


phoneme_cols <- PHOIBLEdf %>% 
  select(-language,-source,-iso6393,-Language_type,-latitude,-longitude)

phoneme_cols <- colnames(phoneme_cols)

PHOIBLEdf_full<- phoible.feature(source='all',na.rm = FALSE)


phoneme_freq <- read_csv(here("data", "phoneme_freq.csv"))



# ------ weighted similarity ---------------------


calculate_language_similarity <- function(PHOIBLEdf, phoneme_freq, phoneme_cols, id_col = "language") {
  
  # Ensure data frames are aligned by phoneme order
  # This is a crucial step to ensure weights match the correct phonemes
  phoneme_freq <- phoneme_freq %>%
    filter(phoneme %in% phoneme_cols) %>%
    arrange(match(phoneme, phoneme_cols))
  
  # Extract binary phoneme data
  binary_data <- PHOIBLEdf %>%
    select(all_of(phoneme_cols))
  
  # Extract language IDs
  language_ids <- PHOIBLEdf[[id_col]]
  
  # Calculate the weights for 1-1 and 0-0 matches
  w_11 <- phoneme_freq$IDF
  
  # Frequencies are needed to calculate w_00.
  # We can get frequencies from IDF: freq = exp(-IDF)
  frequencies <- exp(-w_11)
  w_00 <- -log2(1 - frequencies)
  
  # Calculate the denominator (total possible score) once
  denominator <- sum(w_11 + w_00)
  
  # Get the number of languages
  n_languages <- nrow(binary_data)
  
  # Initialize an empty matrix to store results
  similarity_matrix <- matrix(0, nrow = n_languages, ncol = n_languages,
                              dimnames = list(language_ids, language_ids))
  
  # Loop through all unique pairs of languages to compute scores
  for (i in 1:n_languages) {
    for (j in i:n_languages) {
      
      # Extract the two vectors to compare
      vec_a <- as.numeric(binary_data[i, ])
      vec_b <- as.numeric(binary_data[j, ])
      
      # Calculate the numerator for this pair
      numerator <- sum(
        # 1-1 matches weighted by w_11
        (vec_a == 1 & vec_b == 1) * w_11,
        # 0-0 matches weighted by w_00
        (vec_a == 0 & vec_b == 0) * w_00
      )
      
      # Calculate the normalized score and fill the matrix
      score <- numerator / denominator
      similarity_matrix[i, j] <- score
      similarity_matrix[j, i] <- score # Matrix is symmetric
    }
  }
  
  return(similarity_matrix)
}

similarity_matrix <- calculate_language_similarity(
  PHOIBLEdf = PHOIBLEdf,
  phoneme_freq = phoneme_freq,
  phoneme_cols = phoneme_cols,
  id_col = "language"
)

# Investigate --------------------- 

similarity_matrix['Tagalog','Spanish']
similarity_matrix['Tagalog','Japanese']
similarity_matrix['Tagalog','Tagalog']

sub_matrixspan <- similarity_matrix[ph_lang, 'Spanish']
df_span_ss <- as_tibble(sub_matrixspan, rownames = 'language')
colnames(df_span_ss)[2] <- 'sim'

sub_matrixjap <- similarity_matrix[ph_lang, 'Japanese']
df_jap_ss <- as_tibble(sub_matrixjap, rownames = 'language')
colnames(df_jap_ss)[2] <- 'sim'

sub_matrixeng <- similarity_matrix[ph_lang, 'English']
df_eng_ss <- as_tibble(sub_matrixeng, rownames = 'language')
colnames(df_eng_ss)[2] <- 'sim'

sub_matrixunr <- similarity_matrix[ph_lang, unr_lang]
mean_scores_unr <- rowMeans(sub_matrixunr)
mean_scores_unr_matrix <- as.matrix(mean_scores_unr)
colnames(mean_scores_unr_matrix) <- "Unrelated"


melted_matrix <- melt(similarity_matrix)

# heatmap
ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "Language Similarity Heatmap", x = "", y = "") +
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

ggplot(combined_scores, aes(x = Similarity_Score, fill = Language)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = combined_scores_summary,
    aes(xintercept = mean_score, color = Language),
    linetype = "dashed",
    size = 1.2
  ) +
  labs(
    title = "Weighted Similarity Distribution",
    x = "Similarity Score",
    y = "Density"
  ) +
  theme_bw()










# investigate deeper on a phoneme level ------------------------------------------------------------------

investigate_language_pair <- function(PHOIBLEdf, phoneme_freq, phoneme_cols, id_col, lang1_id, lang2_id) {
  
  # Step 1: Extract and align the two language vectors
  lang1_row <- PHOIBLEdf %>%
    dplyr::filter(!!rlang::sym(id_col) == lang1_id) %>%
    dplyr::select(dplyr::all_of(phoneme_cols))
  
  lang2_row <- PHOIBLEdf %>%
    dplyr::filter(!!rlang::sym(id_col) == lang2_id) %>%
    dplyr::select(dplyr::all_of(phoneme_cols))
  
  # Check if languages were found
  if (nrow(lang1_row) == 0 || nrow(lang2_row) == 0) {
    stop("One or both language IDs not found.")
  }
  
  # Step 2: Prepare the phoneme frequency data
  aligned_freq <- phoneme_freq %>%
    dplyr::filter(phoneme %in% phoneme_cols) %>%
    dplyr::arrange(match(phoneme, phoneme_cols))
  
  w_11 <- aligned_freq$IDF
  frequencies <- exp(-w_11)
  w_00 <- -log2(1 - frequencies)
  
  phoneme_info <- tibble(
    phoneme = phoneme_cols,
    w_11_weight = w_11,
    w_00_weight = w_00
  )
  
  # Step 3: Compare the two language vectors
  comparison_df <- phoneme_info %>%
    dplyr::mutate(
      lang1_value = as.numeric(lang1_row[1, ]),
      lang2_value = as.numeric(lang2_row[1, ])
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      match_type = case_when(
        lang1_value == 1 & lang2_value == 1 ~ "1-1 Match",
        lang1_value == 0 & lang2_value == 0 ~ "0-0 Match",
        TRUE ~ "Mismatch"
      ),
      match_weight = case_when(
        match_type == "1-1 Match" ~ w_11_weight,
        match_type == "0-0 Match" ~ w_00_weight,
        TRUE ~ 0
      ),
      # New column to show the lost weight for mismatches
      mismatch_lost_weight = case_when(
        match_type == "Mismatch" & lang1_value == 1 ~ w_11_weight,
        match_type == "Mismatch" & lang2_value == 1 ~ w_11_weight,
        TRUE ~ 0
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(phoneme, lang1_value, lang2_value, match_type, match_weight, mismatch_lost_weight)
  
  return(comparison_df)
}

comparison_results <- investigate_language_pair(
  PHOIBLEdf = PHOIBLEdf,
  phoneme_freq = phoneme_freq,
  phoneme_cols = phoneme_cols,
  id_col = "language",
  lang1_id = "Tagalog",
  lang2_id = "Spanish"
)

comparison_results2 <- investigate_language_pair(
  PHOIBLEdf = PHOIBLEdf,
  phoneme_freq = phoneme_freq,
  phoneme_cols = phoneme_cols,
  id_col = "language",
  lang1_id = "Tagalog",
  lang2_id = "Japanese"
)

# view the detailed results

sum(comparison_results$match_weight)
sum(comparison_results2$match_weight)

sum(comparison_results$mismatch_lost_weight)
sum(comparison_results2$mismatch_lost_weight)


# analyze the distribution of match types
comparison_results %>%
  dplyr::group_by(match_type) %>%
  dplyr::summarise(
    count = n(),
    total_weight = sum(match_weight)
  )

comparison_results2 %>%
  dplyr::group_by(match_type) %>%
  dplyr::summarise(
    count = n(),
    total_weight = sum(match_weight)
  )

















# ----- cosine similarity ------------------------

calculate_weighted_cosine_similarity <- function(PHOIBLEdf, phoneme_freq, phoneme_cols, id_col = "language") {
  
  # Extract and align the binary data and IDF weights
  # Ensure the phoneme frequencies are in the same order as the phoneme columns
  aligned_freq <- phoneme_freq %>%
    dplyr::filter(phoneme %in% phoneme_cols) %>%
    dplyr::arrange(match(phoneme, phoneme_cols))
  
  idf_weights <- aligned_freq$IDF
  
  # Extract the binary phoneme data
  binary_data <- PHOIBLEdf %>%
    dplyr::select(dplyr::all_of(phoneme_cols)) %>%
    as.matrix() # Convert to a matrix for faster calculations
  
  # Extract language IDs for matrix naming
  language_ids <- PHOIBLEdf[[id_col]]
  
  # Step 2: Create a weighted phoneme matrix
  # Multiply each column of the binary matrix by its corresponding IDF weight
  weighted_data <- sweep(binary_data, 2, idf_weights, FUN = "*")
  
  # Step 3: Calculate the cosine similarity matrix
  n_languages <- nrow(weighted_data)
  cosine_matrix <- matrix(0, nrow = n_languages, ncol = n_languages,
                          dimnames = list(language_ids, language_ids))
  
  # A small epsilon to avoid division by zero for languages with no phonemes
  epsilon <- 1e-9
  
  # Loop through all unique pairs of languages
  for (i in 1:n_languages) {
    for (j in i:n_languages) {
      
      vec_a <- weighted_data[i, ]
      vec_b <- weighted_data[j, ]
      
      # Cosine Similarity Formula: (A . B) / (||A|| * ||B||)
      # Numerator is the dot product
      dot_product <- sum(vec_a * vec_b)
      
      # Denominator is the product of the magnitudes (Euclidean norms)
      magnitude_a <- sqrt(sum(vec_a^2))
      magnitude_b <- sqrt(sum(vec_b^2))
      
      denominator <- (magnitude_a * magnitude_b) + epsilon
      
      score <- dot_product / denominator
      
      cosine_matrix[i, j] <- score
      cosine_matrix[j, i] <- score # Matrix is symmetric
    }
  }
  
  return(cosine_matrix)
}



cosine_matrix <- calculate_weighted_cosine_similarity(
  PHOIBLEdf, 
  phoneme_freq, 
  phoneme_cols, 
  id_col = "language")

# INVESTIGATE --------------------------------------------------

ordered_languages <- PHOIBLEdf %>%
  arrange(Language_type) %>%
  pull(language)


cosine_matrix <- cosine_matrix[ordered_languages, ordered_languages]


cosine_matrix['Tagalog','Spanish'] 
cosine_matrix['Tagalog','Japanese'] 

sub_matrixspan <- cosine_matrix[ph_lang, 'Spanish']
sub_matrixjap <- cosine_matrix[ph_lang, 'Japanese']
sub_matrixeng <- cosine_matrix[ph_lang, 'English']


df_span <- as_tibble(sub_matrixspan, rownames = 'language')
colnames(df_span)[2] <- 'cossim_span'

df_jap <- as_tibble(sub_matrixjap, rownames = 'language')
colnames(df_jap)[2] <- 'cossim_jap'

df_eng <- as_tibble(sub_matrixeng, rownames = 'language')
colnames(df_eng)[2] <- 'cossim_eng'


sub_matrixunr <- cosine_matrix[ph_lang, unr_lang]
mean_scores_unr <- rowMeans(sub_matrixunr)
mean_scores_unr_matrix <- as.matrix(mean_scores_unr)

df_unr <- as_tibble(mean_scores_unr_matrix, rownames = 'language')
colnames(df_unr)[2] <- 'cossim_unr'

colnames(mean_scores_unr_matrix) <- "Unrelated"


melted_matrix <- melt(cosine_matrix)

# heatmap
ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are square



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

cossim_phoneme_density_ridge <- ggplot(combined_scores, aes(x = Similarity_Score, y = Language, fill = Language)) +
  # Changed to geom_density_ridges
  geom_density_ridges(alpha = 0.5, scale = 1.2, color = "black") + 
  
  # Replaced geom_vline with geom_segment to keep lines contained within each ridge
  geom_segment(
    data = combined_scores_summary,
    aes(
      x = mean_score, 
      xend = mean_score, 
      y = as.numeric(factor(Language)), 
      yend = as.numeric(factor(Language)) + 0.9, # Adjust height of the line
      color = Language
    ),
    linetype = "dashed",
    size = 1.2,
    inherit.aes = FALSE # Prevents conflicting with the main plot's y aesthetic mapping
  ) +
  labs(
    title = "Phoneme Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Language"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.05)) +
  scale_y_discrete(expand = c(0.01,0)) +
  theme(legend.position = "none")

#---
cossim_phoneme_density_ridge
#---

ggsave(
  filename = here("figures", "phoneme", "distributions", "phoneme_ridgeplot.png"),
  plot = cossim_phoneme_density_ridge,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)


# Individual Plots ---- 
phoneme_cos_s <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','Spanish')), aes(x = Similarity_Score, fill = Language)) +
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
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.02))

phoneme_cos_e <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','English')), aes(x = Similarity_Score, fill = Language)) +
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
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.02))

phoneme_cos_j <- ggplot(combined_scores %>% filter(Language %in% c('Unrelated','Japanese')), aes(x = Similarity_Score, fill = Language)) +
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
  scale_x_continuous(breaks = seq(0, 0.4, by = 0.02))

phoneme_cos_s + phoneme_cos_e + phoneme_cos_j


# ---- non-independent non-parametric statistical tests for similarity distributions --------------




# ---- distance metric -------------------

df_span<- df_span %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))
df_jap<- df_jap %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))
df_eng<- df_eng %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))


# ---- weighted geo_distance from capital using waypoints -------

library(dplyr)
library(geosphere)
library(purrr)
library(ggplot2)
library(sf)
library(maps)
library(sfheaders)

PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))
#df <- PHONEME_cossim

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
  
  df <- PHONEME_cossim
  land_penalty = 4.44
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

PHONEME_cossim <- compute_shortest_path_df(PHONEME_cossim, ref_coords1, nodes, edges, land_penalty = 44.18)

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

print(PHONEME_cossim$plot[[14]])

land_segments_sf$crosses_land <- TRUE
sea_segments_sf$crosses_land  <- FALSE

main_path_sf <- rbind(land_segments_sf, sea_segments_sf)


connector_segments <- list()

for (i in seq_len(nrow(df))) {
  # Start connector
  if (!is.null(df$land_geom_start[[i]])) {
    connector_segments <- append(connector_segments, list(
      st_sf(geometry = st_sfc(df$land_geom_start[[i]], crs = 4326), crosses_land = TRUE)
    ))
  }
  if (!is.null(df$sea_geom_start[[i]])) {
    connector_segments <- append(connector_segments, list(
      st_sf(geometry = st_sfc(df$sea_geom_start[[i]], crs = 4326), crosses_land = FALSE)
    ))
  }
  
  # End connector
  if (!is.null(df$land_geom_end[[i]])) {
    connector_segments <- append(connector_segments, list(
      st_sf(geometry = st_sfc(df$land_geom_end[[i]], crs = 4326), crosses_land = TRUE)
    ))
  }
  if (!is.null(df$sea_geom_end[[i]])) {
    connector_segments <- append(connector_segments, list(
      st_sf(geometry = st_sfc(df$sea_geom_end[[i]], crs = 4326), crosses_land = FALSE)
    ))
  }
}

connector_sf <- do.call(rbind, connector_segments)


arrow_segments <- list()

for (i in seq_len(nrow(df))) {
  if (!is.null(df$land_geom_start[[i]])) {
    arrow_segments <- append(arrow_segments, list(
      st_sf(geometry = st_sfc(df$land_geom_start[[i]], crs = 4326), crosses_land = TRUE)
    ))
  }
  if (!is.null(df$sea_geom_start[[i]])) {
    arrow_segments <- append(arrow_segments, list(
      st_sf(geometry = st_sfc(df$sea_geom_start[[i]], crs = 4326), crosses_land = FALSE)
    ))
  }
  if (!is.null(df$land_geom_end[[i]])) {
    arrow_segments <- append(arrow_segments, list(
      st_sf(geometry = st_sfc(df$land_geom_end[[i]], crs = 4326), crosses_land = TRUE)
    ))
  }
  if (!is.null(df$sea_geom_end[[i]])) {
    arrow_segments <- append(arrow_segments, list(
      st_sf(geometry = st_sfc(df$sea_geom_end[[i]], crs = 4326), crosses_land = FALSE)
    ))
  }
}

arrow_sf <- do.call(rbind, arrow_segments)


main_path_sf$source   <- "main"
connector_sf$source   <- "connector"

full_tree_sf <- rbind(main_path_sf, connector_sf)

full_tree_lines <- full_tree_sf %>%
  st_cast("MULTILINESTRING") %>%
  st_cast("LINESTRING")

arrow_main <- full_tree_lines %>%
  filter(source == "main") %>%          # ← only main path, not connectors
  mutate(row_id = row_number()) %>%
  group_by(row_id) %>%
  group_modify(~ {
    coords <- st_coordinates(.x)
    if (nrow(coords) < 2) return(.x)
    start_coord <- coords[1, c("X", "Y")]
    end_coord   <- coords[nrow(coords), c("X", "Y")]
    new_geom <- st_sfc(st_linestring(rbind(start_coord, end_coord)), crs = st_crs(.x))
    st_set_geometry(.x, new_geom)
  }) %>%
  ungroup() %>%
  select(-row_id) %>%
  st_as_sf()



ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_sf(data = full_tree_sf, size = 1.5, color = "black") +
  geom_sf(data = arrow_sf,
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "black") +
  geom_sf(data = arrow_main,
        arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
        color = "black") +
  geom_point(data = PHONEME_cossim,aes(x = longitude, y = latitude), 
             size = 3, shape = 21) +
  #geom_text(data = PHONEME_cossim,
   #         aes(x = longitude, y = latitude, label = language),
    #        size = 3, hjust = -0.1, vjust = -0.5) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal() +
  labs(title = "Phonemic Historical Routes Waypoint System",
       color = "Travel Mode",
       x = 'Longitude',
       y = 'Latitude')


arrow_plot <- ggplot() +
  # Main path lines — no arrows
  geom_sf(data = full_tree_sf,
          mapping = aes(geometry = geometry),
          color = "black", linewidth = 1) +
  
  # Connector segments — with arrows only
  geom_sf(data = arrow_sf,
          mapping = aes(geometry = geometry),
          arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
          color = "black", linewidth = 1) +
  
  coord_sf(xlim = c(116, 127), ylim = c(4, 21)) +
  theme_minimal()

print(arrow_plot)
saveRDS(arrow_plot, file = here("data", "phoneme_waypoint_plot.rds"))


# ----- linear model -------------------------------------------

lm_span <- lm(cossim_span ~ geodist_H1_span, data = PHONEME_cossim)
summary_lm_span <- summary(lm_span)

slope <- round(coef(lm_span)[2],5)
coefficient <- round(coef(lm_span)[1],5)
r_squared <- round(summary_lm_span$r.squared,3)
linear_model_equation <- paste0("Y = ", slope, "x + ", coefficient)

ggplot(data = PHONEME_cossim, aes(x = geodist_H1_span, y = cossim_span)) +
  geom_point() +
  geom_smooth(method = 'lm',se = FALSE) +
  theme_bw() +
  labs(title = 'Linear Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity')

# ---- exponential model -------------------------------------

library(rethinking)

# Standardize predictors
PHONEME_cossim$dist_std <- standardize(PHONEME_cossim$geodist_H1_span)

# Model: exponential decay
m_exp <- rethinking::map(
  alist(
    cossim_span ~ dnorm(mu, sigma),
    mu <- a * exp(-b * dist_std),
    a ~ dnorm(0.5, 0.2),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = PHONEME_cossim
)

# Summarize
precis(m_exp)

# Generate predictions
dist_seq <- seq(from = min(PHONEME_cossim$dist_std), to = max(PHONEME_cossim$dist_std), length.out = 100)
preds <- link(m_exp, data = data.frame(dist_std = dist_seq))
mu_mean <- apply(preds, 2, mean)
mu_PI <- apply(preds, 2, PI)
mu_PI_t <- t(mu_PI)
# Plot
scatter_df <- PHONEME_cossim

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
  labs(title = 'Exponential Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity') + 
  theme_bw()


# Extract coefficients
coefs <- coef(m_exp)
a <- round(coefs["a"], 5)
b <- round(coefs["b"], 5)

# Write out the equation
exp_eq <- paste0("μ = ", a, " * exp(-", b, " * x)")
print(exp_eq)

# ----- spline regression --------

library(splines)

spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)
x_vals <- PHONEME_cossim$geodist_H1_span
y_spline <- predict(spline_fit, newdata = data.frame(geodist_H1_span = x_vals))

ggplot(PHONEME_cossim, aes(x = geodist_H1_span, y = cossim_span)) +
  geom_point(color = "gray") +
  geom_line(aes(y = y_spline), color = "blue", linewidth = 1.2) +
  theme_bw() +
  labs(title = 'Spline Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity')

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


spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)

# Linear model predictions
PHONEME_cossim$pred_linear <- predict(lm_span)

# Exponential model predictions using link()
exp_preds <- link(m_exp)
PHONEME_cossim$pred_exp <- apply(exp_preds, 2, mean)

# RMSE calculation
rmse_linear <- sqrt(mean((PHONEME_cossim$cossim_span - PHONEME_cossim$pred_linear)^2))
rmse_exp <- sqrt(mean((PHONEME_cossim$cossim_span - PHONEME_cossim$pred_exp)^2))
rmse_spline <- sqrt(mean((PHONEME_cossim$cossim_span - y_spline)^2))

p_linear <- summary(lm_span)$coefficients[2, 4]  # Slope p-value
p_spline <- pf(summary(spline_fit)$fstatistic[1],
               summary(spline_fit)$fstatistic[2],
               summary(spline_fit)$fstatistic[3],
               lower.tail = FALSE)

exp_post <- extract.samples(m_exp)
pr_b_gt_0 <- mean(exp_post$b > 0)  # Replace 'b' with your actual slope parameter name

# Compare

cossim_span_range <- max(PHONEME_cossim$cossim_span) - min(PHONEME_cossim$cossim_span)

rmse_df <- data.frame(
  Model = c("Linear", "Exponential", "Spline"),
  RMSE = round(c(rmse_linear, rmse_exp, rmse_spline), 5),
  normalized_rmse = round(c(rmse_linear, rmse_exp, rmse_spline) / cossim_span_range, 5),
  p_value = c(round(p_linear, 5), NA, round(p_spline, 5)),
  posterior_prob = c(NA, round(pr_b_gt_0, 3), NA)
)

print(rmse_df)


# ---- dissimilarity matrix ----------------------------------

cosine_matrix_phil <- cosine_matrix[ph_lang, ph_lang]
lang_order <- rownames(cosine_matrix_phil)

cosine_matrix_phil <- 1-cosine_matrix_phil

PHONEME_diss_matrix <- cosine_matrix_phil

# ----- distance matrix --------------------------------------
ph_lang <- PHOIBLEdf %>% 
  filter(Language_type == 'Philippine Language') %>% 
  pull(language)

phil_df <- PHOIBLEdf %>%
  filter(Language_type == "Philippine Language") %>%
  mutate(
    start_coords = map2(longitude, latitude, ~ c(.x, .y)),
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

PHONEME_dist_matrix <- dist_matrix[ph_lang, ph_lang]
PHONEME_dist_matrix[is.na(PHONEME_dist_matrix)] <- 0


# ---- plot matrices -----------------------------------

melt_phoneme_dist_matrix <- melt(PHONEME_dist_matrix)
melt_phoneme_diss_matrix <- melt(PHONEME_diss_matrix)

dist_matrix_p <- ggplot(melt_phoneme_dist_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are square

diss_matrix_p  <- ggplot(melt_phoneme_diss_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") + # Creates the colored tiles
  scale_fill_gradient(low = "yellow", high = "red") + # Customizes the colors
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Rotates x-axis labels
  coord_fixed() # Ensures cells are square

dist_matrix_p + diss_matrix_p 


# ---- mantel test ----------------------------------------
library(vegan)

# Convert to distance objects
x_dist <- as.dist(PHONEME_dist_matrix)
y_dist <- as.dist(PHONEME_diss_matrix)

# Run Mantel test
mantel_result <- mantel(x_dist, y_dist, method = "spearman", permutations = 999)

print(mantel_result)

library(ggplot2)

# Convert matrices to vectors
x_vec <- as.vector(as.dist(PHONEME_dist_matrix))
y_vec <- as.vector(as.dist(PHONEME_diss_matrix))




# Plot
ggplot(data.frame(x = x_vec, y = y_vec), aes(x = x, y = y)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  theme_bw() +
  labs(
    title = "Mantel Test: Phoneme Distance vs. Dissimilarity",
    x = "Relative Migration Pairwise Distance (km)",
    y = "Phonemic Dissimilarity"
  )


# ----- saving dataframes -----------------------------------

write.csv(cosine_matrix, file = here("data", "PHONEME_cosine_matrix.csv"), row.names = TRUE)

PHONEME_cossim <- df_span %>% 
  left_join(df_jap,by = 'language') %>% 
  left_join(df_eng,by = 'language') %>% 
  left_join(df_unr,by = 'language') %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

write.csv(PHONEME_cossim, file = here("data", "PHONEME_cossim.csv"), row.names = TRUE)
