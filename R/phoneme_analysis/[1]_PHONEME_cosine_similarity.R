# =============================================================================
# [1] Phoneme Analysis — Weighted cosine similarity
# Computes the IDF-weighted cosine-similarity matrix over phoneme inventories,
# extracts each Philippine language's similarity to Spanish / Japanese / English
# and its mean similarity to the unrelated controls, plus a 0-to-max min-max
# normalized Spanish similarity (cossim_span_norm) for [4]'s regression, and
# writes the matrix + per-language scores for the downstream [2]–[7] analyses
# and EEMS plotting.
#
# Input:   data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen_austronesian.csv
# Outputs: data/PHONEME_cosine_matrix.csv, data/PHONEME_cossim.csv,
#          data/base_plot_phoneme_cosine.rds (cosine base map, for [7])
# Next:    [2]_PHONEME_cosine_distribution_analysis.R,
#          [3]_PHONEME_network_distance.R (then regression / MMRR)
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(here)

# ---- Prepping Globals -------------------------------------------------------
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"))

ph_lang <- RUHLENdf |>
  filter(Language_type == 'Philippine Language') |>
  pull(language)

int_lang <- RUHLENdf |>
  filter(Language_type == 'Language of Interest') |>
  pull(language)

unr_lang <- RUHLENdf |>
  filter(Language_type == 'Unrelated Language') |>
  pull(language)

phoneme_cols <- RUHLENdf %>%
  dplyr::select(-language, -source, -iso6393, -Language_type, -latitude, -longitude)

phoneme_cols <- colnames(phoneme_cols)

phoneme_freq <- read_csv(here("data", "phoneme_freq_ruhlen_austronesian.csv"))

# ---- sensitivity analysis -----

# Most recent coverage


# Least coverage


# ----- cosine similarity -----------------------------------------------------

calculate_weighted_cosine_similarity <- function(RUHLENdf, phoneme_freq, phoneme_cols, id_col = "language") {

  # Reorder the frequency table to the phoneme-column order so the IDF weights
  # line up with the binary matrix columns.
  aligned_freq <- phoneme_freq %>%
    dplyr::filter(phoneme %in% phoneme_cols) %>%
    dplyr::arrange(match(phoneme, phoneme_cols))

  idf_weights <- aligned_freq$IDF

  binary_data <- RUHLENdf %>%
    dplyr::select(dplyr::all_of(phoneme_cols)) %>%
    as.matrix()

  language_ids <- RUHLENdf[[id_col]]

  # Scale each phoneme column by its IDF weight before the similarity loop.
  weighted_data <- sweep(binary_data, 2, idf_weights, FUN = "*")

  n_languages <- nrow(weighted_data)
  cosine_matrix <- matrix(0, nrow = n_languages, ncol = n_languages,
                          dimnames = list(language_ids, language_ids))

  # epsilon guards against division by zero for a language with no phonemes.
  epsilon <- 1e-9

  for (i in 1:n_languages) {
    for (j in i:n_languages) {

      vec_a <- weighted_data[i, ]
      vec_b <- weighted_data[j, ]

      dot_product <- sum(vec_a * vec_b)

      magnitude_a <- sqrt(sum(vec_a^2))
      magnitude_b <- sqrt(sum(vec_b^2))

      denominator <- (magnitude_a * magnitude_b) + epsilon

      score <- dot_product / denominator

      cosine_matrix[i, j] <- score
      cosine_matrix[j, i] <- score  # symmetric
    }
  }

  return(cosine_matrix)
}

attested_phonemes <- phoneme_freq$phoneme

cosine_matrix <- calculate_weighted_cosine_similarity(
  RUHLENdf,
  phoneme_freq,
  attested_phonemes,
  id_col = "language")

# Order rows/cols by language type for readable matrices/heatmaps
ordered_languages <- RUHLENdf %>%
  arrange(Language_type) %>%
  pull(language)

cosine_matrix <- cosine_matrix[ordered_languages, ordered_languages]


cosine_matrix['Tagalog','Spanish'] 
cosine_matrix['Tagalog','English'] 
cosine_matrix['Tagalog','Japanese'] 

# ---- Per-language similarity to each baseline -------------------------------
df_span <- cosine_matrix[ph_lang, "Spanish"] |>
  enframe(name = "language", value = "cossim_span")

df_jap <- cosine_matrix[ph_lang, "Japanese"] |>
  enframe(name = "language", value = "cossim_jap")

df_eng <- cosine_matrix[ph_lang, "English"] |>
  enframe(name = "language", value = "cossim_eng")

df_unr <- rowMeans(cosine_matrix[ph_lang, unr_lang]) |>
  enframe(name = "language", value = "cossim_unr")

# ---- Diagnostic: full cosine-similarity heatmap -----------------------------
melted_matrix <- melt(cosine_matrix)

ggplot(melted_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

# ----- saving dataframes -----------------------------------------------------
write.csv(cosine_matrix, file = here("data", "PHONEME_cosine_matrix.csv"), row.names = TRUE)

PHONEME_cossim <- RUHLENdf |>
  filter(Language_type == 'Philippine Language') |>
  dplyr::select(language, latitude, longitude) |>
  left_join(df_span, by = 'language') |>
  left_join(df_jap,  by = 'language') |>
  left_join(df_eng,  by = 'language') |>
  left_join(df_unr,  by = 'language')

# Min-max normalized Spanish similarity, floor fixed at 0 (same 0-to-max
# convention as the map colour scale below) rather than the empirical min, so
# 0 reads as "no similarity at all" and 1 as the most Spanish-like PH language
# observed. Used downstream by [4]'s regression so a slope reads as "% of the
# observed range per unit distance" instead of a raw cosine-similarity delta.
cossim_span_max <- max(PHONEME_cossim$cossim_span, na.rm = TRUE)
PHONEME_cossim <- PHONEME_cossim |>
  mutate(cossim_span_norm = cossim_span / cossim_span_max)

write.csv(PHONEME_cossim, file = here("data", "PHONEME_cossim.csv"), row.names = TRUE)

# ---- Cosine base map --------------------------------------------------------
# The per-language scores on a Philippines map, saved for [7] to overlay the
# waypoint routes onto. Built here rather than in [6] so the geographic
# distribution can be read straight off [1], before FEEMS or the arrows exist.
world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

global_lim <- c(0, cossim_span_max)

base_plot_cosine <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray95", color = "gray70") +
  geom_point(data = PHONEME_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 4, alpha = 0.7) +
  geom_point(data = PHONEME_cossim, aes(x = longitude, y = latitude),
             size = 4, shape = 21, color = "black") +
  scale_color_gradient(low = "white", high = "navy", limits = global_lim) +
  guides(color = guide_colorbar(title = "Cosine Similarity",
                                title.position = "top", title.hjust = 0.5)) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_x_continuous(breaks = seq(115, 130, by = 2)) +
  scale_y_continuous(breaks = seq(4, 22, by = 2)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

base_plot_cosine
saveRDS(base_plot_cosine, file = here("data", "base_plot_phoneme_cosine.rds"))
