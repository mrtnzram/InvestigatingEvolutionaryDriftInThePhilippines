# =============================================================================
# [1] Phoneme Analysis — Weighted cosine similarity
# Computes the IDF-weighted cosine-similarity matrix over phoneme inventories,
# extracts each Philippine language's similarity to Spanish / Japanese / English
# and its mean similarity to the unrelated controls, and writes the matrix +
# per-language scores for the downstream [2]–[6] analyses and EEMS plotting.
#
# Input:   data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen.csv
# Outputs: data/PHONEME_cosine_matrix.csv, data/PHONEME_cossim.csv
# Next:    [2]_PHONEME_cosine_distribution_analysis.R,
#          [3]_PHONEME_network_distance.R (then regression / mantel)
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
  select(-language, -source, -iso6393, -Language_type, -latitude, -longitude)

phoneme_cols <- colnames(phoneme_cols)

phoneme_freq <- read_csv(here("data", "phoneme_freq_ruhlen_austronesian.csv"))

# ---- sensitivity analysis -----

# Most recent coverage


# Least coverage



# ----- cosine similarity -----------------------------------------------------

calculate_weighted_cosine_similarity <- function(RUHLENdf, phoneme_freq, phoneme_cols, id_col = "language") {

  # Extract and align the binary data and IDF weights
  # Ensure the phoneme frequencies are in the same order as the phoneme columns
  aligned_freq <- phoneme_freq %>%
    dplyr::filter(phoneme %in% phoneme_cols) %>%
    dplyr::arrange(match(phoneme, phoneme_cols))

  idf_weights <- aligned_freq$IDF

  # Extract the binary phoneme data
  binary_data <- RUHLENdf %>%
    dplyr::select(dplyr::all_of(phoneme_cols)) %>%
    as.matrix() # Convert to a matrix for faster calculations

  # Extract language IDs for matrix naming
  language_ids <- RUHLENdf[[id_col]]

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

attested_phonemes <- phoneme_freq_austronesian$phoneme

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
  select(language, latitude, longitude) |>
  left_join(df_span, by = 'language') |>
  left_join(df_jap,  by = 'language') |>
  left_join(df_eng,  by = 'language') |>
  left_join(df_unr,  by = 'language')

write.csv(PHONEME_cossim, file = here("data", "PHONEME_cossim.csv"), row.names = TRUE)
