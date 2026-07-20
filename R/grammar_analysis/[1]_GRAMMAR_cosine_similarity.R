# =============================================================================
# [1] Grammar Analysis — Weighted cosine similarity
#
# Computes the IDF-weighted cosine-similarity matrix over GRAMBANK feature
# inventories, extracts each Philippine language's similarity to Spanish /
# Japanese / English and its mean similarity to the unrelated controls.
#
# Adaptation note (vs. phoneme [1]): GRAMBANK features are multi-level
# categorical, not binary, so the >2-level columns are one-hot encoded first
# (fastDummies) and IDF weights are looked up per (feature, value) pair from
# gramfeature_freq.csv (unlike phoneme's per-column IDF sweep). Once encoded,
# the weighting (IDF substituted for each present feature-value) and cosine
# calculation are the same as phoneme [1].
#
# Input:   data/GRAMBANKdf_full.csv, data/gramfeature_freq.csv
# Outputs: data/GRAMMAR_cosine_matrix.csv, data/GRAMMAR_cossim.csv
# Next:    [2]_GRAMMAR_cosine_distribution_analysis.R,
#          [3]_GRAMMAR_network_distance.R (then regression / MMRR)
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(ggplot2)
library(here)
library(fastDummies)
library(purrr)
library(stringr)
library(tibble)

# ---- Prepping Globals -------------------------------------------------------
GRAMBANKdf_PH <- read_csv(here("data", "GRAMBANKdf_full.csv"), show_col_types = FALSE)

ph_lang <- GRAMBANKdf_PH |>
  filter(Language_Type == 'Philippine Language') |>
  pull(language)

int_lang <- GRAMBANKdf_PH |>
  filter(Language_Type == 'Language of Interest') |>
  pull(language)

unr_lang <- GRAMBANKdf_PH |>
  filter(Language_Type == 'Unrelated Language') |>
  pull(language)

feature_cols <- GRAMBANKdf_PH %>%
  dplyr::select(-Language_ID, -language, -Family_name, -Macroarea, -Longitude, -Latitude, -Language_Type)

feature_cols <- colnames(feature_cols)

gramfeature_freq <- read_csv(here("data", "gramfeature_freq.csv"), show_col_types = FALSE)

# ---- one-hot encode categorical columns --------------------------------------
n_levels <- map_int(
  feature_cols,
  \(col) GRAMBANKdf_PH[[col]] |> na.omit() |> unique() |> length()
)

binary_cols      <- feature_cols[n_levels <= 2]   # 0/1  -> already valid; keep
categorical_cols <- feature_cols[n_levels >  2]    # >2 levels -> expand to dummies

cat(sprintf(
  "Binary (kept):       %d cols\nCategorical (OHE'd): %d cols\n",
  length(binary_cols), length(categorical_cols)
))

GRAMBANKdf_PH_ohe <- GRAMBANKdf_PH |>
  dummy_cols(
    select_columns          = categorical_cols,
    remove_selected_columns = TRUE,
    remove_first_dummy      = FALSE,   # full OHE - keep all levels
    ignore_na                = TRUE     # NA row -> NA across all dummies
  )

ohe_pattern <- str_c("^(", str_c(categorical_cols, collapse = "|"), ")_")
ohe_cols    <- str_subset(names(GRAMBANKdf_PH_ohe), ohe_pattern)

feature_cols_ohe <- c(binary_cols, ohe_cols)

cat(sprintf(
  "Feature matrix: %d binary + %d OHE cols = %d total features\n",
  length(binary_cols), length(ohe_cols), length(feature_cols_ohe)
))

GRAMBANKdf_PH <- GRAMBANKdf_PH_ohe

# ----- cosine similarity -----------------------------------------------------

calculate_weighted_cosine_similarity <- function(GRAMBANKdf_PH, gramfeature_freq, feature_cols, ohe_cols = character(0), id_col = "language") {

  long_data <- GRAMBANKdf_PH %>%
    select(all_of(feature_cols)) %>%
    mutate(language = GRAMBANKdf_PH[[id_col]]) %>%
    pivot_longer(cols = -language, names_to = "feature", values_to = "value")

  # OHE dummy columns (e.g. "GB065_1") don't exist by that name in
  # gramfeature_freq, which is keyed by the ORIGINAL feature + its raw value
  # (e.g. "GB065", 1) — a naive join here leaves every OHE'd feature 100% NA,
  # silently dropping it from the cosine calculation entirely. Parse the dummy
  # name back to (base feature, level) and weight it as value * IDF(level) —
  # the standard one-hot + IDF scheme: a language's vector gets a weighted
  # spike in the dimension matching its actual category, zero elsewhere.
  # Binary (non-OHE) columns keep the original substitution scheme: IDF of
  # whichever value (0/1) was actually observed.
  is_ohe <- long_data$feature %in% ohe_cols
  long_data <- long_data %>%
    mutate(
      lookup_feature = if_else(is_ohe, str_remove(feature, "_[^_]+$"), feature),
      lookup_value   = if_else(is_ohe, as.numeric(str_extract(feature, "(?<=_)[^_]+$")), value)
    )

  weighted_long <- long_data %>%
    left_join(gramfeature_freq, by = c("lookup_feature" = "feature", "lookup_value" = "value")) %>%
    mutate(weighted_value = if_else(feature %in% ohe_cols, value * IDF, IDF))

  weighted_data <- weighted_long %>%
    select(language, feature, weighted_value) %>%
    pivot_wider(names_from = feature, values_from = weighted_value) %>%
    column_to_rownames("language") %>%
    as.matrix()

  language_ids <- GRAMBANKdf_PH[[id_col]]

  n_languages <- nrow(weighted_data)
  cosine_matrix <- matrix(0, nrow = n_languages, ncol = n_languages,
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

        denominator <- (magnitude_a * magnitude_b) + epsilon
        score <- dot_product / denominator
      }

      cosine_matrix[i, j] <- score
      cosine_matrix[j, i] <- score
    }
  }

  return(cosine_matrix)
}

cosine_matrix <- calculate_weighted_cosine_similarity(
  GRAMBANKdf_PH,
  gramfeature_freq,
  feature_cols_ohe,
  ohe_cols = ohe_cols,
  id_col = "language")

# Order rows/cols by language type for readable matrices/heatmaps
ordered_languages <- GRAMBANKdf_PH %>%
  arrange(Language_Type) %>%
  pull(language)

cosine_matrix <- cosine_matrix[ordered_languages, ordered_languages]

cosine_matrix['Filipino', 'Spanish']
cosine_matrix['Filipino', 'English']
cosine_matrix['Filipino', 'Japanese']

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
write.csv(cosine_matrix, file = here("data", "GRAMMAR_cosine_matrix.csv"), row.names = TRUE)

GRAMMAR_cossim <- GRAMBANKdf_PH |>
  filter(Language_Type == 'Philippine Language') |>
  dplyr::select(language, latitude = Latitude, longitude = Longitude) |>
  left_join(df_span, by = 'language') |>
  left_join(df_jap,  by = 'language') |>
  left_join(df_eng,  by = 'language') |>
  left_join(df_unr,  by = 'language')

write.csv(GRAMMAR_cossim, file = here("data", "GRAMMAR_cossim.csv"), row.names = TRUE)
