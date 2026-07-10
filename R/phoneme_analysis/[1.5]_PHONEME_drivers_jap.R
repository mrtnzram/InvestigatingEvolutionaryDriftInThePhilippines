# =============================================================================
# [1.5] Phoneme drivers for japanese similarity
# investigates the phoneme drivers for Japanese similarity based on 
# IDF weights and ph presence
#
# Input:   data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen.csv
# Outputs: 
# Next:    [2]_PHONEME_cosine_distribution_analysis.R,
#          [3]_PHONEME_network_distance.R (then regression / mantel)
# =============================================================================
library(tidyverse)
library(here)

# --- look into phoneme drivers for japan ----------------

lang_col    <- "language"     # language-name column in RUHLENdf
japanese_id <- "Japanese"     # exact label as it appears in that column

# ---- 1. idf vector keyed by phoneme ---------------------------------------
idf_vec <- phoneme_freq |>
  select(phoneme, IDF) |>
  distinct() |>
  deframe()

# ---- 2. binary phoneme matrix from RUHLENdf -------------------------------
phoneme_mat <- RUHLENdf |>
  column_to_rownames(lang_col) |>
  select(all_of(feat)) |>
  as.matrix()
phoneme_mat[is.na(phoneme_mat)] <- 0
phoneme_mat <- (phoneme_mat > 0) * 1

# ---- 3. IDF-weight rows, unit-normalise -----------------------------------
V         <- sweep(phoneme_mat, 2, idf_vec, `*`)
row_norms <- sqrt(rowSums(V^2))
row_norms[row_norms == 0] <- NA_real_
U         <- V / row_norms

# ---- 4. Philippine centroid + Japanese decomposition ----------------------

# ---- adapter: set to df_unr's actual column names -------------------------
# ---- subset: PH languages within 1 SD of the mean cossim_unr -----------
mean_unr <- mean(df_unr$cossim_unr, na.rm = TRUE)
sd_unr   <- sd(df_unr$cossim_unr,   na.rm = TRUE)

ph_lang_core <- df_unr |>
  filter(
    .data$cossim_unr >= mean_unr - sd_unr,
    .data$cossim_unr <= mean_unr + sd_unr
  ) |>
  pull(.data$language) |>
  intersect(ph_lang)          # guard: keep only rows that are actually PH langs present in U

stopifnot(length(ph_lang_core) > 0, all(ph_lang_core %in% rownames(U)))

n_dropped <- length(ph_lang) - length(ph_lang_core)
message(n_dropped, " of ", length(ph_lang),
        " Philippine languages excluded as cossim_unr outliers (>1 SD from mean)")

# ---- centroid computed only on the trimmed core -------------------------
phil_centroid <- colMeans(U[ph_lang_core, , drop = FALSE])
jp_contrib    <- U[japanese_id, ] * phil_centroid    # sums to mean cos(JP, core PH)

japanese_drivers <- tibble(
  phoneme       = names(jp_contrib),
  IDF           = idf_aligned,
  ph_prevalence = colMeans(phoneme_mat[ph_lang, , drop = FALSE]),  # frac PH langs with p
  jp_has        = phoneme_mat[japanese_id, ] == 1,
  contribution  = jp_contrib,
  share         = jp_contrib / sum(jp_contrib)
) |>
  filter(contribution > 0) |>
  arrange(desc(IDF))

mean_sim_jp <- sum(jp_contrib)

mean_sim_to_ph        <- as.vector(U %*% phil_centroid)
names(mean_sim_to_ph) <- rownames(U)

jp_rank <- tibble(language = names(mean_sim_to_ph),
                  mean_sim = mean_sim_to_ph) |>
  arrange(desc(mean_sim)) |>
  mutate(rank = row_number())

mean_sim_jp
jp_rank |> filter(language == japanese_id)
print(japanese_drivers, n = 25)



# MAX similarity breakdown

# ---- identify the Philippine max-similarity language from df_jap ----------
# adjust the two column names to df_jap's schema:
jap_lang_col <- "language"    # Philippine language id in df_jap
jap_sim_col  <- "cossim_jap"  # Japanese-similarity column in df_jap

ph_max_row  <- df_jap |> slice_max(.data[[jap_sim_col]], n = 1, with_ties = FALSE)
ph_max_lang <- ph_max_row[[jap_lang_col]]
ph_max_val  <- ph_max_row[[jap_sim_col]]

stopifnot(ph_max_lang %in% rownames(U))

# ---- pairwise Japanese-vs-maxlang decomposition ---------------------------
# contribution_p = U[JP, p] * U[maxlang, p]; sums to cos(JP, maxlang)
jp_max_contrib <- U[japanese_id, ] * U[ph_max_lang, ]

japanese_drivers_max <- tibble(
  phoneme      = names(jp_max_contrib),
  IDF          = idf_aligned,
  jp_has       = phoneme_mat[japanese_id, ] == 1,
  maxlang_has  = phoneme_mat[ph_max_lang, ] == 1,
  contribution = jp_max_contrib,
  share        = jp_max_contrib / sum(jp_max_contrib)
) |>
  filter(contribution > 0) |>
  arrange(desc(contribution))

# reconstructed cosine should equal df_jap's stored value for this language
recomputed_max <- sum(jp_max_contrib)

# ---- side-by-side: centroid drivers vs max-pair drivers -------------------
driver_compare <- full_join(
  japanese_drivers     |> select(phoneme, IDF, centroid_contrib = contribution),
  japanese_drivers_max |> select(phoneme,      max_contrib      = contribution),
  by = "phoneme"
) |>
  mutate(across(c(centroid_contrib, max_contrib), \(x) replace_na(x, 0))) |>
  arrange(desc(max_contrib))

# ---- inspect --------------------------------------------------------------
c(ph_max_lang = ph_max_lang, stored = ph_max_val, recomputed = recomputed_max)
print(japanese_drivers_max, n = 25)
print(driver_compare,       n = 25)











# =============================================================================
# [1.5] Phoneme drivers for Spanish similarity
# investigates the phoneme drivers for Japanese similarity based on 
# IDF weights and ph presence
#
# Input:   data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen.csv
# Outputs: 
# Next:    [2]_PHONEME_cosine_distribution_analysis.R,
#          [3]_PHONEME_network_distance.R (then regression / mantel)
# =============================================================================


# --- look into phoneme drivers for japan ----------------

lang_col    <- "language"     # language-name column in RUHLENdf
spanish_id <- "Spanish"     # exact label as it appears in that column

# ---- 1. idf vector keyed by phoneme ---------------------------------------
idf_vec <- phoneme_freq |>
  select(phoneme, IDF) |>
  distinct() |>
  deframe()

# ---- 2. binary phoneme matrix from RUHLENdf -------------------------------
phoneme_mat <- RUHLENdf |>
  column_to_rownames(lang_col) |>
  select(all_of(feat)) |>
  as.matrix()
phoneme_mat[is.na(phoneme_mat)] <- 0
phoneme_mat <- (phoneme_mat > 0) * 1

# ---- 3. IDF-weight rows, unit-normalise -----------------------------------
V         <- sweep(phoneme_mat, 2, idf_vec, `*`)
row_norms <- sqrt(rowSums(V^2))
row_norms[row_norms == 0] <- NA_real_
U         <- V / row_norms

# ---- 4. Philippine centroid + Japanese decomposition ----------------------

# ---- adapter: set to df_unr's actual column names -------------------------
# ---- subset: PH languages within 1 SD of the mean cossim_unr -----------
mean_unr <- mean(df_unr$cossim_unr, na.rm = TRUE)
sd_unr   <- sd(df_unr$cossim_unr,   na.rm = TRUE)

ph_lang_core <- df_unr |>
  filter(
    .data$cossim_unr >= mean_unr - sd_unr,
    .data$cossim_unr <= mean_unr + sd_unr
  ) |>
  pull(.data$language) |>
  intersect(ph_lang)          # guard: keep only rows that are actually PH langs present in U

stopifnot(length(ph_lang_core) > 0, all(ph_lang_core %in% rownames(U)))

n_dropped <- length(ph_lang) - length(ph_lang_core)
message(n_dropped, " of ", length(ph_lang),
        " Philippine languages excluded as cossim_unr outliers (>1 SD from mean)")

# ---- centroid computed only on the trimmed core -------------------------
phil_centroid <- colMeans(U[ph_lang_core, , drop = FALSE])
jp_contrib    <- U[japanese_id, ] * phil_centroid    # sums to mean cos(JP, core PH)

japanese_drivers <- tibble(
  phoneme       = names(jp_contrib),
  IDF           = idf_aligned,
  ph_prevalence = colMeans(phoneme_mat[ph_lang, , drop = FALSE]),  # frac PH langs with p
  jp_has        = phoneme_mat[japanese_id, ] == 1,
  contribution  = jp_contrib,
  share         = jp_contrib / sum(jp_contrib)
) |>
  filter(contribution > 0) |>
  arrange(desc(IDF))

mean_sim_jp <- sum(jp_contrib)

mean_sim_to_ph        <- as.vector(U %*% phil_centroid)
names(mean_sim_to_ph) <- rownames(U)

jp_rank <- tibble(language = names(mean_sim_to_ph),
                  mean_sim = mean_sim_to_ph) |>
  arrange(desc(mean_sim)) |>
  mutate(rank = row_number())

mean_sim_jp
jp_rank |> filter(language == japanese_id)
print(japanese_drivers, n = 25)



# MAX similarity breakdown

# ---- identify the Philippine max-similarity language from df_jap ----------
# adjust the two column names to df_jap's schema:
jap_lang_col <- "language"    # Philippine language id in df_jap
jap_sim_col  <- "cossim_jap"  # Japanese-similarity column in df_jap

ph_max_row  <- df_jap |> slice_max(.data[[jap_sim_col]], n = 1, with_ties = FALSE)
ph_max_lang <- ph_max_row[[jap_lang_col]]
ph_max_val  <- ph_max_row[[jap_sim_col]]

stopifnot(ph_max_lang %in% rownames(U))

# ---- pairwise Japanese-vs-maxlang decomposition ---------------------------
# contribution_p = U[JP, p] * U[maxlang, p]; sums to cos(JP, maxlang)
jp_max_contrib <- U[japanese_id, ] * U[ph_max_lang, ]

japanese_drivers_max <- tibble(
  phoneme      = names(jp_max_contrib),
  IDF          = idf_aligned,
  jp_has       = phoneme_mat[japanese_id, ] == 1,
  maxlang_has  = phoneme_mat[ph_max_lang, ] == 1,
  ph_prevalence = colMeans(phoneme_mat[ph_lang, , drop = FALSE]),  # frac PH langs with p
  contribution = jp_max_contrib,
  share        = jp_max_contrib / sum(jp_max_contrib)
) |>
  filter(contribution > 0) |>
  arrange(desc(contribution))

# reconstructed cosine should equal df_jap's stored value for this language
recomputed_max <- sum(jp_max_contrib)

# ---- side-by-side: centroid drivers vs max-pair drivers -------------------
driver_compare <- full_join(
  japanese_drivers     |> select(phoneme, IDF, centroid_contrib = contribution),
  japanese_drivers_max |> select(phoneme,      max_contrib      = contribution),
  by = "phoneme"
) |>
  mutate(across(c(centroid_contrib, max_contrib), \(x) replace_na(x, 0))) |>
  arrange(desc(max_contrib))

# ---- inspect --------------------------------------------------------------
c(ph_max_lang = ph_max_lang, stored = ph_max_val, recomputed = recomputed_max)
print(japanese_drivers_max, n = 25)
print(driver_compare,       n = 25)











# ------ Cross Reference with phoneme names ---------------

# find the header row (file has a ~13-line preamble), then read the TSV table
raw   <- read_lines(here('data','pnas_1424033112_sd02.txt'))
hdr_i <- which(str_starts(raw, "Column\t"))

pnas_key <- read_tsv(I(raw[hdr_i:length(raw)]), show_col_types = FALSE) |>
  transmute(
    si_column = Column,
    phoneme_id = str_c("phoneme_", str_pad(Column - 9L, 3, pad = "0")),
    ipa        = Phoneme,
    global_n   = Number_of_occurrences,
    class = case_when(
      Consonant          == 1 ~ "consonant",
      Vowel              == 1 ~ "vowel",
      Modified_consonant == 1 ~ "mod_consonant",
      Modified_Vowel     == 1 ~ "mod_vowel",
      Click              == 1 ~ "click"
    )
  )

# attach IPA to Japanese driver table
japanese_drivers_max |>
  left_join(pnas_key, by = c("phoneme" = "phoneme_id")) |>
  select(phoneme, ipa, class, global_n, IDF, ph_prevalence, contribution, share)
