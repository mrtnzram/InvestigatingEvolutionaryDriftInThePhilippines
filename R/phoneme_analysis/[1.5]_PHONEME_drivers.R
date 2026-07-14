# =============================================================================
# [1.5] Phoneme drivers
# Diagnostic: decomposes the IDF-weighted cosine similarity into per-phoneme
# contributions to see which phonemes drive a Philippine language's similarity to
# the interest-language baselines, and which phonemes carry meaningful vs.
# unmeaningful variance across the Philippine languages.
#
# Sections:
#   1. Cross-reference phoneme names (phoneme_id -> IPA / class, from the PNAS key)
#   2. Japanese phoneme drivers (centroid + max-pair decomposition, worked example)
#   3. Final phoneme drivers table (prevalences + baseline indicators + variance;
#      right-skew leverage on the unrelated null + a bidirectional driver flag)
#
# Inputs:  data/RUHLENdf_PH.csv                    (binary phoneme matrix + Language_type)
#          data/phoneme_freq_ruhlen_austronesian.csv (IDF weights + Austronesian freq)
#          data/PHONEME_cossim.csv                 (per-PH-language baseline cosines)
#          data/pnas_1424033112_sd02.txt           (phoneme_id -> IPA key)
# Output:  data/PHONEME_driver_table.csv
# =============================================================================

library(tidyverse)
library(here)

# =============================================================================
# Section 0 — Load & prep (self-contained; no dependency on [1] being sourced)
# =============================================================================
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"), show_col_types = FALSE)

ph_lang  <- RUHLENdf |> filter(Language_type == "Philippine Language") |> pull(language)
unr_lang <- RUHLENdf |> filter(Language_type == "Unrelated Language")  |> pull(language)
int_lang <- RUHLENdf |> filter(Language_type == "Language of Interest") |> pull(language)
feat     <- RUHLENdf |> select(starts_with("phoneme_")) |> names()

# per-PH-language baseline cosines (replaces [1]'s df_unr / df_jap)
PHONEME_cossim <- read_csv(here("data", "PHONEME_cossim.csv"), show_col_types = FALSE) |>
  select(-any_of("...1"))

# IDF weights + Austronesian prevalence (the same weighting [1] uses for the cosine)
phoneme_freq <- read_csv(here("data", "phoneme_freq_ruhlen_austronesian.csv"),
                         show_col_types = FALSE)
idf_vec     <- phoneme_freq |> select(phoneme, IDF) |> distinct() |> deframe()
idf_aligned <- idf_vec[feat]          # aligned to the phoneme-column order of `feat`

# Global prevalence over the full ~1772-language Ruhlen database (the broadest
# reference, above Austronesian). Ultra-rare phonemes (global_n <= 3) are absent
# from this file; their global prevalence is ~0 and coalesced to 0 below.
phoneme_freq_global <- read_csv(here("data", "phoneme_freq_ruhlen.csv"),
                                show_col_types = FALSE)

# binary phoneme matrix (rows = languages, cols = phonemes)
phoneme_mat <- RUHLENdf |>
  column_to_rownames("language") |>
  select(all_of(feat)) |>
  as.matrix()
phoneme_mat[is.na(phoneme_mat)] <- 0
phoneme_mat <- (phoneme_mat > 0) * 1

# IDF-weight each row, then unit-normalise -> U (so U[a,] . U[b,] = cosine(a, b))
V         <- sweep(phoneme_mat, 2, idf_aligned, `*`)
row_norms <- sqrt(rowSums(V^2))
row_norms[row_norms == 0] <- NA_real_
U         <- V / row_norms


# =============================================================================
# Section 1 — Cross-reference phoneme names
# =============================================================================
# The RUHLEN phoneme_### ids map to IPA symbols via the PNAS SI key. The file has
# a ~13-line preamble, so find the header row ("Column\t...") and read from there.
raw   <- read_lines(here("data", "pnas_1424033112_sd02.txt"))
hdr_i <- which(str_starts(raw, "Column\t"))

pnas_key <- read_tsv(I(raw[hdr_i:length(raw)]), show_col_types = FALSE) |>
  transmute(
    si_column  = Column,
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


# =============================================================================
# Section 2 — Japanese phoneme drivers (worked example)
# =============================================================================
japanese_id <- "Japanese"

# ---- 2a. Centroid decomposition -------------------------------------------
# Compare Japanese against the CENTROID of a trimmed "core" of PH languages
# (those within 1 SD of the mean unrelated similarity, to drop outliers). The
# per-phoneme product U[JP, p] * centroid[p] sums to mean cos(JP, core PH).
mean_unr <- mean(PHONEME_cossim$cossim_unr, na.rm = TRUE)
sd_unr   <- sd(PHONEME_cossim$cossim_unr,   na.rm = TRUE)

ph_lang_core <- PHONEME_cossim |>
  filter(cossim_unr >= mean_unr - sd_unr,
         cossim_unr <= mean_unr + sd_unr) |>
  pull(language) |>
  intersect(ph_lang)          # keep only PH langs actually present in U

stopifnot(length(ph_lang_core) > 0, all(ph_lang_core %in% rownames(U)))

message(length(ph_lang) - length(ph_lang_core), " of ", length(ph_lang),
        " Philippine languages excluded as cossim_unr outliers (>1 SD from mean)")

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

# where does Japanese rank among all languages by mean similarity to the core?
mean_sim_jp    <- sum(jp_contrib)
mean_sim_to_ph <- as.vector(U %*% phil_centroid) |> setNames(rownames(U))
jp_rank <- tibble(language = names(mean_sim_to_ph), mean_sim = mean_sim_to_ph) |>
  arrange(desc(mean_sim)) |>
  mutate(rank = row_number())

mean_sim_jp
jp_rank |> filter(language == japanese_id)
print(japanese_drivers, n = 25)

# ---- 2b. Max-pair decomposition -------------------------------------------
# Break down cos(JP, the single most Japanese-similar PH language) phoneme by
# phoneme. contribution_p = U[JP, p] * U[maxlang, p]; sums to cos(JP, maxlang).
ph_max_row  <- PHONEME_cossim |> slice_max(cossim_jap, n = 1, with_ties = FALSE)
ph_max_lang <- ph_max_row$language
ph_max_val  <- ph_max_row$cossim_jap
stopifnot(ph_max_lang %in% rownames(U))

jp_max_contrib <- U[japanese_id, ] * U[ph_max_lang, ]

japanese_drivers_max <- tibble(
  phoneme       = names(jp_max_contrib),
  IDF           = idf_aligned,
  jp_has        = phoneme_mat[japanese_id, ] == 1,
  maxlang_has   = phoneme_mat[ph_max_lang, ] == 1,
  ph_prevalence = colMeans(phoneme_mat[ph_lang, , drop = FALSE]),
  contribution  = jp_max_contrib,
  share         = jp_max_contrib / sum(jp_max_contrib)
) |>
  filter(contribution > 0) |>
  arrange(desc(contribution))

# reconstructed cosine should match PHONEME_cossim's stored value for this language
recomputed_max <- sum(jp_max_contrib)
c(ph_max_lang = ph_max_lang, stored = ph_max_val, recomputed = recomputed_max)

# ---- 2c. Side-by-side: centroid drivers vs max-pair drivers ---------------
driver_compare <- full_join(
  japanese_drivers     |> select(phoneme, IDF, centroid_contrib = contribution),
  japanese_drivers_max |> select(phoneme,      max_contrib      = contribution),
  by = "phoneme"
) |>
  mutate(across(c(centroid_contrib, max_contrib), \(x) replace_na(x, 0))) |>
  arrange(desc(max_contrib))

# ---- 2d. Attach IPA names to the Japanese max-pair drivers ----------------
japanese_drivers_max_named <- japanese_drivers_max |>
  left_join(pnas_key, by = c("phoneme" = "phoneme_id")) |>
  select(phoneme, ipa, class, global_n, IDF, ph_prevalence, contribution, share)

print(japanese_drivers_max_named, n = 25)
print(driver_compare, n = 25)


# =============================================================================
# Section 3 — Final phoneme drivers table (cross-baseline)
# =============================================================================
# One row per phoneme: how prevalent it is across the four reference sets, which
# baselines carry it, its Bernoulli variance per set, and a leverage score +
# bidirectional flag for whether it drives the right skew of the per-language
# unrelated null and/or a baseline's observed similarity.
driver_table <- tibble(
  phoneme              = feat,
  IDF                  = idf_aligned,
  ph_prevalence        = colMeans(phoneme_mat[ph_lang,  feat, drop = FALSE]),  # 58 PH langs
  unrelated_prevalence = colMeans(phoneme_mat[unr_lang, feat, drop = FALSE]),  # 211 controls
  span_has             = phoneme_mat["Spanish",  feat] == 1,
  jap_has              = phoneme_mat["Japanese", feat] == 1,
  eng_has              = phoneme_mat["English",  feat] == 1
) |>
  # Austronesian prevalence (over 284 Austronesian langs) from the freq file
  left_join(phoneme_freq |> select(phoneme, austronesian_prevalence = freq),
            by = "phoneme") |>
  # Global prevalence (over the full ~1772-lang Ruhlen DB); NA -> 0 for ultra-rare
  left_join(phoneme_freq_global |> select(phoneme, global_prevalence = freq),
            by = "phoneme") |>
  mutate(global_prevalence = coalesce(global_prevalence, 0)) |>
  # IPA / class / global occurrence count
  left_join(pnas_key |> select(phoneme = phoneme_id, ipa, class, global_n),
            by = "phoneme") |>
  mutate(
    # Bernoulli variance p(1-p) per reference set: how much the phoneme varies
    # within that set (0 = uniform, i.e. present in ~all or ~none; max 0.25 at p=0.5)
    ph_variance           = ph_prevalence           * (1 - ph_prevalence),
    austronesian_variance = austronesian_prevalence * (1 - austronesian_prevalence),
    unrelated_variance    = unrelated_prevalence    * (1 - unrelated_prevalence),
    global_variance       = global_prevalence       * (1 - global_prevalence),
    # unrelated-vs-baseline overlap ("present in unrelated" = >=1 unrelated lang has it):
    #  - unrelated_not_baseline: in the unrelated set but in NO baseline
    #  - n_baselines_shared_unrelated: if in the unrelated set, how many baselines also
    #    carry it (0-3); 0 when the phoneme is absent from the unrelated set
    unrelated_not_baseline       = unrelated_prevalence > 0 & !(span_has | jap_has | eng_has),
    n_baselines_shared_unrelated = if_else(unrelated_prevalence > 0,
                                           as.integer(span_has + jap_has + eng_has), 0L),
    any_baseline = span_has | jap_has | eng_has,
    # ---- right-skew leverage on the per-language unrelated null --------------
    # Each PH language's null = its cosine to the 211 unrelated languages, and the
    # null is right-skewed. A shared phoneme contributes IDF^2 to the cosine
    # numerator; over the unrelated set that term is Bernoulli(q) with q =
    # unrelated_prevalence, so its (proportional) contribution to the null's 3rd
    # central moment is IDF^6 * q(1-q)(1-2q): positive (right) for q < 0.5, ZERO at
    # q = 0 (a phoneme in ~no unrelated lang can't form a tail) and at q = 0.5,
    # peaking near q ~ 0.21 (a moderate minority). Gated to phonemes present in at
    # least one PH inventory (ph_prevalence > 0), since only those can drive a PH
    # language's null. High leverage = rare-in-Austronesian (high IDF), carried by
    # some PH language, shared with a moderate minority of unrelated languages.
    null_skew_leverage = if_else(
      ph_prevalence > 0,
      IDF^6 * unrelated_prevalence * (1 - unrelated_prevalence) * (1 - 2 * unrelated_prevalence),
      0
    )
  ) |>
  # bidirectional driver flag: a rare, high-leverage phoneme points two ways --
  # toward the unrelated null's right tail (inflating the comparison distribution),
  # and/or toward a baseline (inflating the OBSERVED similarity). "both" phonemes
  # inflate the observed baseline AND the null it is tested against, so they are the
  # ones to scrutinize. The leverage cutoff (top decile of candidates) is tunable.
  mutate(
    drives_unrelated_tail = null_skew_leverage >=
      quantile(null_skew_leverage[ph_prevalence > 0 & unrelated_prevalence > 0], 0.9),
    flag = case_when(
      drives_unrelated_tail & any_baseline ~ "both (baseline + unrelated-tail)",
      drives_unrelated_tail                ~ "unrelated-tail driver",
      any_baseline                         ~ "baseline driver",
      TRUE                                 ~ "neither"
    ) |> factor(levels = c("both (baseline + unrelated-tail)", "unrelated-tail driver",
                           "baseline driver", "neither"))
  ) |>
  select(phoneme, ipa, class, global_n, IDF,
         ph_prevalence,           ph_variance,
         austronesian_prevalence, austronesian_variance,
         unrelated_prevalence,    unrelated_variance,
         global_prevalence,       global_variance,
         span_has, jap_has, eng_has,
         unrelated_not_baseline, n_baselines_shared_unrelated,
         null_skew_leverage, drives_unrelated_tail, flag) |>
  arrange(desc(null_skew_leverage))

print(driver_table, n = 30)
print(count(driver_table, flag))

write.csv(driver_table, here("data", "PHONEME_driver_table.csv"), row.names = FALSE)


# =============================================================================
# Suggested spot-check methods (not implemented — directions for deeper checks)
# =============================================================================
# The driver_table above is descriptive. To confirm which flagged phonemes truly
# drive the baseline-similarity signal, consider:
#
# 1. Point-biserial correlation — for each phoneme, correlate its presence (0/1)
#    across the 58 PH languages with cossim_span / cossim_jap / cossim_eng. Large
#    |r| = the phoneme's presence tracks the similarity signal (a real driver),
#    small |r| = it varies but doesn't move the baseline cosine.
#
# 2. Regularized regression (LASSO, glmnet) — regress a baseline's cossim on the
#    full PH phoneme presence matrix. The non-zero coefficients give a sparse,
#    signed driver set with multivariate control (handles correlated phonemes that
#    marginal correlations double-count).
#
# 3. Enrichment test (Fisher's exact) — for each phoneme, test whether its presence
#    is over-represented in the *_influenced languages from [2] vs. the rest. Ties
#    the drivers back to the population-level influence classification.
#
# 4. Prevalence contrasts — ph_prevalence - unrelated_prevalence flags phonemes that
#    are PH-specific vs. generic; comparing to austronesian_prevalence shows whether
#    a phoneme is Austronesian-characteristic or a local innovation. A phoneme shared
#    with a baseline but also common in the unrelated controls is not baseline-specific.
#
# 5. Leave-one-phoneme-out — drop a phoneme, recompute the weighted cosine, and
#    measure the change in cossim (and, if desired, in [2]'s influence / the Mantel
#    results) to gauge how sensitive the conclusions are to any single phoneme.
#
# 6. Plausibility check — inspect ipa + global_n to catch ultra-rare / coding-artifact
#    phonemes (a high-IDF phoneme with tiny global_n may be a data artifact), and judge
#    whether contact borrowing of that phoneme between the languages is realistic.
