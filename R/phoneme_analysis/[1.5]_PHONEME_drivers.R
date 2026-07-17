# =============================================================================
# [1.5] Phoneme drivers
# Diagnostic: decomposes the IDF-weighted cosine similarity into per-phoneme
# contributions to see which phonemes drive a Philippine language's similarity to
# the interest-language baselines, and which phonemes carry meaningful vs.
# unmeaningful variance across the Philippine languages.
#
# Sections:
#   0b. PHOIBLE Austronesian IDF (from scratch via lingtypology) + Americanist->IPA crosswalk
#   1.  Cross-reference phoneme names (phoneme_id -> IPA / class, from the PNAS key)
#   2.  Japanese phoneme drivers (centroid + max-pair decomposition, worked example)
#   3.  Final phoneme drivers table — prevalences + baseline indicators + variance;
#       right-skew leverage + bidirectional flag (data/PHONEME_driver_table.csv).
#   4.  Ruhlen vs PHOIBLE IDF comparison — IDF_ruhlen vs IDF_phoible (both [0,1]),
#       idf_shift + variance_shift_dir, to see how Ruhlen's phoneme-flag adjustments
#       reweight phonemes. PHOIBLE's individual modified segments are scored against
#       the Ruhlen modification-type flag they belong to (data/PHONEME_driver_phoible.csv).
#
# Inputs:  data/RUHLENdf_PH.csv                    (binary phoneme matrix + Language_type)
#          data/phoneme_freq_ruhlen_austronesian.csv (Ruhlen IDF + Austronesian freq)
#          data/phoneme_freq_ruhlen.csv            (global Ruhlen prevalence)
#          data/PHONEME_cossim.csv                 (per-PH-language baseline cosines)
#          data/pnas_1424033112_sd0{1,2}.txt       (Ruhlen matrix + phoneme_id -> IPA key)
#          PHOIBLE via lingtypology (cached to data/phoible_freq_austronesian.csv, n=88 Austronesian)
# Outputs: data/PHONEME_driver_table.csv  (Section 3, Ruhlen-only, unchanged schema)
#          data/PHONEME_driver_phoible.csv (Section 4, cross-database IDF comparison)
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
# Section 0b — PHOIBLE Austronesian IDF + Americanist->IPA crosswalk
# =============================================================================
# IDF_phoible is computed over PHOIBLE's own Austronesian languages (glottolog
# family) minus the Philippine bounding box (n = 88), with the SAME Laplace IDF
# formula as IDF_ruhlen. Ruhlen uses Americanist symbols (č, š, ñ, dot-under
# retroflex) while PHOIBLE uses strict IPA, so a crosswalk (typology-descr.pdf
# Table 1 / §7 / §19) aligns the BASE phonemes. Ruhlen collapses each modification
# into ONE indicator while PHOIBLE encodes each modified segment individually —
# that asymmetry is handled downstream (Section 3) by matching on modification type.
#
# The PHOIBLE pull uses lingtypology (network) and is cached; delete the CSV to rebuild.
phoible_cache <- here("data", "phoible_freq_austronesian.csv")
if (!file.exists(phoible_cache)) {
  library(lingtypology)
  ph   <- phoible.feature(source = "all", na.rm = TRUE)
  inv1 <- ph |> distinct(glottocode, inventoryid) |> group_by(glottocode) |> slice(1) |> ungroup()
  ph1  <- ph |> semi_join(inv1, by = c("glottocode", "inventoryid")) |>
    filter(is.na(marginal) | !marginal)                       # one inventory/lang, drop marginal
  # global occurrence count over ALL PHOIBLE languages (parallels Ruhlen's global_n)
  global_ph <- ph1 |> distinct(glottocode, phoneme) |> count(phoneme, name = "global_n_phoible")
  # PHOIBLE Austronesian languages (glottolog family) minus the Philippine bounding box
  g <- lingtypology::glottolog |> mutate(is_aus = str_starts(replace_na(affiliation, ""), "Austronesian"))
  aus_glc <- tibble(glottocode = unique(ph1$glottocode)) |>
    left_join(g |> select(glottocode, is_aus, latitude, longitude), by = "glottocode") |>
    filter(is_aus, !(latitude > 4.5 & latitude < 21 & longitude > 115 & longitude < 128)) |>
    pull(glottocode)
  seg <- ph1 |> filter(glottocode %in% aus_glc)
  N   <- n_distinct(seg$glottocode)
  attested <- seg |> distinct(glottocode, phoneme) |> count(phoneme, name = "n_languages") |>
    mutate(n_total = N, IDF_phoible = log((N + 1) / (n_languages + 1)))
  phoible_freq <- bind_rows(attested,
      tibble(phoneme = setdiff(unique(ph1$phoneme), attested$phoneme),
             n_languages = 0L, n_total = N, IDF_phoible = NA_real_)) |>
    left_join(global_ph, by = "phoneme")
  write_csv(phoible_freq, phoible_cache)
}
phoible_freq <- read_csv(phoible_cache, show_col_types = FALSE)

# IPA normalizer (NFC + strip U+0361 tie bar, which PHOIBLE does not use) and the
# Ruhlen(Americanist) -> IPA crosswalk, applied to base classes only.
ipa_norm <- function(x) stringr::str_remove_all(stringi::stri_trans_nfc(x), "͡")
ruhlen_ipa_xwalk <- c("š"="ʃ", "ž"="ʒ", "č"="tʃ", "ǰ"="dʒ", "ǯ"="dʒ", "ĵ"="dʒ", "ñ"="ɲ",
                      "ṭ"="ʈ", "ḍ"="ɖ", "ṇ"="ɳ", "ṣ"="ʂ", "ẓ"="ʐ", "ḷ"="ɭ", "ṛ"="ɽ",
                      "ł"="ɬ", "ƛ"="tɬ")
ruhlen_to_ipa <- function(ipa, class) {
  base <- class %in% c("consonant", "vowel", "click")   # leave mod indicators unmapped
  out  <- stringi::stri_trans_nfc(ipa)                  # compose first so xwalk keys match
  for (k in names(ruhlen_ipa_xwalk))
    out[base] <- stringr::str_replace_all(out[base], stringr::fixed(k), ruhlen_ipa_xwalk[[k]])
  ipa_norm(out)
}
# normalized-IPA -> PHOIBLE IDF lookup (attested segment wins over unattested dupes)
phoible_lookup <- phoible_freq |>
  mutate(ipa_key = ipa_norm(phoneme)) |>
  arrange(desc(n_languages)) |>
  distinct(ipa_key, .keep_all = TRUE) |>
  transmute(ipa_key, IDF_phoible, n_phoible = n_languages, global_n_phoible)

# --- modification handling (for matching PHOIBLE's individual modified segments to
# --- Ruhlen's collapsed modification-type indicators) -------------------------
# Each modification mark -> a canonical modification name (Ruhlen typology §7 / IPA).

mod_map <- c("ʲ"="palatalized","ʷ"="labialized","ʰ"="aspirated","ː"="long","ˑ"="long",
             "̃"="nasalized","ˀ"="glottalized","ʔ"="glottalized","ˠ"="velarized","ˤ"="pharyngealized",
             "ⁿ"="prenasalized","ᵐ"="prenasalized","ᵑ"="prenasalized","ᴺ"="prenasalized",
             "̩"="syllabic","̪"="dental","̥"="voiceless","̤"="breathy","̰"="creaky","̬"="fortis")

# canonical (sorted, +-joined) set of modifications carried by a segment; NA if none
mod_signature <- function(x) {
  x <- stringi::stri_trans_nfc(x)
  vapply(x, function(s) {
    hit <- unique(mod_map[vapply(names(mod_map), \(m) stringr::str_detect(s, stringr::fixed(m)), logical(1))])
    if (length(hit) == 0) NA_character_ else paste(sort(hit), collapse = "+")
  }, character(1), USE.NAMES = FALSE)
}
# is a segment vowel-based? strip modifications/diacritics, test the base letter
ipa_vowels <- strsplit("aeiouyɨʉɯɪʏʊøɘɵɤəɛœɜɞʌɔæɐɑɒɶ", "")[[1]]
is_vowel_seg <- function(x) substr(stringr::str_remove_all(x, "[\\p{M}\\p{Lm}\\p{Sk}]"), 1, 1) %in% ipa_vowels
# min-max normalize to [0, 1]
mm01 <- function(x) { r <- range(x, na.rm = TRUE); (x - r[1]) / (r[2] - r[1]) }


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
# One row per Ruhlen phoneme: how prevalent it is across the four reference sets,
# which baselines carry it, its Bernoulli variance per set, and a leverage score +
# bidirectional flag for whether it drives the right skew of the per-language
# unrelated null and/or a baseline's observed similarity. (The Ruhlen-vs-PHOIBLE IDF
# comparison is a SEPARATE table, built in Section 4 below.)
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
# Section 4 — Ruhlen vs PHOIBLE IDF comparison (separate table)
# =============================================================================
# Compares each phoneme's IDF under Ruhlen vs PHOIBLE (both min-max normalized to
# [0,1]) to see how Ruhlen's phoneme-flag adjustments reweight phonemes. Ruhlen
# collapses each modification into ONE indicator; PHOIBLE encodes each modified
# segment individually. So for a PHOIBLE individual modified segment, idf_shift is
# scored against the Ruhlen modification-type flag it belongs to (matched on
# modification signature), i.e. IDF_ruhlen(mod flag) - IDF_phoible(segment). Filtering
# class == "mod_consonant" therefore shows IDF_ruhlen null + the individual PHOIBLE
# values, with idf_shift carrying the collapsed Ruhlen flag's weight.
#
# NOTE: this table does NOT feed [1]'s cosine (which is computed separately from the
# raw Ruhlen austronesian IDF); the [0,1] normalization here is display-only.

# fixed [0,1] scales (all Ruhlen phonemes / all attested PHOIBLE segments)
ruhlen_idf_range  <- range(idf_aligned, na.rm = TRUE)
phoible_idf_range <- range(phoible_freq$IDF_phoible, na.rm = TRUE)
norm_r <- function(x) (x - ruhlen_idf_range[1]) / diff(ruhlen_idf_range)
norm_p <- function(x) (x - phoible_idf_range[1]) / diff(phoible_idf_range)

# Ruhlen modification-type reference IDF, keyed by (consonant/vowel, modification sig)
ruhlen_mod_ref <- pnas_key |>
  filter(phoneme_id %in% feat, class %in% c("mod_consonant", "mod_vowel")) |>
  transmute(class_cv = if_else(class == "mod_vowel", "v", "c"),
            sig       = mod_signature(ipa),
            idf_raw   = idf_aligned[phoneme_id]) |>
  filter(!is.na(sig)) |>
  summarise(IDF_ruhlen_ref = norm_r(mean(idf_raw, na.rm = TRUE)), .by = c(class_cv, sig))

# Ruhlen BASE phonemes (mod-flag placeholders excluded — they feed idf_shift instead)
ruhlen_base <- pnas_key |>
  filter(phoneme_id %in% feat, !class %in% c("mod_consonant", "mod_vowel")) |>
  transmute(phoneme = phoneme_id, ipa, class, global_n_ruhlen = global_n,
            ruhlen_ipa = ruhlen_to_ipa(ipa, class),
            IDF_ruhlen_raw = idf_aligned[phoneme_id]) |>
  left_join(phoible_lookup, by = c("ruhlen_ipa" = "ipa_key")) |>
  rename(IDF_phoible_raw = IDF_phoible) |>
  mutate(in_ruhlen = TRUE, in_phoible = !is.na(n_phoible))

# PHOIBLE-only segments (attested in the 88 Austronesian langs, unmatched to a base
# Ruhlen phoneme). Individual modified segments are classed here so a mod_consonant
# filter surfaces them.
phoible_only <- phoible_freq |>
  filter(n_languages > 0) |>
  mutate(ruhlen_ipa = ipa_norm(phoneme)) |>
  filter(!ruhlen_ipa %in% ruhlen_base$ruhlen_ipa[ruhlen_base$in_phoible]) |>
  distinct(ruhlen_ipa, .keep_all = TRUE) |>
  transmute(ipa = phoneme, ruhlen_ipa, IDF_phoible_raw = IDF_phoible, global_n_phoible,
            in_ruhlen = FALSE, in_phoible = TRUE,
            sig  = mod_signature(phoneme),
            is_v = is_vowel_seg(phoneme),
            class = case_when(!is.na(sig) &  is_v ~ "mod_vowel",
                              !is.na(sig) & !is_v ~ "mod_consonant",
                              is_v                ~ "vowel",
                              TRUE                ~ "consonant"))

driver_phoible <- bind_rows(ruhlen_base, phoible_only) |>
  mutate(
    IDF_ruhlen  = norm_r(IDF_ruhlen_raw),
    IDF_phoible = norm_p(IDF_phoible_raw),
    class_cv    = if_else(class %in% c("mod_vowel", "vowel"), "v", "c"),
    # canonical modification carried by the segment (NA for plain phonemes) — a
    # grouping key for inspecting the variance shift per consolidated modification
    modification_type = mod_signature(ipa)
  ) |>
  left_join(ruhlen_mod_ref, by = c("class_cv", "sig")) |>
  mutate(
    idf_shift = case_when(
      in_ruhlen & !is.na(IDF_phoible)                          ~ IDF_ruhlen - IDF_phoible,  # matched base
      !in_ruhlen & class %in% c("mod_consonant", "mod_vowel") &
        !is.na(IDF_ruhlen_ref)                                 ~ IDF_ruhlen_ref - IDF_phoible,  # PHOIBLE mod vs Ruhlen flag
      TRUE                                                     ~ NA_real_
    ),
    variance_shift_dir = case_when(is.na(idf_shift) ~ NA_character_,
                                   idf_shift > 0 ~ "+", idf_shift < 0 ~ "-", TRUE ~ "0")
  ) |>
  select(phoneme, ruhlen_ipa, ipa, class, modification_type, global_n_ruhlen, global_n_phoible,
         IDF_ruhlen, IDF_phoible, idf_shift, variance_shift_dir, in_ruhlen, in_phoible) |>
  arrange(class, modification_type, desc(idf_shift))

print(driver_phoible, n = 30)
cat(sprintf("PHOIBLE comparison — both IDFs: %d | Ruhlen-only: %d | PHOIBLE-only: %d | mod-shift scored: %d\n",
            sum(!is.na(driver_phoible$IDF_ruhlen) & !is.na(driver_phoible$IDF_phoible)),
            sum(!is.na(driver_phoible$IDF_ruhlen) &  is.na(driver_phoible$IDF_phoible)),
            sum( is.na(driver_phoible$IDF_ruhlen) & !is.na(driver_phoible$IDF_phoible)),
            sum( is.na(driver_phoible$IDF_ruhlen) & !is.na(driver_phoible$idf_shift))))

write.csv(driver_phoible, here("data", "PHONEME_driver_phoible.csv"), row.names = FALSE)
