# =============================================================================
# [0] Phoneme Analysis — Build the Ruhlen (Creanza et al. 2015 PNAS) feature database
#
# Reads binary phoneme inventories from the Creanza et al. 2015 supplemental
# dataset (pnas_1424033112_sd01.txt), filters to the study languages, and
# writes the feature matrix and IDF frequency table in a schema compatible with
# the original PHOIBLE pipeline.
#
# Source: Creanza et al. 2015 PNAS, SI Dataset 1 — 2082 languages × 728 phonemes
# Phoneme column names (phoneme_001–phoneme_728) correspond to the positional
# indices in SI Dataset 2 (not included here; join on position if IPA names needed).
#
# Outputs: data/PHOIBLEdf_PH.csv, data/phoneme_freq.csv
# =============================================================================

library(lingtypology)
library(tidyverse)
library(readr)
library(here)

# ---- 1. Define column names for the 737-column Ruhlen file ----
# Columns 1–9 are metadata; columns 10–737 are 728 binary phoneme features.
phoneme_cols <- paste0("phoneme_", str_pad(1:728, 3, pad = "0"))

col_names_ruhlen <- c(
  "record_num", "language", "iso6393", "iso_a3", "language_family",
  "population", "region", "latitude", "longitude",
  phoneme_cols
)

# ---- 2. Load Ruhlen data (header = 19 lines, data starts line 20) ----
lines <- readLines(here("data", "pnas_1424033112_sd01.txt"))

ruhlen_raw <- read_tsv(
  I(paste(lines[20:length(lines)], collapse = "\n")),
  col_names      = col_names_ruhlen,
  show_col_types = FALSE
) |>
  mutate(across(all_of(phoneme_cols), \(x) replace_na(as.integer(x), 0L)))

message("Loaded Ruhlen data: ", nrow(ruhlen_raw), " languages × ",
        length(phoneme_cols), " phoneme features")

# filter duplicate English and Spanish

ruhlen_raw <- ruhlen_raw |>
  filter(
    !(language == "Spanish" & region != "Europe"),
    !(language == "English" & region != "AmericasNorthCentral")
  ) |>
  slice_max(population, by = language, n = 1, with_ties = FALSE)

# ---- 3. Define the study language set ----
# Philippine languages: derived from the Ruhlen lat/lon coordinates.
# Bounding box covers the Philippine archipelago (no hardcoded name list needed).
Ph_Languages <- ruhlen_raw |>
  filter(
    latitude  >  4.5 & latitude  < 21 &
      longitude > 115  & longitude < 128
  ) |>
  pull(language)

message("Philippine languages detected by bounding box (n = ", length(Ph_Languages), "):\n",
        paste0("  ", Ph_Languages, collapse = "\n"))

Interest_Languages <- c("English", "Spanish", "Japanese")

# --- Derive unrelated languages set ----------

# Get Philippine families and regions directly from ruhlen
ph_families <- ruhlen_raw |>
  filter(language %in% ph_lang) |>
  pull(language_family) |>
  unique()

ph_regions <- ruhlen_raw |>
  filter(language %in% ph_lang) |>
  pull(region) |>
  unique()

# Ruhlen candidates — same region, different family, not already in study set
ruhlen_candidates <- ruhlen_raw |>
  filter(
    region           %in% ph_regions,
    !language_family %in% ph_families,
    !language        %in% Languages
  )

# Bridge Grambank glottocodes to ISO via lingtypology::languages
grambank_iso <- GRAMBANKdf_unrelated |>
  mutate(iso6393 = iso.gltc(glottocode))

# Inner join on ISO — only keep ruhlen languages with full Grambank coverage
Unrelated_Langauges <- ruhlen_candidates |>
  inner_join(
    grambank_iso |> select(glottocode, iso6393),
    by = join_by(iso6393 == iso6393)
  ) |>
  pull(language)

message("Unrelated languages (n = ", length(unr_lang), "):\n",
        paste0("  ", unr_lang, collapse = "\n"))


# -- All Languages

Languages <- c(Ph_Languages, Interest_Languages, Unrelated_Langauges)

# ---- 4. Filter to study languages and build output schema ----
# No source deduplication needed: Ruhlen has exactly one record per language.
# We add source = "ruhlen" to preserve schema compatibility with the PHOIBLE output.
PHOIBLEdf_PH <- ruhlen_raw |>
  filter(language %in% Languages) |>
  mutate(
    source        = "ruhlen",
    Language_type = case_when(
      language %in% Ph_Languages       ~ "Philippine Language",
      language %in% Interest_Languages ~ "Language of Interest",
      language %in% Unrelated_Langauges     ~ "Unrelated Language"
    )
  ) |>
  select(iso6393, language, source, latitude, longitude, all_of(phoneme_cols), Language_type)

# Diagnostic: report any study languages absent from the Ruhlen database
unmatched <- setdiff(Languages, PHOIBLEdf_PH$language)
if (length(unmatched) > 0) {
  message(
    "\nWarning: ", length(unmatched),
    " study language(s) not found in the Ruhlen database:\n",
    paste0("  - ", unmatched, collapse = "\n"),
    "\nCheck spelling against ruhlen_raw$language or consider alternate names."
  )
} else {
  message("All ", length(Languages), " study languages matched successfully.")
}

# ---- 5. Preview geographic spread ----
map.feature(PHOIBLEdf_PH$language, PHOIBLEdf_PH$Language_type)

write_csv(PHOIBLEdf_PH, here("data", "RUHLENdf_PH"))

# ---- 6. Compute global phoneme frequencies and IDF weights ----
# Denominator: all 2082 Ruhlen languages (one row per language, no deduplication).
# This parallels the PHOIBLE pipeline's PHOIBLEdf_clean step.
n_total_languages <- nrow(ruhlen_raw)

phoneme_freq <- ruhlen_raw |>
  summarise(across(all_of(phoneme_cols), \(x) sum(x, na.rm = TRUE))) |>
  pivot_longer(
    cols      = everything(),
    names_to  = "phoneme",
    values_to = "n_languages"
  ) |>
  filter(n_languages > 0) |>
  mutate(
    freq = n_languages / n_total_languages,
    IDF  = log(1 / freq)
  )

write_csv(phoneme_freq, here("data", "phoneme_freq_ruhlen.csv"))

message(
  "\nDone.",
  "\n  PHOIBLEdf_PH : ", nrow(PHOIBLEdf_PH), " languages × ",
  ncol(PHOIBLEdf_PH) - 4L, " phoneme features",  # minus iso6393, language, source, Language_type
  "\n  phoneme_freq : ", nrow(phoneme_freq), " attested phonemes"
)
