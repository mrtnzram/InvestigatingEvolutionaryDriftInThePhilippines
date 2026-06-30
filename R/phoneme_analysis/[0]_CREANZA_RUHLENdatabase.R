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

# filter further after running [0]_Phylogenetic_Tree.R use Ph_Languages_pruned downstream

Interest_Languages <- c("English", "Spanish", "Japanese")

# --- Derive unrelated languages set ----------

# Get Philippine families and regions directly from ruhlen
ph_families <- ruhlen_raw |>
  filter(language %in% Ph_Languages_pruned) |>
  pull(language_family) |>
  unique()

ph_regions <- ruhlen_raw |>
  filter(language %in% Ph_Languages_pruned) |>
  pull(region) |>
  unique()

# Bounding box buffer (km) — controls for areal contact effects, following
# the convention in typological sampling of excluding/modeling languages
# within ~1000 km to avoid areal (contact-induced) similarity confounding
# genealogical signal. See Jaeger, Graff, Croft & Pontillo (2011),
# "Mixed effect models for genetic and areal dependencies in linguistic
# typology," Linguistic Typology 15(2): 281-319.
buffer_km <- 1000

ph_coords <- ruhlen_raw |>
  filter(language %in% Ph_Languages)

ph_bbox <- list(
  lat_min = min(ph_coords$latitude, na.rm = TRUE),
  lat_max = max(ph_coords$latitude, na.rm = TRUE),
  lon_min = min(ph_coords$longitude, na.rm = TRUE),
  lon_max = max(ph_coords$longitude, na.rm = TRUE)
)

ph_centroid_lat <- mean(c(ph_bbox$lat_min, ph_bbox$lat_max))

# 1 deg latitude ~= 111.32 km everywhere; 1 deg longitude shrinks with cos(lat)
# 111.32 conversion rate
deg_lat_buffer <- buffer_km / 111.32
deg_lon_buffer <- buffer_km / (111.32 * cos(ph_centroid_lat * pi / 180))

ph_bbox_buffered <- list(
  lat_min = ph_bbox$lat_min - deg_lat_buffer,
  lat_max = ph_bbox$lat_max + deg_lat_buffer,
  lon_min = ph_bbox$lon_min - deg_lon_buffer,
  lon_max = ph_bbox$lon_max + deg_lon_buffer
)

# Ruhlen candidates — same region, different family, not already in study
# set, and outside the contact-effect buffer zone around the Philippines
ruhlen_candidates <- ruhlen_raw |>
  filter(
    !region           %in% ph_regions,
    !language_family %in% ph_families,
    !language        %in% c(Ph_Languages,Interest_Languages),
    !(
      latitude  >= ph_bbox_buffered$lat_min & latitude  <= ph_bbox_buffered$lat_max &
        longitude >= ph_bbox_buffered$lon_min & longitude <= ph_bbox_buffered$lon_max
    )
  )

# insert GRAMBANKdf_unrelated here

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

message("Unrelated languages (n = ", length(Unrelated_Langauges), "):\n",
        paste0("  ", Unrelated_Langauges, collapse = "\n"))

# Languages with independent contact to the contrast set (Spanish/English/Japanese),
# which would inject the very signal being measured into the "unrelated" baseline.
# Tiers reflect strength + structurality of contact (see methods notes).
contact_contaminated <- c(
  "Morr",   # Japanese colonial (1910-45) + American English (post-1945) — both axes
  "Ainu",     # Japanese: lexical + grammatical (analytic constructions)
  "Mandarin", # Japanese (wasei-kango back-borrowing) + English lexical
  "Wu",       # as Mandarin (treaty-port Shanghai)
  "Nivkh",    # "Nivkh"(Gilyak): Russian-dominated; Japanese contact thin
  "Burmese"   # "Burmese": British (not American) English lexical contact
)


unrelated_clean <- setdiff(Unrelated_Langauges, contact_contaminated)

# -- All Languages

Languages <- c(Ph_Languages_pruned, Interest_Languages, unrelated_clean)

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

write_csv(PHOIBLEdf_PH, here("data", "RUHLENdf_PH.csv"))

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
