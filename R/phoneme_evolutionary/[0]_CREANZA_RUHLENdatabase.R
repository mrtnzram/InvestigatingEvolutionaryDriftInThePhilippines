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
# Outputs: data/RUHLENdf_PH.csv, data/phoneme_freq_ruhlen.csv
#
# RUN ORDER (this script has a dependency on the phylogenetic tree):
#   1. Run PART A below (down to the `Ph_Languages` definition).
#   2. Run [0]_Phylogenetic_Tree.R — it consumes `Ph_Languages` and produces
#      `Ph_Languages_pruned` (the tree-validated subset used downstream).
#   3. Run PART B below (everything after the banner), which needs
#      `Ph_Languages_pruned`.
# The Grambank unrelated-control set is now built locally (see PART B); no
# grammar script needs to be sourced.
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

# =============================================================================
# >>> END OF PART A <<<
# STOP here and run [0]_Phylogenetic_Tree.R, which consumes `Ph_Languages`
# (defined above) and returns `Ph_Languages_pruned` — the tree-validated subset.
# PART B below requires `Ph_Languages_pruned` to be in the environment.
# =============================================================================
stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `Ph_Languages_pruned` is not defined." =
    exists("Ph_Languages_pruned")
)

# ----------------------------- PART B (run after tree) -----------------------

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

# ---- Build GRAMBANKdf_unrelated locally (no grammar script needed) ----
# The "unrelated" baseline must be languages that are covered in BOTH the Ruhlen
# phoneme data and the Grambank grammar data, so the two analyses share a control
# set. This block reproduces the unrelated-set definition from
# [0]_GRAMBANKdatabase.R, but WITHOUT its iterative matrix-reduction loop: the
# only value that filter consumed from the reduced matrix is the macroarea(s) of
# the Philippine languages (`relatedmacroareas`), which is simply "Papunesia".
languages <- read_csv(here("data", "languages.csv"), show_col_types = FALSE)

# Macroarea(s) of the Philippine languages — the sole input the unrelated-set
# filter needs from Grambank's (here-dropped) reduction step. Resolves to
# "Papunesia"; deriving it from the bbox keeps this self-contained.
relatedmacroareas <- languages |>
  filter(Latitude > 4.5 & Latitude < 21 & Longitude > 115 & Longitude < 128) |>
  pull(Macroarea) |>
  unique()

# Query the same 50 Grambank features used by the grammar analysis.
GRAMBANK_query <- grambank.feature(
  c('gb020','gb021','gb022','gb023','gb028','gb030','gb031','gb035','gb036','gb037','gb042','gb043','gb044','gb051','gb052','gb053',
    'gb054','gb065','gb070','gb071','gb072','gb073','gb079','gb080','gb082','gb083','gb084','gb086','gb089','gb090','gb091','gb092',
    'gb093','gb094','gb107','gb121','gb130','gb131','gb137','gb138','gb171','gb172','gb186','gb192','gb196','gb197','gb316','gb318',
    'gb321','gb415'),
  na.rm = FALSE)

# Keep languages outside the Austronesian family and outside the Philippine
# macroarea(s) — i.e. genealogically and areally unrelated controls.
GRAMBANKdf_unrelated <- GRAMBANK_query |>
  left_join(languages |> dplyr::select(ID, Family_name, Macroarea), by = c("glottocode" = "ID")) |>
  filter(!Family_name %in% "Austronesian",
         !Macroarea   %in% relatedmacroareas) |>
  dplyr::select(glottocode)

# Bridge Grambank glottocodes to ISO via lingtypology::languages
grambank_iso <- GRAMBANKdf_unrelated |>
  mutate(iso6393 = iso.gltc(glottocode))

# Inner join on ISO — only keep ruhlen languages with full Grambank coverage
Unrelated_Langauges <- ruhlen_candidates |>
  inner_join(
    grambank_iso |> dplyr::select(glottocode, iso6393),
    by = join_by(iso6393 == iso6393)
  ) |>
  pull(language)

message("Unrelated languages (n = ", length(Unrelated_Langauges), "):\n",
        paste0("  ", Unrelated_Langauges, collapse = "\n"))

# ---- Exclusions from the 211-language unrelated baseline -------------------
# NOTE: names must match colnames(null_mat) EXACTLY, RTF artifacts included.
# Verify with: setdiff(contact_contaminated, colnames(null_mat))  # must be empty

# Tier A: Indo-European — related to Spanish and/or English by descent.
# This is the primary contaminant. Not contact; phylogeny.
ie_related <- c(
  # Romance — Spanish's own subfamily
  "Latin", "Italian", "Romansch", "French", "Catalan", "Galician", "Portuguese",
  # Germanic — English's own subfamily (Frisian = closest living relative)
  "Dutch", "Frisian", "Faroese",
  # Wider Indo-European
  "Kashmiri", "Marathi", "Konkani", "Punjabi", "Bhojpuri", "Maithili", "Pashto",
  "Albanian", "Classical Greek", "Greek", "Irish", "Breton",
  "Latvian", "Lithuanian", "Russian", "Byelorussian", "Polish", "Czech",
  "Slovene", "Macedonian", "Armenian"
)

# Tier B: Spanish colonial contact with documented inventory-level borrowing
spanish_contact <- c(
  # Mesoamerica (Bennett 2016; Suarez 1983)
  "Huastec", "Yucatec", "Tojolabal", "Chuj", "Acatec", "Quiche",
  "Tlamelula", "Amuzgo", "Huichol", "Tarascan",
  # South America (Adelaar & Muysken 2004; Smeets 2008)
  "Warao", "Itonama", "Mapudungu", "Qawasqar", "Cubeo", "Huambisa",
  "Cocama", "Murui", "Chimane", "Chulupi", "Matses", "Eseejja", "Huarayo",
  # California / SW missions (Bright 1960; Bright & Bright 1959; Bright 1965)
  "Luisen~o", "Cahuilla", "Pima", "Patwin", "Wintun",
  "Central Sierra Miwok", "Southern Sierra Miwok", "Maricopa",
  # Spanish + English, Southwest / Plains (Campbell 1997)
  "Navajo", "Western Apache", "Southern Paiute", "Chemehuevi",
  "Kiowa", "Tonkawa"
)

# Tier B2: Ibero-Romance (Portuguese) contact — same phoneme signature as Spanish
portuguese_contact <- c(
  "Nhengatu",       # also a Tupi-based contact language; excluded on both grounds
  "Pakaasnovos", "Bororo", "Mashakali",
  "Makua", "Lwena"
)

contact_contaminated <- c(ie_related, spanish_contact, portuguese_contact)

unrelated_clean <- setdiff(Unrelated_Langauges, contact_contaminated)
length(unrelated_clean)

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
  dplyr::select(iso6393, language, source, latitude, longitude, all_of(phoneme_cols), Language_type)

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
# Denominator: all Austronesian Languages excluding Philippine Languages 
# to avoid circularity

compute_phoneme_freq <- function(data, n_total = nrow(data)) {
  data |>
    summarise(across(all_of(phoneme_cols), \(x) sum(x, na.rm = TRUE))) |>
    pivot_longer(
      cols      = everything(),
      names_to  = "phoneme",
      values_to = "n_languages"
    ) |>
    mutate(
      n_total = n_total,
      freq    = n_languages / n_total,
      IDF     = log((n_total + 1) / (n_languages + 1))  # Laplace-smoothed, matches your Methods spec
    )
}

ruhlen_austronesian <- ruhlen_raw |>
  filter(
    language_family == "Austronesian",
    !language %in% Philippine_langs
  )

n_total_languages_an <- nrow(ruhlen_austronesian)

phoneme_freq_austronesian <- compute_phoneme_freq(ruhlen_austronesian, n_total_languages_an)

write_csv(phoneme_freq_austronesian, here("data", "phoneme_freq_ruhlen_austronesian.csv"))

message(
  "\nDone.",
  "\n  ruhlen_austronesian : ", n_total_languages_an, " languages × ",
  ncol(ruhlen_austronesian) - 4L, " phoneme features",
  "\n  phoneme_freq_austronesian : ", nrow(phoneme_freq_austronesian), " attested phonemes"
)

# Phoneme names key mapping -------------------------------------

# ---- build the phoneme_n -> IPA keymap -------------------------------------
raw   <- read_lines(here('data','pnas_1424033112_sd02.txt'))
hdr_i <- which(str_starts(raw, "Column\t"))

pnas_key <- read_tsv(I(raw[hdr_i:length(raw)]), show_col_types = FALSE) |>
  transmute(
    phoneme_id = str_c("phoneme_", str_pad(Column - 9L, 3, pad = "0")),  # SI col 10 -> phoneme_001
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

# ---- sanity check: keymap should match your actual phoneme_cols -----------
missing_from_key <- setdiff(phoneme_cols, pnas_key$phoneme_id)
if (length(missing_from_key) > 0)
  warning(length(missing_from_key), " phoneme_cols not found in pnas_key: ",
          paste(head(missing_from_key, 10), collapse = ", "))

pnas_key
