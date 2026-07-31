# =============================================================================
# [0] Grammar Analysis — Build the GRAMBANK feature database
#
# Builds the Philippine + interest-language feature matrix from GRAMBANK using
# the reduction thresholds settled by [0]_prelim.R, hard-codes Spanish (absent
# from GRAMBANK), and derives the unrelated control set as phoneme's
# Ruhlen-derived controls that also have GRAMBANK coverage, so both analyses
# share one control set.
#
# Run order: [0]_prelim.R, PART A, then [0]_Phylogenetic_Tree.R (which turns
# `Ph_Languages` and `GRAMBANKdf_PH_maximized` into `Ph_Languages_pruned`),
# then PART B.
#
# Input:   data/languages.csv, data/values.csv, data/RUHLENdf_PH.csv
#          (phoneme's control-set source, for the dynamic unrelated set)
# Outputs: data/GRAMBANKdf_full.csv, data/gramfeature_freq.csv
# =============================================================================

library(lingtypology)
library(dplyr)
library(tidyverse)
library(readr)
library(glue)
library(here)

# ---- 1. Load GRAMBANK languages + values ----
languages <- read_csv(here("data", "languages.csv"), show_col_types = FALSE)
values <- read_csv(here("data", "values.csv"), show_col_types = FALSE)

# ---- 2. Filter to Philippine languages (by lat/long) plus English/Japanese ----
# Spanish is added by hand below (step 4) — GRAMBANK has no Spanish coverage.
Philippine_langs <- languages %>%
  filter(
    (Latitude > 4.5 & Latitude < 21 &
       Longitude > 115 & Longitude < 128) |
      Name %in% c('English', 'Japanese')
  )

grambank_values_ph <- values %>%
  filter(Language_ID %in% Philippine_langs$ID)

# ---- 3. Pivot to one-row-per-language feature matrix (treat "?" as NA) and attach metadata ----
GRAMBANKdf_PH <- grambank_values_ph %>%
  mutate(Value = na_if(Value, "?"),
         Value = as.numeric(Value)) %>%
  select(Language_ID, Parameter_ID, Value) %>%
  pivot_wider(names_from = Parameter_ID,
              values_from = Value)

GRAMBANKdf_PH <- GRAMBANKdf_PH %>%
  left_join(Philippine_langs %>% select(ID, Name, Longitude, Latitude, Family_name, Macroarea), by = c("Language_ID" = "ID"))
GRAMBANKdf_PH <- GRAMBANKdf_PH %>%
  rename(language = Name)

# Apply [0]_prelim.R's settled thresholds in one deterministic pass. NOTE: the
# loop is deliberately not run to convergence — it over-prunes to nothing by
# iteration 3, so the iter_2 snapshot is taken explicitly below.
metadata_cols <- c("Language_ID", "language", "Family_name", "Macroarea", "Longitude", "Latitude")
feature_thresh <- 80
language_thresh <- 50

for (iteration in 1:2) {
  feature_cols <- names(GRAMBANKdf_PH)[!names(GRAMBANKdf_PH) %in% metadata_cols]
  keep_features <- feature_cols[colSums(!is.na(GRAMBANKdf_PH[feature_cols])) >= feature_thresh]
  GRAMBANKdf_PH <- GRAMBANKdf_PH[, c(metadata_cols, keep_features), drop = FALSE]
  lang_filter <- rowSums(!is.na(GRAMBANKdf_PH[, keep_features, drop = FALSE])) >= language_thresh
  GRAMBANKdf_PH <- GRAMBANKdf_PH[lang_filter, , drop = FALSE]
}

GRAMBANKdf_PH_maximized <- GRAMBANKdf_PH %>%
  mutate(Language_Type = case_when(
    language %in% c("English", "Japanese") ~ "Language of Interest",
    TRUE                                   ~ "Philippine Language"
  ))

message(
  "Retained (feature_thresh = ", feature_thresh, ", language_thresh = ", language_thresh, "): ",
  length(keep_features), " features, ", nrow(GRAMBANKdf_PH_maximized), " languages."
)

Ph_Languages <- GRAMBANKdf_PH_maximized %>%
  filter(Language_Type == "Philippine Language") %>%
  pull(language)

# =============================================================================
# >>> END OF PART A <<<
# STOP here and run [0]_Phylogenetic_Tree.R: it consumes `Ph_Languages` and
# `GRAMBANKdf_PH_maximized` and returns the `Ph_Languages_pruned` PART B needs.
# =============================================================================
stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `Ph_Languages_pruned` is not defined." =
    exists("Ph_Languages_pruned")
)

# ----------------------------- PART B (run after tree) -----------------------

# Drop Philippine languages that didn't survive tree pruning, mirroring
# [0]_CREANZA_RUHLENdatabase.R PART B. Interest languages are untouched, since
# the tree only covers Philippine languages.
n_before_pruning <- sum(GRAMBANKdf_PH_maximized$Language_Type == "Philippine Language")

GRAMBANKdf_PH_maximized <- GRAMBANKdf_PH_maximized %>%
  filter(Language_Type != "Philippine Language" | language %in% Ph_Languages_pruned)

message(
  "\nDropped by tree pruning: ", n_before_pruning - sum(GRAMBANKdf_PH_maximized$Language_Type == "Philippine Language"),
  " Philippine language(s) (", n_before_pruning, " -> ",
  sum(GRAMBANKdf_PH_maximized$Language_Type == "Philippine Language"), ")."
)

# ---- 4. Spanish — hard-coded, not queried ----
# Spanish is absent from GRAMBANK, so these feature values are hand-sourced from
# WALS. NOTE: the source coordinates had Longitude/Latitude swapped; fixed here.
spanish_row <- tibble(
  Language_ID = "stan1288", language = "Spanish",
  Family_name = "Indo-European", Macroarea = "Eurasia",
  Longitude = -4, Latitude = 40,
  GB020 = 1, GB021 = 1, GB022 = 1, GB023 = 0, GB028 = 0, GB030 = 1, GB031 = 0,
  GB035 = 1, GB036 = 0, GB037 = 0, GB042 = 0, GB043 = 0, GB044 = 1, GB051 = 1,
  GB052 = 0, GB053 = 0, GB054 = 0, GB065 = 1, GB070 = 0, GB071 = 1, GB072 = 0,
  GB073 = 0, GB079 = 1, GB080 = 1, GB082 = 1, GB083 = 1, GB084 = 1, GB086 = 1,
  GB089 = 1, GB090 = 1, GB091 = 1, GB092 = 0, GB093 = 1, GB094 = 0, GB107 = 1,
  GB121 = 1, GB130 = 1, GB131 = 0, GB137 = 0, GB138 = 0, GB171 = 1, GB172 = 1,
  GB186 = 1, GB192 = 1, GB196 = 1, GB197 = 1, GB316 = 0, GB318 = 1, GB321 = 1,
  GB415 = NA_real_,
  Language_Type = "Language of Interest"
) %>%
  select(any_of(colnames(GRAMBANKdf_PH_maximized)))

GRAMBANKdf_PH_maximized <- bind_rows(GRAMBANKdf_PH_maximized, spanish_row)

# ---- 5. Query GRAMBANK for the unrelated-control search ----
GRAMBANK_query <- grambank.feature(c('gb020','gb021','gb022','gb023','gb028','gb030','gb031','gb035','gb036','gb037','gb042','gb043','gb044','gb051','gb052','gb053',
                                     'gb054','gb065','gb070','gb071','gb072','gb073','gb079','gb080','gb082','gb083','gb084','gb086','gb089','gb090','gb091','gb092',
                                     'gb093','gb094','gb107','gb121','gb130','gb131','gb137','gb138','gb171','gb172','gb186','gb192','gb196','gb197','gb316','gb318',
                                     'gb321','gb415'), na.rm = FALSE)

# ---- 6. Derive the unrelated control set — mirror phoneme exactly ----
# The control set must be present in BOTH datasets: take phoneme's unrelated set
# and keep those with GRAMBANK coverage, bridging Ruhlen ISO to GRAMBANK glottocode.
RUHLENdf_PH <- read_csv(here("data", "RUHLENdf_PH.csv"), show_col_types = FALSE)

ruhlen_unrelated <- RUHLENdf_PH %>%
  filter(Language_type == "Unrelated Language") %>%
  select(language, iso6393)

grambank_iso <- GRAMBANK_query %>%
  mutate(iso6393 = iso.gltc(glottocode))

unrelated_matched <- ruhlen_unrelated %>%
  inner_join(grambank_iso, by = "iso6393")

dropped_controls <- setdiff(ruhlen_unrelated$language, unrelated_matched$language)

message(
  "\nUnrelated control set (mirrors phoneme's Ruhlen-derived set):",
  "\n  Phoneme unrelated languages : ", nrow(ruhlen_unrelated),
  "\n  With GRAMBANK coverage      : ", nrow(unrelated_matched),
  "\n  Dropped (", length(dropped_controls), "): ",
  paste0(dropped_controls, collapse = ", ")
)

GRAMBANKdf_unrelated <- unrelated_matched %>%
  left_join(languages %>% select(ID, Family_name, Macroarea), by = c("glottocode" = "ID")) %>%
  rename(Language_ID = glottocode, Longitude = longitude, Latitude = latitude) %>%
  mutate(Language_Type = "Unrelated Language",
         across(starts_with("GB"), as.numeric)) %>%
  select(any_of(colnames(GRAMBANKdf_PH_maximized)))

# ---- 7. Combine Philippine + interest (incl. Spanish) + unrelated sets ----
GRAMBANKdf_full <- bind_rows(GRAMBANKdf_PH_maximized, GRAMBANKdf_unrelated)

write_csv(GRAMBANKdf_full, here("data", "GRAMBANKdf_full.csv"))

message(
  "\nFinal GRAMBANKdf_full: ", nrow(GRAMBANKdf_full), " languages ",
  "(", sum(GRAMBANKdf_full$Language_Type == "Philippine Language"), " Philippine, ",
  sum(GRAMBANKdf_full$Language_Type == "Language of Interest"), " interest, ",
  sum(GRAMBANKdf_full$Language_Type == "Unrelated Language"), " unrelated)."
)

# ---- 8. Per-feature frequency + inverse-document-frequency weights ----
feature_cols <- intersect(colnames(GRAMBANKdf_full), colnames(GRAMBANK_query))

GRAMBANK_freq <- GRAMBANK_query %>%
  select(all_of(feature_cols))

n_query_languages <- nrow(GRAMBANK_query)

gramfeature_freq <- GRAMBANK_freq %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature, value) %>%
  summarise(n_languages = n(), .groups = "drop") %>%
  mutate(
    freq = n_languages / n_query_languages,
    IDF = log(1 / freq)
  )

write_csv(gramfeature_freq, here("data", "gramfeature_freq.csv"))
