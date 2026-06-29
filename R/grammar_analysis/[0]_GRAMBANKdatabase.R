# =============================================================================
# [0] Grammar Analysis — Build the GRAMBANK feature database
# Loads GRAMBANK, filters to Philippine + control languages, pivots to a feature
# matrix, applies iterative matrix reduction, adds unrelated controls, and writes
# the per-feature IDF frequency table.
# Inputs:  data/languages.csv, data/values.csv, data/GRAMBANKdf_PH.csv (curated)
# Outputs: data/gramfeature_freq.csv
# =============================================================================

library(lingtypology)
library(dplyr)
library(tidyverse)
library(readr)
library(tidyselect)
library(proxy)
library(glue)
library(here)

# ---- 1. Load GRAMBANK languages + values ----
languages <- read_csv(here("data", "languages.csv"))
values <- read_csv(here("data", "values.csv"))

summary(languages)
summary(values)

UnrelatedLangs <- c("Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Liv",
              "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")

# ---- 2. Filter to Philippine languages (by lat/long) plus the comparison set ----
Philippine_langs <- languages %>%
  filter(
    (Latitude>4.5 & Latitude <21 &
       Longitude > 115 & Longitude <128) |
      Name %in% c('English',"Japanese",UnrelatedLangs)
  )

grambank_values_ph <- values %>%
  filter(Language_ID %in% Philippine_langs$ID)

# ---- 3. Pivot to one-row-per-language feature matrix (treat "?" as NA) and attach metadata ----
GRAMBANKdf_PH <- grambank_values_ph %>%
  mutate(Value = na_if(Value,"?"),
         Value = as.numeric(Value)) %>%
  select(Language_ID,Parameter_ID,Value) %>%
  pivot_wider(names_from = Parameter_ID,
              values_from = Value)

GRAMBANKdf_PH <- GRAMBANKdf_PH %>%
  left_join(Philippine_langs %>% select(ID, Name, Longitude, Latitude, Family_name, Macroarea), by = c("Language_ID" = "ID"))
GRAMBANKdf_PH <- GRAMBANKdf_PH %>%
  rename(language = Name)

# ---- 4. Inspect feature/language coverage to choose reduction thresholds ----
feature_counts <- sort(colSums(!is.na(GRAMBANKdf_PH[ , !(names(GRAMBANKdf_PH) %in% c("Language_ID", "language"))])),decreasing = TRUE)

plot(feature_counts,
     type = "b",
     pch = 19,
     col = "#3182bd",
     main = "Feature Coverage (Sorted)",
     xlab = "Features (ranked)",
     ylab = "Number of Languages with Data")

language_counts <- sort(rowSums(!is.na(GRAMBANKdf_PH[ , !(names(GRAMBANKdf_PH) %in% c("Language_ID", "language"))])),decreasing = TRUE)

plot(language_counts,
     type = "b",
     pch = 19,
     col = "#31a354",
     main = "Language Coverage (Sorted)",
     xlab = "Languages (ranked)",
     ylab = "Number of Features with Data")

# ---- 5. Iterative matrix reduction: prune sparse features/languages until stable ----
GRAMBANKdf_PH_forlooping <- GRAMBANKdf_PH
GRAMBANK_snapshots <- list()
metadata_cols <- c("Language_ID", "language", "Family_name", "Macroarea", "Longitude", "Latitude")

feature_thresh <- 80   # minimum number of languages per feature
language_thresh <- 50  # minimum number of features per language
iteration <- 1

repeat {
  cat(glue("\n--- Iteration {iteration} ---\n"))

  old_dim <- dim(GRAMBANKdf_PH_forlooping)

  feature_cols <- names(GRAMBANKdf_PH)[!names(GRAMBANKdf_PH) %in% c("Language_ID", "language","Family_name","Macroarea","Longitude",'Latitude')]

  # Keep features present in >= feature_thresh languages
  keep_features <- feature_cols[colSums(!is.na(GRAMBANKdf_PH[feature_cols])) >= feature_thresh]
  GRAMBANKdf_PH <- GRAMBANKdf_PH[, c(metadata_cols, keep_features), drop = FALSE]
  cat(glue("Retained {length(keep_features)} features\n"))

  # Keep languages with >= language_thresh non-NA features
  lang_filter <- rowSums(!is.na(GRAMBANKdf_PH[, keep_features, drop = FALSE])) >= language_thresh
  GRAMBANKdf_PH <- GRAMBANKdf_PH[lang_filter, , drop = FALSE]
  cat(glue("Retained {nrow(GRAMBANKdf_PH)} languages\n"))

  GRAMBANK_snapshots[[paste0("iter_", iteration)]] <- GRAMBANKdf_PH

  # Stop once the shape stabilizes
  if (all(dim(GRAMBANKdf_PH) == old_dim)) {
    cat(glue("Matrix stabilized after iteration {iteration}.\n"))
    break
  }

  # Stop if pruning becomes too aggressive
  if (nrow(GRAMBANKdf_PH) < 2 || length(keep_features) < 2) {
    warning("Matrix reduced to nearly nothing — adjust your thresholds.")
    break
  }

  iteration <- iteration + 1
}

# Inspect NA counts per snapshot and select the most stable iteration
colSums(is.na(GRAMBANK_snapshots$iter_1))
colSums(is.na(GRAMBANK_snapshots$iter_2))
colSums(is.na(GRAMBANK_snapshots$iter_200))

GRAMBANKdf_PH_maximized <- GRAMBANK_snapshots$iter_2

# ---- 6. Find unrelated control languages: complete on the retained features ----
GRAMBANK_query <- grambank.feature(c('gb020','gb021','gb022','gb023','gb028','gb030','gb031','gb035','gb036','gb037','gb042','gb043','gb044','gb051','gb052','gb053',
                                     'gb054','gb065','gb070','gb071','gb072','gb073','gb079','gb080','gb082','gb083','gb084','gb086','gb089','gb090','gb091','gb092',
                                     'gb093','gb094','gb107','gb121','gb130','gb131','gb137','gb138','gb171','gb172','gb186','gb192','gb196','gb197','gb316','gb318',
                                     'gb321','gb415'),na.rm = FALSE)

relatedfamilies <- unique(GRAMBANKdf_PH_maximized$Family_name)
relatedmacroareas <- unique(GRAMBANKdf_PH_maximized$Macroarea)

GRAMBANKdf_query_with_familyname <- GRAMBANK_query %>%
  left_join(languages %>% select(ID, Family_name, Macroarea, Name), by = c("glottocode" = "ID"))

# Keep languages outside the related families but within the same macroareas
GRAMBANKdf_unrelated <- GRAMBANKdf_query_with_familyname %>%
  filter(!Family_name %in% c('Austronesian')) %>%
  filter(Macroarea %in% relatedmacroareas) %>%
  select(glottocode,GB020,GB021,GB022,GB023,GB028,GB030,GB031,GB035,GB036,GB037,GB042,GB043,GB044,GB051,GB052,GB053,
         GB054,GB065,GB070,GB071,GB072,GB073,GB079,GB080,GB082,GB083,GB084,GB086,GB089,GB090,GB091,GB092,GB093,GB094,GB107,
         GB121,GB130,GB131,GB137,GB138,GB171,GB172,GB186,GB192,GB196,GB197,GB316,GB318,GB321,GB415,longitude,latitude,Family_name,Name,Macroarea)

GRAMBANKdf_unrelated %>%
  group_by(Macroarea) %>%
  summarise(Total = n()) %>%
  arrange(desc(Total))

GRAMBANKdf_unrelated %>%
  group_by(Family_name) %>%
  summarise(Total = n()) %>%
  arrange(desc(Total))

map.feature(GRAMBANKdf_unrelated$Name)

# Choose the unrelated control set by distance/diversity
UnrelatedLangs <- c("Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Ugaritic","Liv",
                    "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")

GRAMBANKdf_unrelated_sampled <- GRAMBANKdf_unrelated %>%
  filter(Name %in% UnrelatedLangs)

map.feature(GRAMBANKdf_unrelated_sampled$Name)

# ---- 7. Harmonize columns/types and combine Philippine + unrelated sets ----
colnames(GRAMBANKdf_unrelated_sampled)
colnames(GRAMBANKdf_PH_maximized)

GRAMBANKdf_unrelated_sampled <- GRAMBANKdf_unrelated_sampled %>%
  rename(
    Language_ID = glottocode,
    language = Name,
    Longitude = longitude,
    Latitude = latitude
  )

GRAMBANKdf_unrelated_sampled <- GRAMBANKdf_unrelated_sampled %>%
  mutate(Language_Type = "Unrelated Language")

GRAMBANKdf_PH_maximized <- GRAMBANKdf_PH_maximized %>%
  mutate(Language_Type = case_when(
    language %in% c("English", "Japanese") ~ "Language of Interest",
    TRUE                                   ~ "Philippine Language"
  ))

GRAMBANKdf_unrelated_sampled[feature_cols] <- lapply(GRAMBANKdf_unrelated_sampled[feature_cols], as.numeric)

GRAMBANKdf <- bind_rows(GRAMBANKdf_PH_maximized, GRAMBANKdf_unrelated_sampled)

map.feature(GRAMBANKdf$language,
            GRAMBANKdf$Language_Type)

# Write is intentionally disabled: GRAMBANKdf_PH.csv is hand-curated to add
# Spanish (from WALS), so we read the curated file back rather than overwrite it.
# write_csv(GRAMBANKdf, here("data", "GRAMBANKdf_PH.csv"))
GRAMBANKdf <- read.csv(here("data", "GRAMBANKdf_PH.csv"))

# ---- 8. Per-feature frequency + inverse-document-frequency weights ----
feature_cols <- intersect(colnames(GRAMBANKdf), colnames(GRAMBANK_query))

GRAMBANK_freq <- GRAMBANK_query %>%
  select(feature_cols)

# Frequency/IDF per feature-value across the 818 GRAMBANK query languages
gramfeature_freq <- GRAMBANK_freq %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "value") %>%
  group_by(feature, value) %>%
  summarise(n_languages = n(), .groups = "drop") %>%
  mutate(
    freq = n_languages / 818,
    IDF = log(1 / freq)
  )

write_csv(gramfeature_freq, here("data", "gramfeature_freq.csv"))
