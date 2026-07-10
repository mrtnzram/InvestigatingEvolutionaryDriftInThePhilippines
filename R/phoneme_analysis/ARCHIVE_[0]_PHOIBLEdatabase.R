# =============================================================================
# [0] Phoneme Analysis — Build the PHOIBLE feature database
# Fetches PHOIBLE phoneme inventories, filters to the study languages, pivots to
# a binary presence matrix, and writes the feature matrix + IDF frequency table.
# Outputs: data/PHOIBLEdf_PH.csv, data/phoneme_freq.csv
# =============================================================================

library(lingtypology)
library(dplyr)
library(tidyverse)
library(readr)
library(proxy)
library(here)

# ---- 1. Load PHOIBLE and mark every attested phoneme as present (weight = 1) ----
PHOIBLEdf_full <- phoible.feature(source='all',na.rm = FALSE)

phoneme_presence <- PHOIBLEdf_full %>%
  filter(!is.na(phoneme)) %>%
  mutate(weight = 1) %>%
  select(iso6393, language, source, sound = phoneme, weight)

# ---- 2. Define the study language set (Philippine + interest + unrelated controls) ----
Ph_Languages <- c('Tiruray', 'Maranao', 'Central Sama', 'Hiligaynon', 'Aklanon', 'Filipino', 'Coastal-Naga Bikol','Tagalog','Pampanga','Pangasinan','Iloko','Dupaninan Agta',
               'Cebuano','Southern Sama')
Interest_Languages <- c('English','Spanish','Japanese')
UnrelatedLangs <- c("Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Ugaritic","Liv",
                    "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")

Languages <- c('Tiruray', 'Maranao', 'Central Sama', 'Hiligaynon', 'Aklanon', 'Filipino', 'Coastal-Naga Bikol','Tagalog','Pampanga','Pangasinan','Iloko','Dupaninan Agta',
               'Cebuano','Spanish','Japanese','English','Southern Sama',"Basque","Hungarian","Finnish","South Saami","Gagauz","Estonian","Turkish","Moksha","Ugaritic","Liv",
               "Korean","Mandarin Chinese","Udihe","Halh Mongolian","Sedang")

# ---- 3. Pick one source per language (prefer UPSID), then pivot to a binary matrix ----
PHOIBLEdf_PH <- phoneme_presence %>%
  group_by(iso6393, language, source, sound) %>%
  summarise(weight = max(weight), .groups = "drop") %>%
  filter(language %in% Languages) %>%
  group_by(language) %>%
  mutate(
    keep = if ("upsid" %in% source) source == "upsid" else source == first(source)
  ) %>%
  filter(keep) %>%
  ungroup() %>%
  pivot_wider(
    names_from = sound,
    values_from = weight,
    values_fill = list(weight = 0)
  ) %>%
  select(-any_of(c('keep')))

# ---- 4. Label languages by type and preview their geographic spread ----
PHOIBLEdf_PH <- PHOIBLEdf_PH %>%
  mutate(Language_type = case_when(
    language %in% Ph_Languages ~ "Philippine Language",
    language %in% Interest_Languages ~ "Language of Interest",
    language %in% UnrelatedLangs ~ "Unrelated Language"
  ))

map.feature(PHOIBLEdf_PH$language,
            PHOIBLEdf_PH$Language_type)

write_csv(PHOIBLEdf_PH, here("data", "PHOIBLEdf_PH.csv"))

# ---- 5. Compute global phoneme frequencies and inverse-document-frequency weights ----
PHOIBLEdf_clean <- PHOIBLEdf_full %>%
  group_by(language) %>%
  mutate(
    keep = if ("upsid" %in% source) source == "upsid" else source == first(source)
  ) %>%
  filter(keep) %>%
  ungroup()

phoneme_freq <- PHOIBLEdf_clean %>%
  distinct(language, phoneme) %>%
  filter(!is.na(phoneme)) %>%
  group_by(phoneme) %>%
  summarise(n_languages = n(), .groups = "drop") %>%
  mutate(
    freq = n_languages / n_distinct(PHOIBLEdf_clean$language),
    IDF = log(1 / freq)
  )

write_csv(phoneme_freq, here("data", "phoneme_freq.csv"))
