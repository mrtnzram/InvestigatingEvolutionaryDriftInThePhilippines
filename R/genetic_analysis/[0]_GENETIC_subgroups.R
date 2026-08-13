# =============================================================================
# [0] Genetic Analysis — Map populations to Glottolog subgroups
#
# GENETIC_final.csv carries no linguistic column at all, unlike the phoneme/
# grammar/cognate pipelines, which each get a language -> subgroup -> colour
# lookup from their own [0]_Phylogenetic_Tree.R. This produces the genetic
# equivalent — one Glottolog language and its level-4 subgroup per population —
# without a tree: genetic has no ABVD mcc.tree bridge (only 70/115 populations
# have any tree tip at all, and the missing 45 skew toward the strongest
# Spanish-admixture outliers), so [4]_GENETIC_PVR.R stays blocked on that
# separately. This script only feeds the FEEMS map colouring in [6]/[7].
#
# `population` names are compressed CamelCase / concatenated-word forms
# (`ManoboAgusan`, `AytaMagbukon`) that rarely match a Glottolog language name
# outright. Resolution runs five ordered tiers, each only applied to what the
# previous tier left unresolved:
#   1. direct match to a Glottolog language name
#   2. de-camelCase (insert a space before each internal capital; treat "-"
#      as a space too) — turns `SamaBangingi` into `Sama Bangingi`
#   3. word-order swap for exactly-two-word de-camelCased names — Glottolog
#      uses "Location Group" order, several population names use the reverse
#      ("Group Location"): `ManoboAgusan` -> `Agusan Manobo`
#   4. a hand-built rename table for exact-but-irregular correspondences
#      (endonym/orthography variants, or a named subgroup member standing in
#      for its Glottolog parent language)
#   5. a hand-built rename table for approximate correspondences — either a
#      Glottolog dialect-level node (no language-level entry exists) or the
#      nearest attested relative where no Glottolog node exists at all. These
#      are flagged (match_tier = 5) rather than silently mixed in with tier 4,
#      because they are genuine judgment calls, not lookups.
#
# Every resolved name is looked up through one helper (resolve_glottolog)
# rather than bare lingtypology::aff.lang(), because aff.lang() silently
# matches non-PH homonyms (there is an Indonesian "Batak" distinct from the
# Palawan "Batak") and family-level nodes, and it returns no glottocode or
# level to audit the match with.
#
# Input:   data/GENETIC_final.csv (population), data/subgroup_palette.csv
#          (canonical subgroup -> colour, from R/shared/subgroup_palette.R —
#          run that first)
# Outputs: data/GENETIC_subgroup_lookup.csv (population, glottocode,
#          glottolog_name, level, subgroup, colour, match_tier, note)
# Next:    [6]_feems_plot_GENETIC_span.R
# =============================================================================

library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(here)
library(lingtypology)

GENETIC_final <- read.csv(here("data", "GENETIC_final.csv"), stringsAsFactors = FALSE)
pops <- GENETIC_final$population

stopifnot(
  "GENETIC_final.csv$population has duplicates." = !anyDuplicated(pops)
)

# ── 1. Glottolog resolver ────────────────────────────────────────────────
# One row per Glottolog language name, PH-scoped rows preferred over any
# homonym elsewhere; falls back to the full table for dialect-level nodes
# (Umayam, Tigwa, Kabinga'an), whose `countries` field is blank.
gl_all <- lingtypology::glottolog %>%
  mutate(ph = str_detect(coalesce(countries, ""), "PH")) %>%
  arrange(desc(ph)) %>%
  distinct(language, .keep_all = TRUE)

sub4 <- function(affiliation, name, level) {
  mapply(function(a, n, l) {
    parts <- str_trim(str_split(a, ",")[[1]])
    if (l == "family") parts <- c(parts, n)
    if (length(parts) >= 4) parts[4] else tail(parts, 1)
  }, affiliation, name, level)
}

resolve_glottolog <- function(names) {
  tibble(glottolog_name = names) %>%
    left_join(
      gl_all %>% transmute(glottolog_name = language, glottocode, level,
                            subgroup = sub4(affiliation, language, level)),
      by = "glottolog_name"
    )
}

# ── 2. Automated tiers 1-3 ───────────────────────────────────────────────
decamel  <- str_replace_all(str_replace_all(pops, "-", " "), "(?<=[a-z])(?=[A-Z])", " ")
n_words  <- str_count(decamel, "[A-Za-z]+")
swapped  <- ifelse(n_words == 2,
                    vapply(strsplit(decamel, " "), \(x) paste(x[2], x[1]), character(1)),
                    NA_character_)

hit1 <- pops %in% gl_all$language
hit2 <- !hit1 & decamel %in% gl_all$language
hit3 <- !hit1 & !hit2 & !is.na(swapped) & swapped %in% gl_all$language

auto_name <- case_when(
  hit1 ~ pops,
  hit2 ~ decamel,
  hit3 ~ swapped,
  TRUE ~ NA_character_
)
auto_tier <- case_when(hit1 ~ 1L, hit2 ~ 2L, hit3 ~ 3L, TRUE ~ NA_integer_)

stopifnot(
  "Automated tiers 1-3 no longer resolve the expected 46 populations — Glottolog build may have changed." =
    sum(!is.na(auto_name)) == 46
)

# ── 3. Manual tiers 4-5 ──────────────────────────────────────────────────
# Tier 4: exact Glottolog language, reached via an orthographic/endonym
# variant or by naming a member language that stands in for its Glottolog
# parent. Tier 5: approximate — a dialect-level Glottolog node with no
# language-level entry, or (marked with a note) the nearest attested relative
# where the Glottolog build has no node for the population at all.
manual <- tribble(
  ~population,              ~glottolog_name,             ~match_tier, ~note,
  "Itbayaten",               "Itbayat",                    4L, NA,
  "Ivatan",                  "Ibatan",                      4L, "no 'Ivatan' language row in this Glottolog build, only Basco/Southern Ivatan dialects",
  "Ilocano",                 "Iloko",                       4L, NA,
  "Gddang",                  "Ga'dang",                     4L, "distinct from Gaddang (resolved separately by tier 1)",
  "Itawis",                  "Itawit",                      4L, NA,
  "Apayao",                  "Isnag",                       4L, NA,
  "Bontoc",                  "Central Bontoc",              4L, NA,
  "Bugkalot",                "Ilongot",                     4L, NA,
  "Sambal",                  "Tinà Sambal",                 4L, NA,
  "Kapampangan",             "Pampanga",                    4L, NA,
  "Bicolano",                "Coastal-Naga Bikol",          4L, NA,
  "Hilgaynon",               "Hiligaynon",                  4L, NA,
  "Sulodnon",                "Sulod",                       4L, NA,
  "Agutaya",                 "Agutaynen",                   4L, NA,
  "Palawano",                "Central Palawano",            4L, NA,
  "Meranao",                 "Maranao",                     4L, NA,
  "Maguindanaon",            "Maguindanao",                 4L, NA,
  "Teduray",                 "Tiruray",                     4L, NA,
  "BagoboKlata",             "Giangan",                     4L, NA,
  "Tagakaolo",               "Tagakaulu Kalagan",           4L, NA,
  "Mamanwa",                 "Minamanwa",                   4L, NA,
  "Subanon",                 "Western Subanon",             4L, NA,
  "Cinamiguin",              "Cinamiguin Manobo",           4L, NA,
  "Dibabaon",                "Dibabawon Manobo",            4L, NA,
  "Obo",                     "Obo Manobo",                  4L, NA,
  "OvuManuvo",               "Obo Manobo",                  4L, NA,
  "ManoboRK",                "Rajah Kabunsuwan Manobo",     4L, NA,
  "ManoboTagabawa",          "Tagabawa",                    4L, NA,
  "ManoboDulangan",          "Cotabato Manobo",             4L, NA,
  "AtiPanay",                "Ati",                         4L, NA,
  "AtiNegros",               "Ati",                         4L, NA,
  "BinukidnonNegros",        "Southern Binukidnon",         4L, NA,
  "MangyanBuhid",            "Buhid",                       4L, NA,
  "MangyanHanunuo",          "Hanunoo",                     4L, NA,
  "MangyanIraya",            "Iraya",                       4L, NA,
  "BukidnonBinukid",         "Talaandig-Binukid",           4L, NA,
  "BukidnonTalaandig",       "Talaandig-Binukid",           4L, NA,
  "BukidnonHigaonon",        "Higaonon",                    4L, NA,
  "BukidnonManobo",          "Western Bukidnon Manobo",     4L, NA,
  "BukdinonMatigsalug",      "Matigsalug Manobo",           4L, NA,
  "DavaoMatigsalug",         "Matigsalug Manobo",           4L, NA,
  "SamaDilautTaluksangay",   "Southern Sama",               4L, NA,
  "SamaDilautMampang",       "Southern Sama",               4L, NA,
  "SamaDilautBonggao",       "Southern Sama",               4L, NA,
  "SamaDeyaBonggao",         "Central Sama",                4L, NA,
  "SamaBangingi",            "Balangingi",                  4L, NA,
  "AytaMagbukon",            "Bataan Ayta",                 4L, NA,
  "AytaMag-indi",            "Mag-Indi Ayta",               4L, NA,
  "AytaMag-antsi",           "Mag-Anchi Ayta",               4L, NA,
  "AytaSambal",              "Botolan Sambal",              4L, NA,
  "AgtaDupaningan",          "Dupaninan Agta",              4L, NA,
  "AgtaLabin",               "Central Cagayan Agta",        4L, NA,
  "AttaRizal",               "Faire Atta",                  4L, "Rizal, Cagayan is the former municipality of Faire",
  "AgtaMaddela",             "Casiguran-Nagtipunan Agta",   4L, NA,
  "AgtaCasiguran",           "Casiguran-Nagtipunan Agta",   4L, NA,
  "AgtaDumagat",             "Umiray Dumaget Agta",         4L, NA,
  "AgtaRemontado",           "Hatang Kayi",                 4L, NA,
  "AgtaLopez",               "Alabat Island Agta",          4L, NA,
  "AgtaManide",              "Camarines Norte Agta",        4L, NA,
  "AgtaIriga",               "Mt. Iriga Agta",              4L, NA,
  "AgtaIraya",               "Mt. Iraya Agta",              4L, NA,
  "BukidnonUmayamnon",       "Umayam",                      5L, "Glottolog dialect-level node, no language-level entry",
  "BukidnonTigwahanon",      "Tigwa",                       5L, "Glottolog dialect-level node, no language-level entry",
  "SamaKabingaan",           "Kabinga'an",                  5L, "Glottolog dialect-level node, no language-level entry",
  "MangyanBangon",           "Eastern Tawbuid",             5L, "no Bangon Mangyan row; Bangon groups with Tawbuid, same province",
  "Manguangan",              "Dibabawon Manobo",            5L, "no Mangguangan row in this Glottolog build; nearest attested relative, Davao de Oro neighbour",
  "Lambanguian",             "Tiruray",                     5L, "Lambangian are a Teduray-Dulangan community, Maguindanao",
  "AgtaMatnog",              "Mt. Iraya Agta",              5L, "Sorsogon Ayta exists but is level='unattested' with a truncated affiliation ('Austronesian (Unattested)')",
  "AgtaBulusan",             "Mt. Iraya Agta",              5L, "Sorsogon Ayta exists but is level='unattested' with a truncated affiliation ('Austronesian (Unattested)')"
)

stopifnot(
  "manual tribble has a population not present in GENETIC_final.csv." =
    all(manual$population %in% pops),
  "manual tribble duplicates a population." =
    !anyDuplicated(manual$population)
)

# ── 4. Assemble ───────────────────────────────────────────────────────────
resolved_auto <- resolve_glottolog(auto_name) %>%
  mutate(population = pops, match_tier = auto_tier, note = NA_character_) %>%
  filter(!is.na(match_tier))

resolved_manual <- resolve_glottolog(manual$glottolog_name) %>%
  mutate(population = manual$population, match_tier = manual$match_tier, note = manual$note)

GENETIC_subgroup_lookup <- bind_rows(resolved_auto, resolved_manual) %>%
  right_join(tibble(population = pops), by = "population") %>%
  left_join(read.csv(here("data", "subgroup_palette.csv")) %>% select(subgroup, colour),
            by = "subgroup") %>%
  select(population, glottocode, glottolog_name, level, subgroup, colour, match_tier, note) %>%
  arrange(population)

stopifnot(
  "GENETIC_subgroup_lookup does not cover all 115 populations 1:1." =
    nrow(GENETIC_subgroup_lookup) == length(pops) &&
    setequal(GENETIC_subgroup_lookup$population, pops),
  "Some populations have no subgroup." =
    !anyNA(GENETIC_subgroup_lookup$subgroup),
  "Some subgroups have no colour — rerun R/shared/subgroup_palette.R." =
    !anyNA(GENETIC_subgroup_lookup$colour)
)

write.csv(GENETIC_subgroup_lookup, here("data", "GENETIC_subgroup_lookup.csv"), row.names = FALSE)

# Tally for the run log: expect tiers 1-5 = 36 / 3 / 7 / 61 / 8.
print(count(GENETIC_subgroup_lookup, match_tier))
