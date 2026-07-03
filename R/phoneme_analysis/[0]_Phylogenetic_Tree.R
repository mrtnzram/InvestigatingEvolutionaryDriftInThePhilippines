# =============================================================================
# [0] Phoneme Analysis — Prune the phylogenetic tree to the study languages
#
# Reads the MCC tree (data/mcc.tree), reconciles its tip labels with the
# Ruhlen/phoneme language names, and prunes it to the Philippine study set.
#
# DEPENDENCY / RUN ORDER:
#   - REQUIRES `Ph_Languages` from PART A of [0]_CREANZA_RUHLENdatabase.R
#     (run that first, down to its END OF PART A banner).
#   - PRODUCES `Ph_Languages_pruned`, consumed by PART B of
#     [0]_CREANZA_RUHLENdatabase.R.
#   - ALSO PRODUCES `tree_pruned` and `tree_df_matched`, consumed by
#     [6]_PHONEME_PGLS.R.
# =============================================================================

library(tibble)
library(tidyverse)
library(ape)
library(stringr)
library(ggplot2)
library(here)

stopifnot(
  "Run PART A of [0]_CREANZA_RUHLENdatabase.R first: `Ph_Languages` is not defined." =
    exists("Ph_Languages")
)

tree <- read.nexus(here('data','mcc.tree'))

tree_df <- tibble(
  original = tree$tip.label
) %>%
  mutate(
    sanitized = str_remove(original, "_[0-9]+$") %>%
      str_replace_all("_", " ") %>%
      str_squish()
  )

lookup <- tribble(
  ~tree, ~ph,
  
  # exact matches
  "Agta", "Agta",
  "Gaddang", "Gaddang",
  "Ibanag", "Ibanag",
  "Ilokano", "Ilokano",
  "Balangaw", "Balangaw",
  "Inibaloi", "Inibaloi",
  "Iraya", "Iraya",
  "Binukid", "Binukid",
  "Mamanwa", "Mamanwa",
  "Hanunoo", "Hanunoo",
  "Buhid", "Buhid",
  "Tagalog", "Tagalog",
  "Kalagan", "Kalagan",
  "Cebuano", "Cebuano",
  "Hiligaynon", "Hiligaynon",
  "Maranao", "Maranao",
  "Tiruray", "Tiruray",
  "Pangasinan", "Pangasinan",
  "Kapampangan", "Kapampangan",
  "Yogad", "Yogad",
  
  # variants
  "Maguindanaon", "Magindanao",
  "Cuyonon", "Kuyunon",
  "Mansaka", "Mansakan",
  "Sangil Saragani Islands", "Sangil",
  "Tausug Jolo Dialect", "Tausug",
  "Samar-Leyte", "Waray",
  
  # Agta / Dumagat
  "Atta Pamplona", "Atta",
  "Dumagat Casiguran", "Casiguran Dumagat",
  
  # Isneg
  "Isneg Dibagat-Kabugao-Isneg", "Isnag",
  
  # Itneg
  "Itneg Binongan", "Itneg",
  
  # Kalinga
  "Kalinga Guinaang Lubuagan Dialect", "North Kalinga",
  "Kalinga Minangali", "Central Kalinga",
  "Kalinga Southern", "South Kalinga",
  
  # Ifugao
  "Ifugao Amganad", "Central Ifugao",
  "Ifugao Batad", "East Ifugao",
  
  # Bontok
  "Bontok Guina-ang", "Central Bontok",
  
  # Kankanaey
  "Kankanay Northern", "North Kankanay",
  
  # Kallahan
  "Kallahan Kayapa Proper", "North Kallahan",
  
  # Ilongot
  "Ilongot Kakidugen", "Ilongot",
  
  # Ivatan
  "Ivatan Basco Dialect", "Ivatanen",
  
  # Sambal
  "Sambal", "Tina",
  
  # Manobo
  "Manobo Ata up-river", "Ata",
  "Manobo Ata down-river", "Ata",
  "Manobo Dibabawon", "Dibabawon",
  "Manobo Ilianen Kibudtungan Dialect", "Ilianen",
  "Manobo Tigwa Iglogsad Dialect", "Tigwa",
  "Manobo Western Bukidnon", "Western Bukidnon",
  "Manobo Sarangani Kayaponga Dialect", "Sarangani Manobo",
  "Manobo Kalamansig Cotabato Paril Dialect", "Kalamansig",
  
  # Bilaan
  "Bilaan Koronadal", "Cotabato Bilaan",
  "Bilaan Sarangani", "Sarangani Bilaan",
  
  # Batak / Tagbanwa
  "Palawan Batak", "Batak",
  "Batak Palawan", "Batak",
  "Tagbanwa Kalamian Coron Island Dialect", "Kalamianen",
  "Tagbanwa Aborlan Dialect", "Tagbanwa",
  
  # Bikol
  "Naga Bikol", "Naga",
  
  # Subanon
  "Subanon Siocon", "Siocon Subanon",
  "Subanun Sindangan", "Sindangan",
  
  # Tboli
  "Tboli Tagabili", "Tboli"
)

Ph_Languages <- setdiff(Ph_Languages, "Ga-dang")

tree_df_matched <- tree_df %>%
  left_join(lookup, by = c("sanitized" = "tree")) %>%
  mutate(ph = coalesce(ph, sanitized))

setdiff(Ph_Languages, tree_df_matched$ph)

# --- Pruning -----

tips_keep <- tree_df_matched %>%
  filter(ph %in% Ph_Languages) %>%
  pull(original)

tree_pruned <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_keep)
)

tip_df <- tree_df_matched %>%
  select(original, ph)


Ph_Languages_pruned <- tree_df_matched %>%
  filter(ph %in% Ph_Languages) %>%
  pull(ph)


plot(tree_pruned, cex = 0.4)

