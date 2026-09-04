# ================================================================
# [0] Genetic Analysis — qpAdm / MALDER results
# Loads the qpAdm and MALDER runs for the 112 Philippine populations, attaches
# the ethnolinguistic coordinates, and reports which MALDER admixture dates fall
# inside the Spanish colonial window. Metadata columns come from Larena et al. SI.
#
# Inputs:  data/Phil_admixture.tsv, data/philippine_ethnolinguistic_coords.csv
# Outputs: data/GENETIC_admixture.csv
# ================================================================

library(here)
library(tidyverse)

PH_genetic <- read_tsv(here('data', 'genetic', 'Phil_admixture.tsv'))
coords <- read_csv(here('data', 'genetic', 'philippine_ethnolinguistic_coords.csv'))

# override negative estimates to 0
PH_genetic <- PH_genetic |>
  mutate(balangao_admx = ifelse(balangao_admx < 0, 0, balangao_admx),
         aytamagbukon_admx = ifelse(aytamagbukon_admx < 0, 0, aytamagbukon_admx),
         span_admx = ifelse(span_admx < 0, 0, span_admx),
         eastasian_admx = ifelse(eastasian_admx < 0, 0, eastasian_admx))

# ManoboRajahKabunsuwan is ManoboRK in the coords table, and `region` there
# means island-group rather than qpAdm model block.
PH_genetic <- PH_genetic |>
  mutate(population = ifelse(population == "ManoboRajahKabunsuwan",
                             "ManoboRK", population)) |>
  rename(model_region = region)

stopifnot(
  "Population names in Phil_admixture.tsv are missing from the coords table." =
    all(PH_genetic$population %in% coords$population)
)

# merge metadata columns
PH_genetic <- left_join(
  coords, PH_genetic,
  by = "population"
  )

# The three qpAdm source populations have no target row of their own.
SOURCE_POPS <- c("AtaManobo", "AytaMagbukon", "Balangao")
stopifnot(
  "Unexpected populations lack admixture estimates — check the name reconciliation above." =
    setequal(PH_genetic$population[is.na(PH_genetic$model)], SOURCE_POPS)
)

# Admixture dates 

COLONIAL_LOWER <- 128   # 2026 - 1898 (end of Spanish colonial period)
COLONIAL_UPPER <- 461   # 2026 - 1565 (start of Spanish colonial period)
# "present" taken as 2026

categories <- c("span", "eastasian", "joint_ceu_eas")

for (category in categories) {
  years_col <- paste0("malder_", category, "_years")
  se_col    <- paste0("malder_", category, "_se_years")
  
  sub <- PH_genetic[!is.na(PH_genetic[[years_col]]),
                    c("population", years_col, se_col)]
  sub <- sub[order(sub[[years_col]]), ]
  
  cat(sprintf("\n=== %s (n = %d) ===\n", category, nrow(sub)))
  for (i in seq_len(nrow(sub))) {
    pop    <- sub$population[i]
    yrs    <- sub[[years_col]][i]
    se     <- sub[[se_col]][i]
    lo_est <- yrs - se
    hi_est <- yrs + se
    
    point_in_window   <- yrs >= COLONIAL_LOWER & yrs <= COLONIAL_UPPER
    interval_overlaps <- lo_est <= COLONIAL_UPPER & COLONIAL_LOWER <= hi_est
    
    flag <- if (point_in_window) {
      "overlaps"
    } else if (interval_overlaps) {
      "estimate out of bounds but standard error overlaps"
    } else {
      "estimate out of bounds"
    }
    
    cat(sprintf("  %-22s %8.1f +/- %7.1f y BP   [%s]\n", pop, yrs, se, flag))
  }
}


write.csv(PH_genetic, here("data", "genetic", "GENETIC_admixture.csv"), row.names = FALSE)
