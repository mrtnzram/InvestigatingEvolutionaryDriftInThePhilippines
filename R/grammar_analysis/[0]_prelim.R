# =============================================================================
# [0] Grammar Analysis — Preliminary coverage inspection & matrix reduction
#
# Diagnostic script, run once before [0]_GRAMBANKdatabase.R. Builds a
# provisional Philippine + interest-language GRAMBANK feature matrix (Spanish
# excluded — GRAMBANK has no Spanish coverage, see [0]_GRAMBANKdatabase.R),
# inspects feature/language coverage, and iteratively prunes sparse
# features/languages until the matrix stabilizes. The settled thresholds
# (feature_thresh, language_thresh) and the resulting retained feature set are
# what [0]_GRAMBANKdatabase.R applies (as a single deterministic pass, not a
# re-run of this exploration) to build the pipeline's real feature matrix.
#
# Input:   data/languages.csv, data/values.csv
# Outputs: none (diagnostic; console/plot output only)
# Next:    [0]_GRAMBANKdatabase.R
# =============================================================================

library(dplyr)
library(tidyverse)
library(readr)
library(glue)
library(here)

# ---- 1. Load GRAMBANK languages + values ----
languages <- read_csv(here("data", "languages.csv"))
values <- read_csv(here("data", "values.csv"))

# ---- 2. Filter to Philippine languages (by lat/long) plus English/Japanese ----
# Spanish is added later, hard-coded from WALS (see [0]_GRAMBANKdatabase.R) —
# it is absent from GRAMBANK, so it contributes nothing to this coverage pass.
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

# ---- 4. Inspect feature/language coverage to choose reduction thresholds ----
feature_counts <- sort(colSums(!is.na(GRAMBANKdf_PH[ , !(names(GRAMBANKdf_PH) %in% c("Language_ID", "language"))])), decreasing = TRUE)

plot(feature_counts,
     type = "b",
     pch = 19,
     col = "#3182bd",
     main = "Feature Coverage (Sorted)",
     xlab = "Features (ranked)",
     ylab = "Number of Languages with Data")

language_counts <- sort(rowSums(!is.na(GRAMBANKdf_PH[ , !(names(GRAMBANKdf_PH) %in% c("Language_ID", "language"))])), decreasing = TRUE)

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

  feature_cols <- names(GRAMBANKdf_PH)[!names(GRAMBANKdf_PH) %in% c("Language_ID", "language", "Family_name", "Macroarea", "Longitude", 'Latitude')]

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

# ---- 6. Report the settled reduction ----
# NA counts per snapshot, to confirm iter_2 is the stable choice reused by
# [0]_GRAMBANKdatabase.R (feature_thresh = 80, language_thresh = 50).
sapply(names(GRAMBANK_snapshots), function(nm) max(colSums(is.na(GRAMBANK_snapshots[[nm]]))))

GRAMBANKdf_PH_maximized <- GRAMBANK_snapshots$iter_2

message(
  "\nSettled reduction (feature_thresh = ", feature_thresh, ", language_thresh = ", language_thresh, "):",
  "\n  Retained features  : ", length(keep_features),
  "\n  Retained languages : ", nrow(GRAMBANKdf_PH_maximized), " (of ", nrow(grambank_values_ph %>% distinct(Language_ID)), " candidates)"
)
