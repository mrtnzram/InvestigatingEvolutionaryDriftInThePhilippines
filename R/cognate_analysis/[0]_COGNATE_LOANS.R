# ================================================================
# [0] Cognate Analysis — Spanish loanword counts
# Loads the per-language ABVD loanword table, checks it against the 210 per-gloss
# loan flags, drops those flag columns, and min-max normalizes the loan count
# into the response the rest of the pipeline regresses on.
#
# Inputs:  data/cognate/PH_final.csv
# Outputs: data/cognate/COGNATE_loans.csv
# ================================================================

library(here)
library(tidyverse)

# read_csv, not read.csv: several gloss headers contain commas and slashes that
# check.names would mangle.
PH_cognate <- read_csv(here("data", "cognate", "initial_datasets", "PH_final.csv"),
                       show_col_types = FALSE)

META_COLS <- c("language_id", "language", "latitude", "longitude",
               "glottocode", "author", "number_of_entries", "source_url")

# Every remaining column is a per-gloss 0/1 loan flag.
gloss_cols <- setdiff(names(PH_cognate), c(META_COLS, "number_of_loans"))

stopifnot(
  "PH_final.csv is missing an expected metadata column." =
    all(META_COLS %in% names(PH_cognate)),
  "number_of_loans does not match the row sum of the per-gloss loan flags." =
    all(rowSums(PH_cognate[gloss_cols]) == PH_cognate$number_of_loans)
)

message(sprintf("%d languages x %d glosses; %d loan flags total.",
                nrow(PH_cognate), length(gloss_cols),
                sum(PH_cognate$number_of_loans)))


# ── 1. Key and coordinate checks ────────────────────────────────────────────
# Coordinates are required by [3]; glottocode is the join key to the phoneme and
# grammar tables.
stopifnot(
  "Some languages have no coordinates." =
    !any(is.na(PH_cognate$latitude) | is.na(PH_cognate$longitude)),
  "Glottocodes are not unique — merged wordlists should already be collapsed." =
    !any(duplicated(PH_cognate$glottocode))
)


# ── 2. Response variable ────────────────────────────────────────────────────
# Min-max normalized so a slope reads as fraction of the observed borrowing range
# per unit distance, matching [1]'s cosine convention. No standard error exists
# on this side, so unlike the genetic pipeline there are no per-language weights.
loans_min <- min(PH_cognate$number_of_loans)
loans_max <- max(PH_cognate$number_of_loans)

PH_cognate <- PH_cognate |>
  mutate(loans_norm = (number_of_loans - loans_min) / (loans_max - loans_min))

message(sprintf("number_of_loans range [%d, %d] -> loans_norm [0, 1]; %d languages at zero.",
                loans_min, loans_max, sum(PH_cognate$number_of_loans == 0)))

cat("\n=== Loan-count distribution ===\n")
print(table(PH_cognate$number_of_loans))

cat("\n=== Most-borrowing languages ===\n")
PH_cognate |>
  select(language, glottocode, number_of_entries, number_of_loans, loans_norm) |>
  arrange(desc(number_of_loans)) |>
  print(n = 20)


# number_of_entries is carried through unused: the count is kept raw by design,
# but it is what an exposure-offset model would need.
PH_cognate |>
  select(all_of(META_COLS), number_of_loans, loans_norm) |>
  write.csv(file = here("data", "cognate", "initial_datasets", "COGNATE_loans.csv"), row.names = FALSE)
