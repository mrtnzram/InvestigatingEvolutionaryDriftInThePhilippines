# ================================================================
# [0] Cognate Analysis — Spanish loanword counts + NEXUS cognate-coding matrix
# Two independent sections. (1) Loads the per-language ABVD loanword table,
# checks it against the 210 per-gloss loan flags, drops those flag columns,
# and min-max normalizes the loan count into the response the rest of the
# pipeline regresses on. (2) Extracts the binary cognate-set coding matrix
# out of the raw ABVD NEXUS file and collapses it to one row per glottocode,
# for [9]_COGNATE_PROCRUSTES.R's PCA.
#
# Inputs:  data/cognate/PH_final.csv,
#          data/cognate/abvd_philippines.nex,
#          data/cognate/PH_df.csv
# Outputs: data/cognate/COGNATE_loans.csv,
#          data/cognate/COGNATE_NEXUS_matrix.csv
# ================================================================

library(here)
library(tidyverse)

# read_csv, not read.csv: several gloss headers contain commas and slashes that
# check.names would mangle.
PH_cognate <- read_csv(here("data", "cognate", "PH_final.csv"),
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
  write.csv(file = here("data", "cognate", "COGNATE_loans.csv"), row.names = FALSE)


# ── 3. NEXUS cognate-coding matrix ──────────────────────────────────────────
# Extracts the 202-taxon x 9422-character binary cognate-set matrix out of the
# ABVD/King et al. 2024 NEXUS file (nobody else in this pipeline reads it),
# joins it down to glottocode via each taxon's trailing ABVD numeric ID —
# which matches PH_df.csv's language_id column exactly, so this is an exact
# ID join, not name-string matching — and writes a glottocode-keyed matrix for
# [9]_COGNATE_PROCRUSTES.R's PCA. Independent of sections 1-2 above: shares no
# inputs or outputs with the loanword-count logic.
MAX_MISSING <- 0.20   # drop characters missing in more than this fraction of languages

nex_lines <- readLines(here("data", "cognate", "abvd_philippines.nex"))
matrix_start <- which(nex_lines == "MATRIX")
stopifnot("Expected exactly one MATRIX block in the NEXUS file." = length(matrix_start) == 1)

ntax_line <- nex_lines[str_detect(nex_lines, "NTAX=")][1]
ntax <- as.integer(str_match(ntax_line, "NTAX=([0-9]+)")[, 2])
data_lines <- nex_lines[(matrix_start + 1):(matrix_start + ntax)]

# Each row is "<taxon label><whitespace><0/1/? code string>", no embedded
# whitespace in the code string itself.
taxon_rows <- str_match(data_lines, "^(\\S+)\\s+(\\S+)$")
stopifnot("A MATRIX row didn't parse as <taxon> <code string>." = !anyNA(taxon_rows))
taxa    <- taxon_rows[, 2]
strings <- taxon_rows[, 3]
stopifnot("Character strings aren't all the same width." = length(unique(nchar(strings))) == 1)

char_mat <- do.call(rbind, strsplit(strings, ""))
char_mat[char_mat == "?"] <- NA
mode(char_mat) <- "integer"
rownames(char_mat) <- taxa
colnames(char_mat) <- paste0("char_", seq_len(ncol(char_mat)))   # stable across the missingness filter below

PH_df <- read_csv(here("data", "cognate", "PH_df.csv"), show_col_types = FALSE)
id_to_glottocode <- PH_df |>
  select(language_id, glottocode) |>
  separate_longer_delim(language_id, ", ") |>
  mutate(language_id = str_trim(language_id))

taxon_id <- str_extract(taxa, "[0-9]+$")
taxon_glottocode <- tibble(taxon = taxa, language_id = taxon_id) |>
  left_join(id_to_glottocode, by = "language_id") |>
  filter(!is.na(glottocode))

message(sprintf("%d of %d NEXUS taxa matched a PH_df glottocode (%d distinct languages).",
                nrow(taxon_glottocode), length(taxa), n_distinct(taxon_glottocode$glottocode)))

# Collapse multiple taxa per glottocode with 1 > 0 > NA priority: a language's
# cell is 1 if any matched taxon attests it, else 0 if any explicitly rejects
# it, else NA only if every matched taxon is missing there.
merge_priority <- function(x) {
  if (any(x == 1L, na.rm = TRUE)) return(1L)
  if (any(x == 0L, na.rm = TRUE)) return(0L)
  NA_integer_
}

glottocodes <- unique(taxon_glottocode$glottocode)
nexus_mat <- matrix(NA_integer_, nrow = length(glottocodes), ncol = ncol(char_mat),
                    dimnames = list(glottocodes, colnames(char_mat)))
for (g in glottocodes) {
  rows <- char_mat[taxon_glottocode$taxon[taxon_glottocode$glottocode == g], , drop = FALSE]
  nexus_mat[g, ] <- if (nrow(rows) == 1) rows[1, ] else apply(rows, 2, merge_priority)
}

# Drop characters missing in too many languages, then impute the rest to 0 —
# matches the "absence" convention the wordform matrix already uses: a
# character with scattered "?"s but no observed "1" stays 0, not invented.
col_missing <- colMeans(is.na(nexus_mat))
nexus_mat <- nexus_mat[, col_missing <= MAX_MISSING, drop = FALSE]
nexus_mat[is.na(nexus_mat)] <- 0L

message(sprintf("NEXUS matrix: %d languages x %d characters kept (%d dropped above %.0f%% missing).",
                nrow(nexus_mat), ncol(nexus_mat), sum(col_missing > MAX_MISSING), MAX_MISSING * 100))

nexus_mat |>
  as.data.frame() |>
  rownames_to_column("glottocode") |>
  write.csv(file = here("data", "cognate", "COGNATE_NEXUS_matrix.csv"), row.names = FALSE)
