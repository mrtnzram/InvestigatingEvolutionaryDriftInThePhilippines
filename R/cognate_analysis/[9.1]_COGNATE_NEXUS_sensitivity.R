# =============================================================================
# [9.1] Cognate Analysis — Sensitivity: alternate NEXUS cognate-coding matrices
#
# Internal robustness check for [9]_COGNATE_PROCRUSTES.R's PCA-vs-geography
# result. Re-runs the identical extraction (locate the MATRIX block, join via
# each taxon's trailing ABVD numeric ID, collapse multi-taxon glottocodes with
# 1 > 0 > NA priority, drop >20%-missing columns, impute the rest to 0) and
# PCA/Procrustes pipeline against three alternate ABVD NEXUS exports that
# differ only in upstream ascertainment-correction/character-combining
# choices, alongside the same baseline file [0]/[9] use, to see how much the
# correlation depends on that choice.
#
# Variants (all four share the same 202 taxa in the same order):
#   baseline              abvd_philippines.nex                     (NCHAR=9422)
#   asc                   abvd_philippines.asc.nex                 (NCHAR=9423)
#   ascwords              abvd_philippines.ascwords.nex             (NCHAR=9632)
#   nocombining_ascwords  abvd_philippines.nocombining.ascwords.nex (NCHAR=8894)
#
# Run order: independent of [0]_COGNATE_LOANS.R and [9]_COGNATE_PROCRUSTES.R —
# re-derives its own matrices rather than reading COGNATE_NEXUS_matrix.csv, so
# all four variants go through an identical code path for a fair comparison.
#
# Input:   data/cognate/abvd_philippines*.nex (all four),
#          data/cognate/PH_df.csv,
#          data/network_distance/COGNATE_final.csv
# Outputs: data/procrustes/COGNATE_NEXUS_sensitivity.csv,
#          figures/procrustes/COGNATE_NEXUS_sensitivity.png
# =============================================================================

library(adegenet)   # loads ade4 (dudi.pca, procuste, procuste.randtest)
library(tidyverse)
library(here)

dir.create(here("figures", "procrustes"), recursive = TRUE, showWarnings = FALSE)

MAX_MISSING <- 0.20   # same cutoff [0]_COGNATE_LOANS.R uses
N_PERM      <- 999

COGNATE_final <- read.csv(here("data", "network_distance", "COGNATE_final.csv"))
analysis_glottocodes <- COGNATE_final$glottocode

PH_df <- read_csv(here("data", "cognate", "PH_df.csv"), show_col_types = FALSE)
id_to_glottocode <- PH_df |>
  select(language_id, glottocode) |>
  separate_longer_delim(language_id, ", ") |>
  mutate(language_id = str_trim(language_id))

# ── 1. Extraction (mirrors [0]_COGNATE_LOANS.R's NEXUS section, parameterized
#      by file path so it can run over all four variants identically) ───────
extract_nexus_matrix <- function(nex_path, max_missing = MAX_MISSING) {
  nex_lines <- readLines(nex_path)
  matrix_start <- which(nex_lines == "MATRIX")
  stopifnot("Expected exactly one MATRIX block in the NEXUS file." = length(matrix_start) == 1)

  ntax_line <- nex_lines[str_detect(nex_lines, "NTAX=")][1]
  ntax <- as.integer(str_match(ntax_line, "NTAX=([0-9]+)")[, 2])
  data_lines <- nex_lines[(matrix_start + 1):(matrix_start + ntax)]

  taxon_rows <- str_match(data_lines, "^(\\S+)\\s+(\\S+)$")
  stopifnot("A MATRIX row didn't parse as <taxon> <code string>." = !anyNA(taxon_rows))
  taxa    <- taxon_rows[, 2]
  strings <- taxon_rows[, 3]
  stopifnot("Character strings aren't all the same width." = length(unique(nchar(strings))) == 1)

  char_mat <- do.call(rbind, strsplit(strings, ""))
  char_mat[char_mat == "?"] <- NA
  mode(char_mat) <- "integer"
  rownames(char_mat) <- taxa
  colnames(char_mat) <- paste0("char_", seq_len(ncol(char_mat)))

  taxon_id <- str_extract(taxa, "[0-9]+$")
  taxon_glottocode <- tibble(taxon = taxa, language_id = taxon_id) |>
    left_join(id_to_glottocode, by = "language_id") |>
    filter(!is.na(glottocode))

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

  col_missing <- colMeans(is.na(nexus_mat))
  nexus_mat <- nexus_mat[, col_missing <= max_missing, drop = FALSE]
  nexus_mat[is.na(nexus_mat)] <- 0L

  list(matrix = nexus_mat, n_taxa_matched = nrow(taxon_glottocode), n_taxa_total = length(taxa),
       n_chars_total = ncol(char_mat), n_chars_kept = ncol(nexus_mat))
}

# ── 2. PCA + Procrustes fit (mirrors [9]_COGNATE_PROCRUSTES.R's §1-4) ───────
run_pca_procrustes <- function(nexus_mat, label) {
  X <- nexus_mat[rownames(nexus_mat) %in% analysis_glottocodes, , drop = FALSE]
  X <- X[, apply(X, 2, var) > 0, drop = FALSE]

  pca_fit <- dudi.pca(as.data.frame(X), center = TRUE, scale = FALSE, scannf = FALSE, nf = 2)
  scores_df <- tibble(glottocode = rownames(X), PC1 = pca_fit$li[, 1], PC2 = pca_fit$li[, 2])

  analysis_df <- scores_df |>
    left_join(COGNATE_final |> select(glottocode, longitude, latitude), by = "glottocode") |>
    arrange(glottocode)
  stopifnot("Some PCA rows are missing a coordinate." = !anyNA(analysis_df$longitude))

  pc_mat  <- as.data.frame(analysis_df[, c("PC1", "PC2")])
  geo_mat <- as.data.frame(analysis_df[, c("longitude", "latitude")])
  rownames(pc_mat) <- rownames(geo_mat) <- analysis_df$glottocode

  pr  <- procuste(pc_mat, geo_mat, scale = TRUE, nf = 2)
  ss  <- sum((as.matrix(pr$rotX) - as.matrix(pr$tabY))^2)
  prt <- procuste.randtest(pc_mat, geo_mat, nrepet = N_PERM)

  tibble(variant = label, n = nrow(X), n_features = ncol(X),
        ss_gower = ss, correlation = prt$obs, p_value = prt$pvalue, n_perm = N_PERM)
}

# ── 3. Run all four variants ─────────────────────────────────────────────────
variants <- tribble(
  ~label,                  ~file,
  "baseline",              "abvd_philippines.nex",
  "asc",                   "abvd_philippines.asc.nex",
  "ascwords",              "abvd_philippines.ascwords.nex",
  "nocombining_ascwords",  "abvd_philippines.nocombining.ascwords.nex"
)

results <- purrr::pmap(variants, function(label, file) {
  ex <- extract_nexus_matrix(here("data", "cognate", file))
  message(sprintf("[%s] %d/%d taxa matched -> %d languages; %d/%d characters kept (>%.0f%% missing dropped).",
                  label, ex$n_taxa_matched, ex$n_taxa_total, nrow(ex$matrix),
                  ex$n_chars_kept, ex$n_chars_total, MAX_MISSING * 100))
  run_pca_procrustes(ex$matrix, label)
}) |> bind_rows()

print(results)
write_csv(results, here("data", "procrustes", "COGNATE_NEXUS_sensitivity.csv"))

# ── 4. Comparison plot ───────────────────────────────────────────────────────
p_sensitivity <- ggplot(results, aes(x = fct_reorder(variant, correlation), y = correlation)) +
  geom_col(fill = "#4C72B0", width = 0.6) +
  geom_text(aes(label = sprintf("r = %.3f\np = %.3f", correlation, p_value)),
            vjust = -0.3, size = 3.2) +
  scale_y_continuous(limits = c(0, max(results$correlation) * 1.3)) +
  labs(x = NULL, y = "Procrustes correlation",
       title = "Cognate PCA vs. geography - sensitivity to NEXUS matrix variant") +
  theme_minimal() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7))
print(p_sensitivity)
ggsave(here("figures", "procrustes", "COGNATE_NEXUS_sensitivity.png"), p_sensitivity,
       width = 8, height = 6, units = "in", dpi = 300, bg = "white")
