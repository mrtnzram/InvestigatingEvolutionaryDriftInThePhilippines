# =============================================================================
# [5] All domains — MMRR, consolidated
#
# Runs each domain's [5]_*_MMRR.R and stacks their summaries into two tables:
# the joint two-predictor model and, for the three domains that have a
# phylogeny, the single-predictor (Mantel-equivalent) models.
#
# The four scripts are self-contained — they read only CSVs and define no
# globals — so each is sourced on its own into its own environment.
#
# Two schema notes:
#   - Genetic has no tree, so its joint model is single-predictor: it
#     contributes NA for beta_phylo/p_phylo and writes no single-model table at
#     all. `mmrr_single_results.csv` therefore covers three domains.
#   - The single-model `model` column is the plot title, so it reads
#     "Linguistic similarity vs ..." for phoneme/grammar and "Cognate
#     similarity vs ..." for cognate. It is recoded to a canonical key below.
#
# Permutation p-values floor at 1/(nperm+1) = 1e-04, so `p_geo = 1e-04` means
# "<= 0.0001", not an exact value. The domain scripts set.seed(1), so a re-run
# reproduces these exactly.
#
# Inputs:  R/{phoneme,grammar,cognate,genetic}_analysis/[5]_*_MMRR.R
# Outputs: data/mmrr/mmrr_results.csv        (joint model, 4 domains)
#          data/mmrr/mmrr_single_results.csv (single models, 3 domains)
#          (the per-domain *_sim_matrix.csv and the MMRR figures are still
#           written by the domain scripts themselves)
# =============================================================================

library(here)
library(dplyr)
source(here("R", "shared", "consolidate_domains.R"))

SCRIPTS <- list(
  phoneme = here("R", "phoneme_analysis", "[5]_PHONEME_MMRR.R"),
  grammar = here("R", "grammar_analysis", "[5]_GRAMMAR_MMRR.R"),
  cognate = here("R", "cognate_analysis", "[5]_COGNATE_MMRR.R"),
  genetic = here("R", "genetic_analysis", "[5]_GENETIC_MMRR.R")
)

# Only these three build a single-predictor table; genetic has no phylogeny.
SINGLE_DOMAINS <- c("phoneme", "grammar", "cognate")

# Canonical keys for the free-text plot titles, matched on the pair of terms
# rather than the leading noun so a retitled figure doesn't silently drop a row.
canonical_model <- function(x) {
  key <- dplyr::case_when(
    grepl("^Geograph.*vs Phylogenetic", x) ~ "geography_vs_phylogeny",
    grepl("similarity vs Geographic",   x) ~ "similarity_vs_geography",
    grepl("similarity vs Phylogenetic", x) ~ "similarity_vs_phylogeny",
    .default = NA_character_
  )
  if (any(is.na(key))) {
    stop("Unrecognised single-model title(s): ",
         paste(unique(x[is.na(key)]), collapse = "; "), call. = FALSE)
  }
  key
}

joint  <- list()
single <- list()

for (domain in names(SCRIPTS)) {
  message("\n[", domain, "]")
  env <- run_domain(SCRIPTS[[domain]])
  upper <- toupper(domain)

  joint[[domain]] <- harvest(env, paste0(upper, "_mmrr_results"), domain)

  if (domain %in% SINGLE_DOMAINS) {
    s <- harvest(env, paste0(upper, "_mmrr_single"), domain) %>%
      mutate(model = canonical_model(model))
    # Assert on content, not row order: all three models present exactly once.
    stopifnot(
      "single-model rows are not the expected three" =
        setequal(s$model, c("similarity_vs_geography",
                            "similarity_vs_phylogeny",
                            "geography_vs_phylogeny")) &&
        nrow(s) == 3L
    )
    single[[domain]] <- s
  }
}

mmrr_results <- write_consolidated(
  joint, here("data", "mmrr", "mmrr_results.csv")
)
mmrr_single_results <- write_consolidated(
  single, here("data", "mmrr", "mmrr_single_results.csv")
)
