# =============================================================================
# [9] All domains — Procrustes: PCA vs. geography, consolidated
#
# Runs each domain's [9]_*_PROCRUSTES.R and stacks their one-row summaries into
# a single table. The four scripts are self-contained — they read only CSVs and
# define no globals — so each is sourced on its own into its own environment.
# All four name their table `results_df`, which is exactly why run_domain()
# isolates them.
#
# Genetic contributes an extra `n_individual` column (its PCA figure plots 1028
# individuals while the fit itself runs at population grain); the other three
# get NA there.
#
# Inputs:  R/{phoneme,grammar,cognate,genetic}_analysis/[9]_*_PROCRUSTES.R
# Outputs: data/procrustes/procrustes_results.csv
#          (the per-domain *_procrustes_scores.csv and the Procrustes figures
#           are still written by the domain scripts themselves)
# =============================================================================

library(here)
source(here("R", "shared", "consolidate_domains.R"))

SCRIPTS <- list(
  phoneme = here("R", "phoneme_analysis", "[9]_PHONEME_PROCRUSTES.R"),
  grammar = here("R", "grammar_analysis", "[9]_GRAMMAR_PROCRUSTES.R"),
  cognate = here("R", "cognate_analysis", "[9]_COGNATE_PROCRUSTES.R"),
  genetic = here("R", "genetic_analysis", "[9]_GENETIC_PROCRUSTES.R")
)

rows <- lapply(names(SCRIPTS), function(domain) {
  message("\n[", domain, "]")
  env <- run_domain(SCRIPTS[[domain]])
  harvest(env, "results_df", domain)
})

procrustes_results <- write_consolidated(
  rows, here("data", "procrustes", "procrustes_results.csv")
)
