# =============================================================================
# [4] All domains — geography vs. ancestry, consolidated
#
# Runs each domain's [4]_*_PVR.R and stacks their results tables into one.
#
# Run order matters per domain, not across them. The phoneme, grammar and
# cognate [4] scripts open with stopifnot(exists("tree_pruned")), so this
# sources that domain's [0]_Phylogenetic_Tree.R into the same environment first.
# Phoneme's and grammar's tree scripts in turn need `Ph_Languages` (and, for
# grammar, `GRAMBANKdf_PH_maximized`) from PART A of their database script, so
# those are sourced ahead of the tree via part_a() — PART A only, since PART B
# asserts the tree has already run. Cognate's tree script reads CSVs and needs
# no predecessor. [4]_GENETIC_PVR.R needs no tree at all: it sources
# R/shared/select_moran_eigenvectors.R itself and reads PLINK principal
# components in place of phylogenetic eigenvectors.
#
# NOTE: sourcing a [0]_Phylogenetic_Tree.R rewrites that domain's
# *_subgroup_lookup.csv and *_phylo_dist_matrix.csv and re-saves its tree
# figure. Idempotent from the committed inputs, but it is not a no-op.
#
# Two schema notes:
#   - Cognate and genetic fit the model twice, with the zero-similarity rows
#     removed and kept, so they contribute 4 rows each and carry a `zeros`
#     column; phoneme and grammar contribute 2 rows and get NA there.
#   - Genetic's `predictor` reads "PC eigenvectors, migration distance
#     scrubbed" where the others read "phylo eigenvectors, ...". That column is
#     left as written, so the distinction survives the merge.
#
# Inputs:  R/phoneme_analysis/[0]_CREANZA_RUHLENdatabase.R  (PART A)
#          R/grammar_analysis/[0]_GRAMBANKdatabase.R        (PART A)
#          R/{phoneme,grammar,cognate}_analysis/[0]_Phylogenetic_Tree.R
#          R/{phoneme,grammar,cognate,genetic}_analysis/[4]_*_PVR.R
# Outputs: data/pvr/geo_vs_ancestry.csv
#          (the regression figures are still written by the domain scripts)
# =============================================================================

library(here)
source(here("R", "shared", "consolidate_domains.R"))

# Each entry is sourced in order into one environment: prerequisites first.
SCRIPTS <- list(
  phoneme = list(part_a(here("R", "phoneme_analysis", "[0]_CREANZA_RUHLENdatabase.R")),
                 here("R", "phoneme_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "phoneme_analysis", "[4]_PHONEME_PVR.R")),
  grammar = list(part_a(here("R", "grammar_analysis", "[0]_GRAMBANKdatabase.R")),
                 here("R", "grammar_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "grammar_analysis", "[4]_GRAMMAR_PVR.R")),
  cognate = list(here("R", "cognate_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "cognate_analysis", "[4]_COGNATE_PVR.R")),
  genetic = list(here("R", "genetic_analysis", "[4]_GENETIC_PVR.R"))
)

rows <- lapply(names(SCRIPTS), function(domain) {
  message("\n[", domain, "]")
  env <- run_domain(SCRIPTS[[domain]])
  harvest(env, paste0(toupper(domain), "_geo_vs_ancestry"), domain)
})

geo_vs_ancestry <- write_consolidated(
  rows, here("data", "pvr", "geo_vs_ancestry.csv")
)
