# =============================================================================
# [4] All domains — geography vs. pure phylogeny, consolidated
#
# Runs each domain's [4]_*_PVR.R and stacks their `pvr_results_table`s into one.
#
# Run order matters per domain, not across them. The phoneme, grammar and
# cognate [4] scripts open with stopifnot(exists("tree_pruned")), so this
# sources that domain's [0]_Phylogenetic_Tree.R into the same environment first.
# Phoneme's and grammar's tree scripts in turn need `Ph_Languages` (and, for
# grammar, `GRAMBANKdf_PH_maximized`) from PART A of their database script, so
# those are sourced ahead of the tree via part_a() — PART A only, since PART B
# asserts the tree has already run. Cognate's tree script reads CSVs and needs
# no predecessor. [4]_GENETIC_PVR.R needs no tree: its axis is PLINK PC1.
#
# NOTE: sourcing a [0]_Phylogenetic_Tree.R rewrites that domain's
# *_subgroup_lookup.csv and *_phylo_dist_matrix.csv and re-saves its tree
# figure. Idempotent from the committed inputs, but it is not a no-op.
#
# Schema: four rows per fit — phylogeny_vs_geography twice (one per term:
# signed geographic axis, cost to Manila; F/p/R^2 are the joint model's),
# A_y_vs_geography, B_y_vs_pure_phylogeny. Cognate and genetic fit twice, zeros
# removed and kept, so they contribute 8 rows and carry a `zeros` column.
# Genetic's `axis1_var_share` is NA (PC1 is not re-derived from a distance matrix).
#
# Inputs:  R/phoneme_analysis/[0]_CREANZA_RUHLENdatabase.R  (PART A)
#          R/grammar_analysis/[0]_GRAMBANKdatabase.R        (PART A)
#          R/{phoneme,grammar,cognate}_analysis/[0]_Phylogenetic_Tree.R
#          R/{phoneme,grammar,cognate,genetic}_analysis/[4]_*_PVR.R
# Outputs: data/pvr/geo_vs_ancestry.csv
#          data/pvr/pvr_alignment_checks.csv (correlation at each step; status FLAG = failed check)
#          (the figures are written by the domain scripts)
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

envs <- lapply(setNames(names(SCRIPTS), names(SCRIPTS)), function(domain) {
  message("\n[", domain, "]")
  run_domain(SCRIPTS[[domain]])
})

geo_vs_ancestry <- write_consolidated(
  lapply(names(envs), \(d) harvest(envs[[d]], "pvr_results_table", d)),
  here("data", "pvr", "geo_vs_ancestry.csv")
)
pvr_alignment_checks <- write_consolidated(
  lapply(names(envs), \(d) harvest(envs[[d]], "pvr_checks", d)),
  here("data", "pvr", "pvr_alignment_checks.csv")
)
