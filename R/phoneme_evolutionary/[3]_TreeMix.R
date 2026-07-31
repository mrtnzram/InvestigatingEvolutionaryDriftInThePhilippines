# =================================================================
# Creates TreeMix inputs using the fixed Gray et al. phylogeny
#
# Phoneme filtering:
#   n_languages calculated across the 284-language Ruhlen set
#   retain phonemes present in > 2 languages
#
# TreeMix backbone:
#   172 tips from phoneme_evolutionary_tree.nwk
#
# Per-donor analysis:
#   172 Austronesian tree tips + 1 donor = 173 populations
#
# Inputs:
#   data/RUHLENdf_PH_evolutionary.csv
#   data/RUHLENdf_PH.csv
#   data/phoneme_freq_ruhlen_evolutionary.csv
#   data/phoneme_evolutionary_analysis_tree.rds
#   data/phoneme_evolutionary_tip_mapping.csv

# Outputs:
#   data/treemix/treemix_feature_set.csv
#
#   Per donor:
#     <donor>treemix_phonemes96.gz
#     <donor>treemix_pop_order.txt
#     <donor>treemix_language_map.csv
#     <donor>treemix_fixed_tree.nwk
#
#   run_treemix_donors.sh
# =================================================================

library(here)
library(tidyverse)
library(ape)

RUHLENdf <- read_csv(here("data", "RUHLENdf_PH_evolutionary.csv"), show_col_types = FALSE)

RUHLENdf_full <- read_csv(here("data", "RUHLENdf_PH.csv"), show_col_types = FALSE)

phoneme_freq <- read_csv(here("data", "phoneme_freq_ruhlen_evolutionary.csv"), show_col_types = FALSE)
tip_mapping <- read_csv(here("data", "phoneme_evolutionary_tip_mapping.csv"), show_col_types = FALSE)
analysis_objects <- readRDS(here("data", "phoneme_evolutionary_analysis_tree.rds"))
evolutionary_tree <- analysis_objects$analysis_tree

stopifnot(inherits(evolutionary_tree, "phylo"))
stopifnot(Ntip(evolutionary_tree) == 172)
stopifnot(!anyDuplicated(evolutionary_tree$tip.label))

message("Fixed TreeMix backbone contains ", Ntip(evolutionary_tree), " tips")

DONORS <- c("Spanish", "English", "Japanese")
OUT_DIR <- here("data", "treemix")
RESULT_DIR <- here("data", "treemix_results")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)



# Grab phoneme columns


phoneme_cols <- names(RUHLENdf)[str_detect(names(RUHLENdf), "^phoneme_[0-9]{3}$")]

stopifnot(length(phoneme_cols) == 728)



# Filter phonemes using 284-language prevalence calculation

stopifnot(all(c("phoneme", "n_languages") %in% names(phoneme_freq)))
stopifnot(!anyDuplicated(phoneme_freq$phoneme))
stopifnot(all(phoneme_freq$n_languages >= 0 & phoneme_freq$n_languages <= 284))

phoneme_filter <- phoneme_freq |>
  transmute(
    phoneme,
    n_present_prevalence_set = n_languages,
    n_absent_prevalence_set = 284L - n_languages,
    keep = n_languages > 2 & n_languages < 284
  )

eligible_phonemes <- phoneme_filter |> filter(keep) |> pull(phoneme)
keep_phonemes <- phoneme_cols[phoneme_cols %in% eligible_phonemes]

stopifnot(length(keep_phonemes) == 96)

message("Retained ", length(keep_phonemes), " / ", length(phoneme_cols), " phonemes")

phoneme_filter |>
  filter(keep) |>
  mutate(feature_order = match(phoneme, keep_phonemes)) |>
  arrange(feature_order) |>
  write_csv(file.path(OUT_DIR, "treemix_feature_set.csv"))



# Sanitize TreeMix population names

sanitize_name <- function(x) {
  x |> str_trim() |> str_replace_all("[^A-Za-z0-9_.-]+", "_") |> str_replace_all("_+", "_") |> str_remove("^_") |> str_remove("_$")
}


# Match the 172 Gray tree tips to Ruhlen inventories
# Mapping file:
#   Name     = tree tip
#   language = Ruhlen language
#
# Tree tip names will be used as TreeMix population names.
# This avoids unnecessary renaming of the established phylogeny.


stopifnot(all(c("tip", "ruhlen_language") %in% names(tip_mapping)))

tree_mapping <- tip_mapping |>
  filter(tip %in% evolutionary_tree$tip.label) |>
  select(tree_tip = tip, language = ruhlen_language) |>
  distinct()

stopifnot(n_distinct(tree_mapping$tree_tip) == 172)
stopifnot(setequal(tree_mapping$tree_tip, evolutionary_tree$tip.label))



# Join the 172 tree tips to Ruhlen phoneme inventories

# The final row order is forced to match evolutionary_tree$tip.label.

backbone <- tibble(tree_tip = evolutionary_tree$tip.label) |>
  left_join(tree_mapping, by = "tree_tip") |>
  left_join(RUHLENdf, by = "language")

stopifnot(nrow(backbone) == 172)
stopifnot(!any(is.na(backbone$language)))
stopifnot(!anyDuplicated(backbone$tree_tip))



# Verify all 96 phonemes are available in the backbone


missing_backbone_phonemes <- setdiff(keep_phonemes, names(backbone))

if (length(missing_backbone_phonemes) > 0) stop("Missing backbone phonemes: ", paste(missing_backbone_phonemes, collapse = ", "))



# Create TreeMix names from the existing tree-tip labels

backbone <- backbone |> mutate(treemix_name = sanitize_name(tree_tip))

stopifnot(!anyDuplicated(backbone$treemix_name))


# Create fixed TreeMix backbone topology

# TreeMix receives the Gray topology but estimates its own branch lengths.

treemix_backbone_tree <- evolutionary_tree
treemix_backbone_tree$tip.label <- sanitize_name(treemix_backbone_tree$tip.label)
treemix_backbone_tree$edge.length <- NULL

stopifnot(Ntip(treemix_backbone_tree) == 172)
stopifnot(setequal(treemix_backbone_tree$tip.label, backbone$treemix_name))

write.tree(treemix_backbone_tree, file = file.path(OUT_DIR, "treemix_backbone_fixed_tree.nwk"))


# Validate donor dataset


missing_donor_phonemes <- setdiff(keep_phonemes, names(RUHLENdf_full))

if (length(missing_donor_phonemes) > 0) stop("Full Ruhlen donor dataset is missing phonemes: ", paste(missing_donor_phonemes, collapse = ", "))

for (donor in DONORS) {
  n_donor <- RUHLENdf_full |> filter(language == donor) |> nrow()
  if (n_donor != 1) stop("Expected exactly one ", donor, " row in RUHLENdf_full; found ", n_donor)
}


# Create donor-specific TreeMix inputs

write_donor_treemix <- function(donor) {
  
  donor_row <- RUHLENdf_full |> filter(language == donor)
  stopifnot(nrow(donor_row) == 1)
  
  donor_name <- sanitize_name(donor)
  donor_id <- str_to_lower(donor)
  
  # ---------------------------------------------------------------
  # 172 tree-backed inventories
  # ---------------------------------------------------------------
  
  X_backbone <- backbone |> select(all_of(keep_phonemes)) |> as.matrix()
  
  # ---------------------------------------------------------------
  # Donor inventory
  # ---------------------------------------------------------------
  
  X_donor <- donor_row |> select(all_of(keep_phonemes)) |> as.matrix()
  
  # ---------------------------------------------------------------
  # Final TreeMix matrix: 173 populations × 96 phonemes
  # ---------------------------------------------------------------
  
  X <- rbind(X_backbone, X_donor)
  
  stopifnot(nrow(X) == 173)
  stopifnot(ncol(X) == 96)
  stopifnot(all(!is.na(X)))
  stopifnot(all(X %in% c(0, 1)))
  
  # ---------------------------------------------------------------
  # Population names
  # ---------------------------------------------------------------
  
  treemix_names <- c(backbone$treemix_name, donor_name)
  
  stopifnot(length(treemix_names) == 173)
  stopifnot(!anyDuplicated(treemix_names))
  
  # ---------------------------------------------------------------
  # Convert binary presence/absence to TreeMix allele-count format
  #
  # present = 1,0
  # absent  = 0,1
  #
  # Rows = phonemes
  # Columns = populations
  # ---------------------------------------------------------------
  
  treemix_matrix <- t(ifelse(X == 1, "1,0", "0,1"))
  
  matrix_file <- file.path(OUT_DIR, paste0(donor_id, "treemix_phonemes96.gz"))
  con <- gzfile(matrix_file, open = "wt")
  
  writeLines(paste(treemix_names, collapse = " "), con)
  writeLines(apply(treemix_matrix, 1, paste, collapse = " "), con)
  
  close(con)
  
  # ---------------------------------------------------------------
  # Population order file for residual plotting
  # ---------------------------------------------------------------
  
  pop_order_file <- file.path(OUT_DIR, paste0(donor_id, "treemix_pop_order.txt"))
  write_lines(treemix_names, pop_order_file)
  
  # ---------------------------------------------------------------
  # Language map
  # ---------------------------------------------------------------
  
  backbone_map <- backbone |>
    transmute(
      treemix_name,
      tree_tip,
      language,
      iso6393,
      Language_type
    )
  
  donor_map <- donor_row |>
    transmute(
      treemix_name = donor_name,
      tree_tip = NA_character_,
      language,
      iso6393,
      Language_type = "Language of Interest"
    )
  
  language_map <- bind_rows(backbone_map, donor_map)
  
  stopifnot(nrow(language_map) == 173)
  
  map_file <- file.path(OUT_DIR, paste0(donor_id, "treemix_language_map.csv"))
  write_csv(language_map, map_file)
  
  # ---------------------------------------------------------------
  # Fixed topology
  #
  # Attach donor as sister to the entire Austronesian backbone:
  #
  #       donor
  #         |
  #         +------ entire Gray Austronesian tree
  #
  # ---------------------------------------------------------------
  
  backbone_newick <- write.tree(treemix_backbone_tree) |> str_remove(";$")
  fixed_tree_text <- paste0("(", donor_name, ",", backbone_newick, ");")
  fixed_tree <- read.tree(text = fixed_tree_text)
  
  # ---------------------------------------------------------------
  # Mandatory TreeMix/tree consistency checks
  # ---------------------------------------------------------------
  
  stopifnot(Ntip(fixed_tree) == 173)
  stopifnot(!anyDuplicated(fixed_tree$tip.label))
  stopifnot(setequal(fixed_tree$tip.label, treemix_names))
  
  missing_from_tree <- setdiff(treemix_names, fixed_tree$tip.label)
  extra_in_tree <- setdiff(fixed_tree$tip.label, treemix_names)
  
  if (length(missing_from_tree) > 0) stop(donor, ": populations missing from fixed tree: ", paste(missing_from_tree, collapse = ", "))
  if (length(extra_in_tree) > 0) stop(donor, ": extra populations in fixed tree: ", paste(extra_in_tree, collapse = ", "))
  
  fixed_tree_file <- file.path(OUT_DIR, paste0(donor_id, "treemix_fixed_tree.nwk"))
  write.tree(fixed_tree, file = fixed_tree_file)
  
  # ---------------------------------------------------------------
  # Summary
  # ---------------------------------------------------------------
  
  message(
    "\n====================================\n",
    donor, "\n",
    "TreeMix populations: ", nrow(X), "\n",
    "Austronesian backbone: 172\n",
    "Donor populations: 1\n",
    "Phonemes: ", ncol(X), "\n",
    "Fixed-tree tips: ", Ntip(fixed_tree), "\n",
    "====================================\n"
  )
  
  tibble(
    donor = donor,
    n_backbone = 172,
    n_donor = 1,
    n_populations = nrow(X),
    n_phonemes = ncol(X),
    matrix_file = matrix_file,
    fixed_tree_file = fixed_tree_file,
    pop_order_file = pop_order_file,
    language_map_file = map_file
  )
}


# 
# Generate Spanish, English and Japanese inputs


manifest <- map_dfr(DONORS, write_donor_treemix)

write_csv(manifest, file.path(OUT_DIR, "treemix_manifest.csv"))

print(manifest)



# Script for treemix run 
# m = 0:
#   fixed Gray topology, no migration
# m = 1:
#   fixed Gray topology + one migration edge
# Donor is used as the outgroup/root.


shell_lines <- c(
  "#!/usr/bin/env bash",
  "",
  "set -euo pipefail",
  "",
  'mkdir -p "results/treemix"',
  ""
)

for (donor in DONORS) {
  
  donor_id <- str_to_lower(donor)
  donor_name <- sanitize_name(donor)
  
  matrix_file <- paste0("data/treemix/", donor_id, "treemix_phonemes96.gz")
  tree_file <- paste0("data/treemix/", donor_id, "treemix_fixed_tree.nwk")
  
  shell_lines <- c(
    shell_lines,
    paste0("# ==================== ", donor, " ===================="),
    "",
    paste0(
      "treemix -i \"", matrix_file,
      "\" -tf \"", tree_file,
      "\" -root \"", donor_name,
      "\" -o \"results/treemix/", donor_id,
      "_m0\" -m 0 -k 1 -noss"
    ),
    "",
    paste0(
      "treemix -i \"", matrix_file,
      "\" -tf \"", tree_file,
      "\" -root \"", donor_name,
      "\" -o \"results/treemix/", donor_id,
      "_m1\" -m 1 -k 1 -noss"
    ),
    "",
    ""
  )
}

run_script <- here("run_treemix_donors.sh")

writeLines(shell_lines, run_script)
Sys.chmod(run_script, mode = "0755")


# summary statistics

message(
  "\nTreeMix preparation complete.\n\n",
  "Feature prevalence set: 284 languages\n",
  "Retained phonemes: ", length(keep_phonemes), "\n",
  "Fixed phylogenetic backbone: ", Ntip(treemix_backbone_tree), " tips\n",
  "Per-donor TreeMix dataset: 173 populations × ", length(keep_phonemes), " phonemes\n\n",
  "Run from the project root with:\n",
  "./run_treemix_donors.sh\n"
)