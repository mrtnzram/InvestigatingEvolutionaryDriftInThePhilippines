# =============================================================================
# [0] Grammar Evolutionary — Match the pulses tree to the GRAMBANK grammar data
#
# Pulls every Austronesian language in GRAMBANK for the shared 50-feature set,
# joins it to the pulses tree on Glottocode, prunes the tree to matched tips,
# resolves near-zero-divergence sibling conflicts, and paints the Philippine
# regime for [1]/[2]. The feature matrix stays MULTISTATE — GB065 and GB130 keep
# their {1,2,3} levels rather than being expanded into dummies the way
# [1]_GRAMMAR_cosine_similarity.R does for cosine similarity.
#
# Self-contained: unlike the phoneme and grammar database scripts there is no
# PART A / PART B handshake, because the CLDF dump is read directly off disk.
#
# Inputs:  data/shared/values.csv, data/shared/languages.csv,
#          data/phoneme/pulses_summary.trees, data/phoneme/pulses_languages.csv
#
# Outputs: data/grammar/grammar_evolutionary_tip_matrix.csv
#          data/grammar/grammar_evolutionary_tree.nwk
#          data/grammar/grammar_evolutionary_analysis_tree.rds
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ape)
library(here)
library(phytools)   # multi2di, paintSubTree, for the [1]/[2] analysis tree

source(here("R", "shared", "grambank_features.R"))

# Languages answering fewer than this many of the 50 features are dropped; the
# handful of cells still missing above the floor become ambiguous tip priors in
# [1]/[2] rather than being imputed to a state.
MIN_ANSWERED <- 48L

# Philippine bounding box, matching the phoneme and grammar database scripts.
PH_LAT <- c(4.5, 21)
PH_LON <- c(115, 128)

# Sabah languages: inside the bounding box, but North Bornean rather than
# Philippine, so including them drags the Philippine MRCA down to the
# Malayo-Polynesian root. Same reasoning as the phoneme script's
# 'Yakan'/'Timugon'/'Sangil' exclusion.
EXCLUDED_PH_TIPS <- c("Idaan", "TimugonMurut")

# Sibling tips below this divergence force runaway rates in [1]/[2] whenever
# they disagree on a feature.
MIN_SIBLING_DIVERGENCE <- 0.001

# ---- 1. GRAMBANK: Austronesian languages on the shared feature set ----------
# The cached CLDF dump, not lingtypology::grambank.feature() — that call hits
# the network and inner-joins down to languages with complete coverage on all
# 50 features, which would discard two thirds of the Austronesian set.
grambank_languages <- read_csv(
  here("data", "shared", "languages.csv"),
  show_col_types = FALSE
) |>
  filter(Family_name == "Austronesian") |>
  dplyr::select(glottocode = ID, gb_name = Name, latitude = Latitude, longitude = Longitude)

grambank_values <- read_csv(
  here("data", "shared", "values.csv"),
  col_select = c(Language_ID, Parameter_ID, Value),
  show_col_types = FALSE
) |>
  filter(
    Parameter_ID %in% GRAMBANK_FEATURE_IDS,
    Language_ID %in% grambank_languages$glottocode
  ) |>
  mutate(Value = as.numeric(na_if(Value, "?")))

# "?" is GRAMBANK's explicit not-known code and is already NA above; rows that
# never appear at all are equally unknown, so both arrive as NA after the pivot.
feature_matrix <- grambank_values |>
  pivot_wider(names_from = Parameter_ID, values_from = Value) |>
  dplyr::select(glottocode = Language_ID, any_of(GRAMBANK_FEATURE_IDS))

missing_features <- setdiff(GRAMBANK_FEATURE_IDS, names(feature_matrix))
if (length(missing_features) > 0) {
  feature_matrix[missing_features] <- NA_real_
}
feature_matrix <- feature_matrix |> dplyr::select(glottocode, all_of(GRAMBANK_FEATURE_IDS))

feature_matrix <- feature_matrix |>
  mutate(n_answered = rowSums(!is.na(across(all_of(GRAMBANK_FEATURE_IDS)))))

message(
  "GRAMBANK Austronesian languages: ", nrow(grambank_languages), " total, ",
  nrow(feature_matrix), " with at least one of the ", length(GRAMBANK_FEATURE_IDS),
  " features answered."
)

# ---- 2. Coverage filter ----------------------------------------------------
covered <- feature_matrix |> filter(n_answered >= MIN_ANSWERED)

message(
  "Coverage filter (>= ", MIN_ANSWERED, " of ", length(GRAMBANK_FEATURE_IDS),
  " features answered): ", nrow(covered), " of ", nrow(feature_matrix), " languages kept."
)

# ---- 3. Glottocode join to the pulses tree ---------------------------------
# GRAMBANK's languages.csv ID *is* the Glottocode, and pulses_languages.csv
# carries one, so this is a direct join. The phoneme pipeline needs gltc.iso()
# to bridge Ruhlen's iso6393 first; grammar does not.
tree <- read.nexus(here("data", "phoneme", "pulses_summary.trees"))

pulses_languages <- read_csv(
  here("data", "phoneme", "pulses_languages.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(tip = Name, glottocode = Glottocode) |>
  filter(tip %in% tree$tip.label)

matched <- pulses_languages |>
  inner_join(covered, by = "glottocode") |>
  inner_join(grambank_languages, by = "glottocode") |>
  mutate(
    is_philippine =
      latitude  > PH_LAT[1] & latitude  < PH_LAT[2] &
      longitude > PH_LON[1] & longitude < PH_LON[2] &
      !tip %in% EXCLUDED_PH_TIPS
  )

message(
  "\nPulses tips matched to a covered Austronesian GRAMBANK language: ", nrow(matched),
  " (", n_distinct(matched$glottocode), " distinct Glottocodes)."
)

excluded_present <- intersect(EXCLUDED_PH_TIPS, matched$tip)
if (length(excluded_present) > 0) {
  message(
    "  Bounding-box tips held out of the Philippine set as North Bornean: ",
    paste(excluded_present, collapse = ", ")
  )
}

# ---- 4. Collapse duplicate Glottocodes -------------------------------------
# GRAMBANK is keyed by Glottocode, so several tips under one Glottocode would
# carry byte-identical feature rows. Keeping them would pseudo-replicate the
# data on near-zero branches, so one representative tip is kept per Glottocode.
# `_D`-suffixed labels are the tree's dialect variants and are deprioritised.
duplicated_glottocodes <- matched |> count(glottocode) |> filter(n > 1)

if (nrow(duplicated_glottocodes) > 0) {
  dup_log <- matched |>
    semi_join(duplicated_glottocodes, by = "glottocode") |>
    arrange(glottocode, str_detect(tip, "_D$"), tip)
  message(
    "\n", nrow(duplicated_glottocodes), " Glottocode(s) carry more than one tip ",
    "(identical GRAMBANK rows; one representative kept each):\n",
    paste0(
      "  ", dup_log$glottocode, ": ", dup_log$tip,
      collapse = "\n"
    )
  )
}

matched_unique <- matched |>
  arrange(glottocode, str_detect(tip, "_D$"), tip) |>
  group_by(glottocode) |>
  slice(1) |>
  ungroup()

stopifnot(
  "one tip per Glottocode after collapsing" =
    !anyDuplicated(matched_unique$glottocode),
  "tip labels must stay unique" =
    !anyDuplicated(matched_unique$tip)
)

# ---- 5. Prune ---------------------------------------------------------------
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, matched_unique$tip))
tree_pruned$root.edge <- 0
tree_pruned <- multi2di(tree_pruned, random = FALSE)

# ---- 6. Resolve near-zero-divergence sibling conflicts ---------------------
# The keep-rule must be a total order, not just a ranking: the Lenakel /
# TannaSouthwest cherry is a zero-length pair whose members both answer all 50
# features, so coverage alone ties and the phoneme script's tie assertion would
# abort. Coverage decides first, then a Philippine tip is preferred over a
# background one, then the label breaks any remainder.
short_nodes <- unique(
  tree_pruned$edge[tree_pruned$edge.length < MIN_SIBLING_DIVERGENCE, 1]
)

sibling_clusters <- lapply(short_nodes, function(nd) {
  ii <- which(tree_pruned$edge[, 1] == nd)
  ch <- tree_pruned$edge[ii, 2]
  short_tips <- ch[
    tree_pruned$edge.length[ii] < MIN_SIBLING_DIVERGENCE & ch <= Ntip(tree_pruned)
  ]
  if (length(short_tips) < 2) return(NULL)
  tree_pruned$tip.label[short_tips]
})
sibling_clusters <- Filter(Negate(is.null), sibling_clusters)

sibling_conflicts <- bind_rows(lapply(seq_along(sibling_clusters), function(i) {
  matched_unique |>
    filter(tip %in% sibling_clusters[[i]]) |>
    dplyr::select(tip, n_answered, is_philippine) |>
    arrange(desc(n_answered), !is_philippine, tip) |>
    mutate(pair_id = i, keep = row_number() == 1L)
}))

if (nrow(sibling_conflicts) > 0) {
  stopifnot(
    "the sibling keep-rule must retain exactly one tip per cluster" =
      all(tapply(sibling_conflicts$keep, sibling_conflicts$pair_id, sum) == 1)
  )
  message(
    "\nZero/near-zero-divergence sibling conflicts resolved by coverage, then ",
    "Philippine status, then label:\n",
    paste0(
      "  ", sibling_conflicts$tip, " (", sibling_conflicts$n_answered,
      if_else(sibling_conflicts$is_philippine, ", Philippine", ""), ")",
      if_else(sibling_conflicts$keep, " <- kept", ""),
      collapse = "\n"
    )
  )
}

CONFLICTING_TIPS_DROPPED <- sibling_conflicts$tip[!sibling_conflicts$keep]

analysis_tree <- drop.tip(tree_pruned, CONFLICTING_TIPS_DROPPED)
analysis_tree$root.edge <- 0
analysis_tree <- multi2di(analysis_tree, random = FALSE)

stopifnot("analysis tree must be binary" = is.binary(analysis_tree))

# ---- 7. Philippine regime ---------------------------------------------------
philippine_tips  <- matched_unique$tip[matched_unique$is_philippine]
philippine_tips  <- intersect(philippine_tips, analysis_tree$tip.label)
philippine_node  <- getMRCA(analysis_tree, philippine_tips)
philippine_clade <- extract.clade(analysis_tree, philippine_node)

# `regime_tips` is the clade paintSubTree actually paints and the set [2] must
# split the tree on — NOT `philippine_tips`. The two differ whenever the MRCA
# sweeps in a non-Philippine tip, and splitting on the smaller set would leave
# the swept-in tip on both halves of the split.
regime_tips   <- philippine_clade$tip.label
swept_in_tips <- setdiff(regime_tips, philippine_tips)

# Philippine languages need not be monophyletic, so this is reported rather than
# asserted; the bound only catches a regime that has stopped being Philippine.
message(
  "\nPhilippine-language clade: ", length(regime_tips), " tips (",
  length(philippine_tips), " Philippine, ", length(swept_in_tips),
  " non-Philippine swept in by the smallest common ancestor)."
)
if (length(swept_in_tips) > 0) {
  message("  Non-Philippine tips in the clade: ", paste(swept_in_tips, collapse = ", "))
}

stopifnot(
  "the Philippine MRCA must not sweep in more non-Philippine tips than Philippine ones" =
    length(swept_in_tips) < length(philippine_tips)
)

regime_tree <- paintSubTree(
  tree      = analysis_tree,
  node      = philippine_node,
  state     = "Philippine",
  anc.state = "Background",
  stem      = FALSE
)

# ---- 8. Export for [1]/[2] --------------------------------------------------
# The tip matrix is written whole so [1]/[2] never re-derive it from the CLDF.
tip_matrix <- matched_unique |>
  filter(tip %in% analysis_tree$tip.label) |>
  dplyr::select(
    tip, glottocode, gb_name, is_philippine, n_answered,
    all_of(GRAMBANK_FEATURE_IDS)
  ) |>
  arrange(tip)

stopifnot(
  "tip matrix must be one row per analysis-tree tip" =
    nrow(tip_matrix) == Ntip(analysis_tree) && !anyDuplicated(tip_matrix$tip),
  "regime tree must share the analysis tree's edge matrix" =
    identical(regime_tree$edge, analysis_tree$edge)
)

n_missing_cells <- sum(is.na(tip_matrix[GRAMBANK_FEATURE_IDS]))

write_csv(tip_matrix, here("data", "grammar", "grammar_evolutionary_tip_matrix.csv"))
write.tree(tree_pruned, here("data", "grammar", "grammar_evolutionary_tree.nwk"))
saveRDS(
  list(
    analysis_tree   = analysis_tree,
    regime_tree     = regime_tree,
    philippine_tips = philippine_tips,
    regime_tips     = regime_tips,
    swept_in_tips   = swept_in_tips
  ),
  here("data", "grammar", "grammar_evolutionary_analysis_tree.rds")
)

message(
  "\n=== Final counts ===\n",
  "Pulses tree tips:                      ", Ntip(tree), "\n",
  "Matched to covered GRAMBANK languages: ", nrow(matched), "\n",
  "After collapsing duplicate Glottocodes:", nrow(matched_unique), "\n",
  "Sibling-conflict tips dropped:         ", length(CONFLICTING_TIPS_DROPPED), "\n",
  "\nFinal analysis tree:\n",
  "  ", Ntip(analysis_tree), " tips (", length(philippine_tips), " Philippine, ",
  Ntip(analysis_tree) - length(philippine_tips), " background)\n",
  "  regime clade: ", length(regime_tips), " tips\n",
  "  missing cells: ", n_missing_cells, " of ",
  nrow(tip_matrix) * length(GRAMBANK_FEATURE_IDS),
  " (", round(100 * n_missing_cells / (nrow(tip_matrix) * length(GRAMBANK_FEATURE_IDS)), 2), "%)\n",
  "\nExported for [1]/[2]:\n",
  "  data/grammar/grammar_evolutionary_tip_matrix.csv\n",
  "  data/grammar/grammar_evolutionary_tree.nwk\n",
  "  data/grammar/grammar_evolutionary_analysis_tree.rds"
)
