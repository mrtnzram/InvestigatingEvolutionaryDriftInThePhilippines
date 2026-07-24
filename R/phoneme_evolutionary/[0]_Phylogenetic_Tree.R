# =============================================================================
# [0] Phoneme Evolutionary — Match the pulses tree to the Ruhlen phoneme data
#
# Reads the pulses tree (data/pulses_summary.trees, 400 tips), and joins its
# tips to the Ruhlen/Creanza phoneme languages via Glottocode: Ruhlen's
# iso6393 -> Glottocode (lingtypology::gltc.iso()), then Glottocode is the
# shared key with data/pulses_languages.csv (whose Name column supplies the
# tree's own tip labels). Unlike the old mcc.tree script, no hand-built
# tip-name lookup table is needed — the pulses tree tips already correspond
# almost 1:1 to Ruhlen's own language names (e.g. "AttaPamplona",
# "DumagatCasiguran"), and the join is done at the Glottocode level instead.
#
# DEPENDENCY / RUN ORDER:
#   - REQUIRES `Ph_Languages` and `ruhlen_raw` from PART A of
#     [0]_CREANZA_RUHLENdatabase.R (run that first, down to its END OF PART A
#     banner).
#   - PRODUCES `Ph_Languages_pruned`, consumed by PART B of
#     [0]_CREANZA_RUHLENdatabase.R.
#   - PRODUCES `tree_pruned` (Austronesian-wide, tips kept under the pulses
#     tree's OWN original labels — not renamed to Ruhlen language names, so
#     tip labels stay guaranteed-unique) and a simple ape::plot.phylo() sanity
#     plot, plus a second plot of just the Philippine-language clade.
#   - WRITES data/phoneme_evolutionary_tree.nwk (Newick, tree_pruned) and
#     data/phoneme_evolutionary_tip_mapping.csv (one row per tip: tip ->
#     Ruhlen language -> iso6393/glottocode/inventory_size), consumed by
#     [1]_PHONEME_evolutionary_ML_model.R.
#
# Data quirks handled here:
#   - `Baliledo` in pulses_languages.csv ships with a blank Glottocode. It is
#     manually assigned "bali1288" (Glottolog's "Baliledu", a Central-East
#     Sumbanese dialect with no ISO code of its own — confirmed by the tree
#     placing Baliledo next to Pondok in the Sumba clade). Baliledo has no
#     corresponding entry anywhere in the Ruhlen dataset, so this fix does not
#     make it match a phoneme record — it is still dropped when the tree is
#     pruned to matched tips, same as every other unmatched tip. The manual
#     Glottocode is assigned purely so that omission is visible/intentional
#     rather than a silent NA.
#   - A handful of Glottocodes are shared by two distinct (non-identical)
#     Ruhlen phoneme records under the same ISO code (e.g. ivat1242 covers
#     both "Itbayaten" and "Ivatanen" — 8 of 728 phonemes differ between
#     them). Tips are kept under their own original pulses label (so this
#     never produces a duplicate tip label), but each such tip's Ruhlen match
#     is ambiguous between the two records until step 4b resolves it: keep
#     whichever record has the larger phoneme inventory (Kapuas/Katingan:
#     24 vs. 19, clean). Itbayaten/Ivatanen tie exactly at 26 — that specific
#     tie is broken manually in Itbayaten's favor, not derived from the data.
# =============================================================================

library(tibble)
library(tidyverse)
library(ape)
library(stringr)
library(lingtypology)   # gltc.iso() bridges Ruhlen's iso6393 -> Glottocode
library(here)

stopifnot(
  "Run PART A of [0]_CREANZA_RUHLENdatabase.R first: `Ph_Languages` is not defined." =
    exists("Ph_Languages"),
  "Run PART A of [0]_CREANZA_RUHLENdatabase.R first: `ruhlen_raw` is not defined." =
    exists("ruhlen_raw")
)

tree <- read.nexus(here("data", "pulses_summary.trees"))

# ---- 1. Ruhlen (Austronesian) -> Glottocode -------------------------------
ruhlen_austronesian_gltc <- ruhlen_raw |>
  filter(language_family == "Austronesian") |>
  mutate(
    glottocode      = gltc.iso(iso6393),
    inventory_size  = rowSums(across(all_of(phoneme_cols)))
  )

unresolved <- ruhlen_austronesian_gltc |> filter(is.na(glottocode))
if (nrow(unresolved) > 0) {
  message(
    nrow(unresolved), " Austronesian Ruhlen language(s) have no resolvable Glottocode ",
    "(dropped from the tree match):\n",
    paste0("  ", unresolved$language, " (iso6393 = ", unresolved$iso6393, ")", collapse = "\n")
  )
}

# ---- 2. pulses_languages.csv, with the Baliledo Glottocode patched --------
pulses_languages <- read_csv(here("data", "pulses_languages.csv"), show_col_types = FALSE) |>
  mutate(Glottocode = if_else(Name == "Baliledo", "bali1288", Glottocode))

# ---- 3. Mapping table: Glottocode join between Ruhlen and the pulses tree -
mapping_table <- ruhlen_austronesian_gltc |>
  filter(!is.na(glottocode)) |>
  inner_join(
    pulses_languages |> dplyr::select(Name, Glottocode),
    by = c("glottocode" = "Glottocode"),
    relationship = "many-to-many"
  )

# ---- 4. Diagnostics: matched / dropped / duplicated ------------------------
dropped_tips <- setdiff(tree$tip.label, mapping_table$Name)
matched_tips <- unique(mapping_table$Name)

message(
  "Baliledo dropped as expected (no Ruhlen match despite manual Glottocode fix): ",
  "Baliledo" %in% dropped_tips
)

dup_language <- mapping_table |> count(language, sort = TRUE) |> filter(n > 1)
message(
  nrow(dup_language), " Ruhlen language(s) matched to more than one tip ",
  "(dialect-level fan-out on the tree side):\n",
  paste0("  ", dup_language$language, " (", dup_language$n, " tips)", collapse = "\n")
)

# Tips are kept under their own original pulses label (no renaming — see
# header), so duplicate tip labels are impossible by construction. What CAN
# happen is the reverse: one tip's Glottocode matches more than one distinct
# Ruhlen phoneme record, leaving that tip's phoneme match ambiguous.
ambiguous_tips <- mapping_table |>
  distinct(Name, language) |>
  count(Name, sort = TRUE) |>
  filter(n > 1)
message(
  nrow(ambiguous_tips), " tip(s) matched more than one Ruhlen language under the same ",
  "Glottocode (phoneme match ambiguous between them):\n",
  paste0(
    "  ", ambiguous_tips$Name, " -> ",
    map_chr(ambiguous_tips$Name, \(nm) {
      mapping_table |> filter(Name == nm) |> pull(language) |> unique() |> sort() |> paste(collapse = " / ")
    }),
    collapse = "\n"
  )
)

# ---- 4b. Resolve ambiguous tips: keep the larger phoneme inventory ---------
# Each ambiguous tip's candidate languages are ranked by inventory_size and
# the top one kept. Itbayaten/Ivatanen tie exactly (26 phonemes each) — the
# general rule can't break that one, so it's resolved manually in Itbayaten's
# favor (documented decision, not derived from the data).
mapping_resolved <- mapping_table |>
  group_by(Name) |>
  slice_max(order_by = inventory_size, n = 1, with_ties = TRUE) |>
  filter(!(n() > 1 & language == "Ivatanen")) |>
  ungroup()

resolved_log <- mapping_table |>
  semi_join(ambiguous_tips, by = "Name") |>
  distinct(Name, language, inventory_size) |>
  left_join(
    mapping_resolved |> distinct(Name, language) |> mutate(kept = TRUE),
    by = c("Name", "language")
  ) |>
  mutate(kept = coalesce(kept, FALSE)) |>
  arrange(Name, desc(inventory_size))
message(
  "\nAmbiguous tips resolved by phoneme inventory size (ties broken manually):\n",
  paste0(
    "  ", resolved_log$Name, ": ", resolved_log$language, " (", resolved_log$inventory_size, ")",
    if_else(resolved_log$kept, " <- kept", ""),
    collapse = "\n"
  )
)

# ---- 5. Prune (original tip labels kept as-is) ------------------------------
tree_pruned <- drop.tip(tree, dropped_tips)

# `Ph_Languages_pruned` — the Philippine-language subset validated against the
# pulses tree, consumed by PART B of [0]_CREANZA_RUHLENdatabase.R.
Ph_Languages_pruned <- mapping_resolved |>
  filter(language %in% Ph_Languages) |>
  pull(language) |>
  unique()

# Original pulses tip labels (not Ruhlen language names) that carry at least
# one matched Philippine language — used below for the clade plot.
ph_tip_names <- mapping_resolved |>
  filter(language %in% Ph_Languages &
        !language %in% c('Yakan','Timugon','Sangil')) |>
  pull(Name) |>
  unique()

message(
  "\nPhilippine languages matched to the pulses tree: ", length(Ph_Languages_pruned),
  " of ", length(Ph_Languages), " detected by the bounding box.\n",
  "  Dropped: ", paste0(setdiff(Ph_Languages, Ph_Languages_pruned), collapse = ", ")
)

# ---- Final counts ------------------------------------------------------------
message(
  "\n=== Final counts ===\n",
  "Total tips (pulses tree):              ", length(tree$tip.label), "\n",
  "Dropped tips:                          ", length(dropped_tips), "\n",
  "Total Austronesian languages (Ruhlen): ", nrow(ruhlen_austronesian_gltc), "\n",
  "Dropped Austronesian languages:        ", nrow(ruhlen_austronesian_gltc) - n_distinct(mapping_resolved$language), "\n",
  "  (of which lost the inventory-size tie-break rather than lacking a tip: ",
  n_distinct(mapping_table$language) - n_distinct(mapping_resolved$language), ")\n",
  "Total Philippine languages:            ", length(Ph_Languages), "\n",
  "Dropped Philippine languages:          ", length(Ph_Languages) - length(Ph_Languages_pruned)
)

# ---- 6. Simple sanity plot ---------------------------------------------------
plot.phylo(tree_pruned, cex = 0.4, no.margin = TRUE)

# ---- 7. Philippine-language clade plot --------------------------------------
# Smallest monophyletic subtree containing every matched Philippine language
# (matched on the tree's own original tip labels via `ph_tip_names`).
ph_mrca  <- getMRCA(tree_pruned, ph_tip_names)
ph_clade <- extract.clade(tree_pruned, ph_mrca)

# The Philippine languages need not be monophyletic in this tree, so the
# smallest common ancestor can sweep in non-Philippine tips too.
non_ph_in_clade <- setdiff(ph_clade$tip.label, ph_tip_names)
message(
  "Philippine-language clade: ", length(ph_clade$tip.label), " tips (",
  length(ph_tip_names), " Philippine, ", length(non_ph_in_clade),
  " non-Philippine swept in by the smallest common ancestor)."
)
if (length(non_ph_in_clade) > 0) {
  message("  Non-Philippine tips in the clade: ", paste(non_ph_in_clade, collapse = ", "))
}

plot.phylo(ph_clade, cex = 0.6, no.margin = TRUE)

# ---- 8. Export tree + tip mapping for [1]_PHONEME_evolutionary_ML_model.R --
# tree_pruned's tips are still the pulses tree's own original labels (step 5).
# The ambiguity resolved in step 4b means this is now a clean one-to-one
# lookup: exactly one Ruhlen phoneme record per tip.
tip_mapping <- mapping_resolved |>
  dplyr::select(tip = Name, ruhlen_language = language, iso6393, glottocode, inventory_size) |>
  arrange(tip)

stopifnot("tip_mapping should be one row per tip after ambiguity resolution" =
            !anyDuplicated(tip_mapping$tip))

write_csv(tip_mapping, here("data", "phoneme_evolutionary_tip_mapping.csv"))
write.tree(tree_pruned, here("data", "phoneme_evolutionary_tree.nwk"))

message(
  "\nExported for [1]: ", length(tree_pruned$tip.label), "-tip tree -> ",
  "data/phoneme_evolutionary_tree.nwk\n",
  nrow(tip_mapping), "-row tip mapping (one row per tip, ambiguity resolved) -> ",
  "data/phoneme_evolutionary_tip_mapping.csv"
)
