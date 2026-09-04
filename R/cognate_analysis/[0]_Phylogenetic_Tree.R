# =============================================================================
# [0] Cognate Analysis — Prune the phylogenetic tree to the study languages
#
# Reads the MCC tree (data/mcc.tree, ABVD / King et al. 2024 — shared with the
# phoneme and grammar analyses) and reduces it to the 108 cognate languages.
#
# Two things differ from the phoneme/grammar versions of this script:
#
#   1. Matching is exact. Tip labels are <Name>_<ABVD language_id> and
#      COGNATE_loans.csv$language_id *is* that ABVD id, so no spelling/dialect
#      lookup table is needed on this side.
#   2. The tree is reduced to one tip per language. The scrape merged every ABVD
#      wordlist sharing a glottocode into a single row, so a language is one
#      observation here; phoneme/grammar instead keep all dialect tips and
#      duplicate the language's data row once per tip.
#
# The tree is Philippine-archipelago/Formosan-focused, so the Sabah/Borneo and
# Sama-Bajau languages have no tip and are dropped: 108 -> 94.
#
# Run order: after [0]_COGNATE_LOANS.R, before [3]_COGNATE_network_distance.R.
#
# Input:   data/mcc.tree, data/cognate/COGNATE_loans.csv, data/shared/subgroup_palette.csv
#          (canonical subgroup -> colour, from R/shared/subgroup_palette.R — run
#          that first)
# Outputs: data/cognate/COGNATE_tree_languages.csv (the 94-language analysis set;
#          [3] restricts to it, so everything below [3] inherits the 94),
#          data/cognate/COGNATE_phylo_dist_matrix.csv (94x94 patristic distances,
#          keyed by glottocode — written for parity with the phoneme/grammar
#          pipelines and for the collinearity check in [4]; there is no cognate
#          [5]_MMRR, since the response is the per-language scalar loans_norm),
#          data/cognate/COGNATE_tip_selection.csv (representative-tip audit),
#          data/cognate/COGNATE_subgroup_lookup.csv (language -> subgroup -> colour),
#          figures/phylogenetic_tree/cognate_phylogenetic_tree.png.
# In-env:  tree_pruned, tip_map — consumed by [4]_COGNATE_PVR.R.
# Next:    [3]_COGNATE_network_distance.R
# =============================================================================

library(tibble)
library(tidyverse)
library(ape)
library(stringr)
library(ggplot2)
library(lingtypology)   # glottolog affiliation -> family subgroup for tip colours
library(here)

COGNATE_loans <- read.csv(here("data", "cognate", "COGNATE_loans.csv"))

tree <- read.nexus(here("data", "shared", "mcc.tree"))


# ── 1. Match cognate languages to tree tips ─────────────────────────────────
# The numeric suffix on each tip label is the ABVD language id, which is the
# exact join key to COGNATE_loans.csv; language_id is comma-joined where several
# wordlists were merged onto one glottocode, so it is split back out to give one
# row per (language, candidate tip).
tree_df <- tibble(original = tree$tip.label) %>%
  mutate(abvd_id = str_extract(original, "[0-9]+$"))

tip_candidates <- COGNATE_loans %>%
  select(glottocode, language, language_id) %>%
  mutate(language_id = as.character(language_id)) %>%
  separate_longer_delim(language_id, ",") %>%
  mutate(language_id = str_trim(language_id)) %>%
  inner_join(tree_df, by = c("language_id" = "abvd_id"))

unmatched <- COGNATE_loans %>% filter(!glottocode %in% tip_candidates$glottocode)
message(
  "Cognate languages with no tree tip (", nrow(unmatched), " of ",
  nrow(COGNATE_loans), "):\n  ", paste0(unmatched$language, collapse = ", ")
)


# ── 2. Reduce to one representative tip per language ────────────────────────
# Medoid selection needs pairwise distances, so the tree is first pruned to every
# candidate tip and only then reduced.
tree_all <- drop.tip(tree, setdiff(tree$tip.label, tip_candidates$original))
D_all    <- cophenetic.phylo(tree_all)

# The representative is the medoid — the tip with the smallest mean patristic
# distance to the others of its glottocode. Most multi-tip groups are clades, but
# a third of them are paraphyletic, so there is no MRCA to collapse to, and
# collapsing one would cost the tree its ultrametricity. drop = FALSE keeps the
# single-tip groups' 1x1 submatrix a matrix, where rowMeans() gives 0.
tip_map <- tip_candidates %>%
  mutate(mean_within = rowMeans(D_all[original, original, drop = FALSE]),
         .by = glottocode) %>%
  slice_min(mean_within, n = 1, with_ties = FALSE, by = glottocode)

tree_pruned <- drop.tip(tree_all, setdiff(tree_all$tip.label, tip_map$original))
# Tips carry the glottocode, not the ABVD name: it is unique across all 108,
# whereas `language` is "; "-joined on merge. [4] joins its data frame on it.
tree_pruned$tip.label <- tip_map$glottocode[match(tree_pruned$tip.label, tip_map$original)]

stopifnot(
  "Tree tips and cognate languages are not 1:1." =
    length(tree_pruned$tip.label) == nrow(tip_map) &&
    setequal(tree_pruned$tip.label, tip_map$glottocode),
  "Reduction left the tree non-ultrametric." = is.ultrametric(tree_pruned)
)

message(
  "\nKept after reduction: ", nrow(tip_map), " of ", nrow(COGNATE_loans),
  " languages (", nrow(tip_candidates), " candidate tips -> ", nrow(tip_map), ")."
)

# The canonical analysis set. [3] reads this to restrict the distance table and
# the waypoint connectors, so everything downstream of it inherits the 94.
tip_map %>%
  select(glottocode, language, tip = original) %>%
  write.csv(here("data", "cognate", "COGNATE_tree_languages.csv"), row.names = FALSE)


# ── 3. What the restriction costs the response ──────────────────────────────
# This is the only script holding both sets, so the loanword accounting for the
# 108 -> 94 cut is reported here.
kept    <- COGNATE_loans %>% filter(glottocode %in% tip_map$glottocode)
dropped <- COGNATE_loans %>% filter(!glottocode %in% tip_map$glottocode)

cat("\n=== Loanwords lost to the tree restriction ===\n")
print(dropped %>%
        filter(number_of_loans > 0) %>%
        arrange(desc(number_of_loans)) %>%
        select(language, number_of_loans, loans_norm))
message(
  sum(dropped$number_of_loans > 0), " of ", nrow(dropped),
  " dropped languages borrow at all; ", sum(dropped$number_of_loans), " of ",
  sum(COGNATE_loans$number_of_loans), " loans (",
  round(100 * sum(dropped$number_of_loans) / sum(COGNATE_loans$number_of_loans), 1),
  "%) leave with them.\nZero-loan share: ",
  sum(COGNATE_loans$number_of_loans == 0), "/", nrow(COGNATE_loans), " -> ",
  sum(kept$number_of_loans == 0), "/", nrow(kept), "."
)

# loans_norm is min-max normalized over all 108 in [0]_COGNATE_LOANS.R. That
# stays valid on the reduced set only while the restriction keeps both extremes;
# if a future ABVD refresh drops the top borrower, [0] has to renormalize.
stopifnot(
  "Restriction changed the loan range — loans_norm must be recomputed on the kept set." =
    all(range(kept$number_of_loans) == range(COGNATE_loans$number_of_loans))
)


# ── 4. Representative-tip audit ─────────────────────────────────────────────
# Which wordlists the merge lumped together, how far apart the tree puts them,
# and which tip was kept. Within-group spread is a check on the upstream
# glottocode merge: a group spanning more than the tree's median root-to-tip
# depth is two near-separate languages sharing one data row.
root_to_tip <- median(diag(vcv(tree_all)))

tip_selection <- tip_candidates %>%
  filter(n() > 1, .by = glottocode) %>%
  group_split(glottocode) %>%
  map_dfr(function(g) {
    sub <- D_all[g$original, g$original]
    kept <- tip_map$original[match(g$glottocode[1], tip_map$glottocode)]
    tibble(
      glottocode     = g$glottocode[1],
      language       = g$language[1],
      n_tips         = nrow(g),
      monophyletic   = is.monophyletic(tree_all, g$original),
      max_within_dist = max(sub),
      chosen_tip     = kept,
      dropped_tips   = paste(setdiff(g$original, kept), collapse = "; ")
    )
  }) %>%
  arrange(desc(max_within_dist))

cat("\n=== Languages merged from several ABVD wordlists ===\n")
print(tip_selection %>% select(-glottocode, -dropped_tips), n = Inf)
message(
  sum(!tip_selection$monophyletic), " of ", nrow(tip_selection),
  " multi-wordlist languages are paraphyletic on the tree."
)

deep_split <- tip_selection %>% filter(max_within_dist > root_to_tip)
if (nrow(deep_split) > 0) {
  message(
    "\nMerged wordlists further apart than the median root-to-tip depth (",
    round(root_to_tip), ") — the tree treats these as separate languages:\n  ",
    paste0(deep_split$language, " (", round(deep_split$max_within_dist), ")",
           collapse = "\n  ")
  )
}

write.csv(tip_selection,
          here("data", "cognate", "COGNATE_tip_selection.csv"), row.names = FALSE)


# ── 5. Coloured phylogram (tips coloured by Glottolog family subgroup) ──────
# The MCC/nexus tree carries no clade annotations, so subgroups come from
# Glottolog. The phoneme and grammar versions of this script reach it through
# ISO 639-3 because that is the key their source tables carry; cognate has the
# glottocode itself, so it joins Glottolog directly — the ISO bridge returns ""
# for some codes, which would then match every iso-less Glottolog row.
# Display labels use the first merged wordlist name plus a "+k" marker; the full
# "; "-joined names run to seven wordlists and would swamp the plot.
tip_subgroup <- tibble(glottocode = tree_pruned$tip.label) %>%
  left_join(tip_map %>% select(glottocode, language), by = "glottocode") %>%
  mutate(
    n_merged = str_count(language, "; ") + 1,
    label = if_else(n_merged > 1,
                    paste0(str_split_i(language, "; ", 1), " (+", n_merged - 1, ")"),
                    language)
  ) %>%
  left_join(
    lingtypology::glottolog %>% select(glottocode, affiliation) %>% distinct(),
    by = "glottocode", relationship = "one-to-one"
  ) %>%
  mutate(subgroup = map_chr(str_split(affiliation, ","), \(x) {
    x <- str_trim(x)
    if (length(x) >= 4) x[4] else tail(x, 1)   # 4th level = finer PH subgroup
  }))

# The layout coordinates read out below are indexed by tip number, so the row
# order here must still be the tree's tip order.
stopifnot(
  "tip_subgroup no longer matches the tree's tips 1:1." =
    identical(tip_subgroup$glottocode, tree_pruned$tip.label)
)

# Palette assigned in order of clade size (largest first).
subgroup_levels <- tip_subgroup %>% count(subgroup, sort = TRUE) %>% pull(subgroup)

# Colours come from the canonical cross-domain palette (R/shared/subgroup_palette.R),
# not a locally-derived one, so the same subgroup is the same colour in every
# domain's figures (phoneme/grammar/cognate/genetic).
subgroup_pal <- read_csv(here("data", "shared", "subgroup_palette.csv"), show_col_types = FALSE) %>%
  select(subgroup, colour) %>%
  deframe()
stopifnot(
  "Some subgroups are missing from data/shared/subgroup_palette.csv — rerun R/shared/subgroup_palette.R." =
    all(subgroup_levels %in% names(subgroup_pal))
)

# Island group, assigned from each language's COORDINATE (not its ethnonym) so
# the split matches the geography every spatial step downstream already uses.
# Keyed by glottocode, not language: cognate's language strings are compound
# doculect labels ("Tagabili; Tboli (Tagabili)") and unstable as a join key.
# MIMAROPA (Mindoro, Romblon, Palawan incl. Cuyo/Cagayancillo) -> Luzon;
# Sulu -> Mindanao.
# Note: this `region` is the Luzon/Visayas/Mindanao sense, i.e. genetics'
# `island_group` column — NOT GENETIC_final.csv's `region` (admin regions).
# Two coordinates that override their names: cebu1242 (8.39, 124.37) is a
# northern-Mindanao ABVD variety, not Cebu; mama1275 (Mamanwa, 9.45, 125.55)
# is Surigao, so Mindanao despite a Visayas-range latitude.
region_lookup <- tribble(
  ~glottocode,  ~region,
  "aben1249",   "Luzon",
  "agut1237",   "Luzon",
  "alab1246",   "Luzon",
  "alan1249",   "Luzon",
  "amba1267",   "Luzon",
  "amga1235",   "Luzon",
  "arta1239",   "Luzon",
  "babu1242",   "Luzon",
  "bala1310",   "Luzon",
  "bant1288",   "Luzon",
  "bata1298",   "Luzon",
  "bata1301",   "Luzon",
  "bino1237",   "Luzon",
  "boli1256",   "Luzon",
  "boto1242",   "Luzon",
  "buhi1243",   "Luzon",
  "buhi1245",   "Luzon",
  "cala1258",   "Luzon",
  "cama1250",   "Luzon",
  "casi1235",   "Luzon",
  "cent2084",   "Luzon",
  "cent2087",   "Luzon",
  "cent2090",   "Luzon",
  "cent2292",   "Luzon",
  "cuyo1237",   "Luzon",
  "dupa1235",   "Luzon",
  "gadd1244",   "Luzon",
  "guin1256",   "Luzon",
  "hanu1241",   "Luzon",
  "ibal1244",   "Luzon",
  "iban1267",   "Luzon",
  "ilok1237",   "Luzon",
  "ilon1239",   "Luzon",
  "iray1237",   "Luzon",
  "irig1242",   "Luzon",
  "isin1239",   "Luzon",
  "isna1241",   "Luzon",
  "itaw1240",   "Luzon",
  "ivat1242",   "Luzon",
  "kaga1256",   "Luzon",
  "kank1243",   "Luzon",
  "kaya1320",   "Luzon",
  "kele1259",   "Luzon",
  "limo1248",   "Luzon",
  "lowe1412",   "Luzon",
  "lubu1243",   "Luzon",
  "maga1263",   "Luzon",
  "magi1241",   "Luzon",
  "nort2877",   "Luzon",
  "pamp1243",   "Luzon",
  "pamp1244",   "Luzon",
  "pang1290",   "Luzon",
  "rata1245",   "Luzon",
  "remo1247",   "Luzon",
  "romb1245",   "Luzon",
  "sout2908",   "Luzon",
  "tady1237",   "Luzon",
  "taga1270",   "Luzon",
  "tina1248",   "Luzon",
  "yoga1237",   "Luzon",
  "atam1240",   "Mindanao",
  "binu1244",   "Mindanao",
  "butu1244",   "Mindanao",
  "cebu1242",   "Mindanao",
  "cent2089",   "Mindanao",
  "cota1241",   "Mindanao",
  "diba1242",   "Mindanao",
  "gian1241",   "Mindanao",
  "ilia1236",   "Mindanao",
  "kala1388",   "Mindanao",
  "kama1363",   "Mindanao",
  "koro1310",   "Mindanao",
  "magu1243",   "Mindanao",
  "mama1275",   "Mindanao",
  "mans1262",   "Mindanao",
  "mara1404",   "Mindanao",
  "mati1250",   "Mindanao",
  "sang1337",   "Mindanao",
  "sara1326",   "Mindanao",
  "sara1327",   "Mindanao",
  "sout2916",   "Mindanao",
  "suri1273",   "Mindanao",
  "taus1251",   "Mindanao",
  "tbol1240",   "Mindanao",
  "tiru1241",   "Mindanao",
  "west2555",   "Mindanao",
  "west2557",   "Mindanao",
  "akla1241",   "Visayas",
  "atii1237",   "Visayas",
  "capi1239",   "Visayas",
  "hili1240",   "Visayas",
  "kina1250",   "Visayas",
  "masb1238",   "Visayas",
  "wara1300",   "Visayas"
)

cognate_lookup_out <- tip_subgroup %>%
  distinct(glottocode, language, subgroup) %>%
  mutate(colour = unname(subgroup_pal[subgroup])) %>%
  left_join(region_lookup, by = "glottocode")

# An exact partition — a typo or a newly-added language must fail here rather
# than silently drop out of [4.3]'s regional split.
stopifnot(
  "region_lookup does not cover every study language." =
    !anyNA(cognate_lookup_out$region),
  "region_lookup names a glottocode that is not in the study set." =
    all(region_lookup$glottocode %in% cognate_lookup_out$glottocode)
)
message("Regions: ", paste(names(table(cognate_lookup_out$region)),
                           table(cognate_lookup_out$region),
                           sep = "=", collapse = ", "))

write.csv(cognate_lookup_out,
          here("data", "cognate", "COGNATE_subgroup_lookup.csv"), row.names = FALSE)

# Let ape compute the rectangular phylogram layout, then read the tip/node
# coordinates back out instead of re-deriving them.
grDevices::pdf(NULL)
plot.phylo(tree_pruned, plot = FALSE)
grDevices::dev.off()
pp     <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
n_tip  <- length(tree_pruned$tip.label)
x_tip  <- max(pp$xx)   # ultrametric tree: every tip sits at this x

edge_df <- tibble(P = tree_pruned$edge[, 1], C = tree_pruned$edge[, 2]) %>%
  mutate(x0 = pp$xx[P], y0 = pp$yy[P], x1 = pp$xx[C], y1 = pp$yy[C])

tip_plot_df <- tip_subgroup %>%
  mutate(
    x = pp$xx[seq_len(n_tip)], y = pp$yy[seq_len(n_tip)],
    subgroup = factor(subgroup, levels = subgroup_levels)
  )

bar_x0 <- x_tip * 1.37
bar_w  <- x_tip * 0.035
grp_x  <- bar_x0 + bar_w + x_tip * 0.02

bar_runs   <- tip_plot_df %>% arrange(y)
run_rle    <- rle(as.character(bar_runs$subgroup))
bar_runs$run <- rep(seq_along(run_rle$lengths), run_rle$lengths)
grp_lab_df <- bar_runs %>%
  group_by(run) %>%
  summarise(subgroup = first(subgroup), y = mean(y), .groups = "drop")

phylo_tree_plot <- ggplot() +
  geom_segment(data = edge_df, aes(x = x0, xend = x1, y = y1, yend = y1),
               linewidth = 0.3, colour = "grey30") +
  geom_segment(data = edge_df, aes(x = x0, xend = x0, y = y0, yend = y1),
               linewidth = 0.3, colour = "grey30") +
  geom_text(data = tip_plot_df, aes(x, y, label = label),
            hjust = 0, nudge_x = x_tip * 0.015, size = 1.5, colour = "grey15") +
  geom_rect(data = tip_plot_df,
            aes(xmin = bar_x0, xmax = bar_x0 + bar_w,
                ymin = y - 0.5, ymax = y + 0.5, fill = subgroup)) +
  geom_text(data = grp_lab_df, aes(x = grp_x, y = y, label = subgroup),
            hjust = 0, size = 2, colour = "grey15") +
  scale_fill_manual(values = subgroup_pal, guide = "none") +
  scale_x_continuous(limits = c(-x_tip * 0.02, x_tip * 1.72), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Philippine study languages (Cognate) - phylogeny (ABVD / King et al. 2024)") +
  theme_void() +
  theme(
    plot.background     = element_rect(fill = "white", colour = NA),
    plot.title.position = "plot",
    plot.title          = element_text(size = 11, hjust = 0.5),
    plot.margin         = margin(6, 6, 6, 6)
  )
print(phylo_tree_plot)

dir.create(here("figures", "phylogenetic_tree"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figures", "phylogenetic_tree", "cognate_phylogenetic_tree.png"),
       phylo_tree_plot, width = 8, height = 18, units = "in", dpi = 300)


# ── 6. Pairwise phylogenetic (patristic) distance matrix ────────────────────
# One tip per language, so this needs no collapsing step: the matrix comes
# straight off the reduced tree, already keyed by glottocode.
COGNATE_phylo_dist_matrix <- cophenetic.phylo(tree_pruned)
diag(COGNATE_phylo_dist_matrix) <- 0

write.csv(COGNATE_phylo_dist_matrix,
          file = here("data", "cognate", "COGNATE_phylo_dist_matrix.csv"),
          row.names = TRUE)
