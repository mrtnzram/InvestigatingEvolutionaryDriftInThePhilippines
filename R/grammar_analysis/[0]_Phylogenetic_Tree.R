# =============================================================================
# [0] Grammar Analysis — Prune the phylogenetic tree to the study languages
#
# Reads the MCC tree (data/mcc.tree, ABVD / King et al. 2024 — shared with the
# phoneme analysis), reconciles its tip labels with the GRAMBANK Philippine
# language names, and prunes it to the grammar Philippine study set.
#
# The tree is Philippine-archipelago/Formosan-focused: several grammar Sabah/
# Sama-Bajau languages (Balangingi, Bonggi, Pangutaran Sama, Rungus, Sabah
# Malay, Southern Sama, Tambunan Dusun, Tatana, Timugon Murut, West Coast
# Bajau, Yakan) plus Karao and Tagabawa have no tree tip at all and are
# dropped here; everything else matches exactly or via a documented
# spelling/dialect-level correspondence in `lookup` below.
#
# Input:   data/mcc.tree, data/GRAMBANKdf_full.csv (for the Philippine
#          language list and ISO codes, via [0]_GRAMBANKdatabase.R)
# Outputs: data/GRAMMAR_phylo_dist_matrix.csv (patristic distance matrix,
#          consumed by [5]_GRAMMAR_MMRR.R), data/GRAMMAR_subgroup_lookup.csv
#          (language -> subgroup -> colour, consumed by [8]), and
#          figures/shared/grammar_phylogenetic_tree.png.
# In-env:  tree_pruned, tree_df_matched — consumed by [4]_GRAMMAR_PGLS.R.
# Next:    [1]_GRAMMAR_cosine_similarity.R
# =============================================================================

library(tibble)
library(tidyverse)
library(ape)
library(stringr)
library(ggplot2)
library(lingtypology)   # glottolog affiliation -> family subgroup for tip colours
library(here)

GRAMBANKdf_full <- read_csv(here("data", "GRAMBANKdf_full.csv"), show_col_types = FALSE)

Ph_Languages <- GRAMBANKdf_full %>%
  filter(Language_Type == "Philippine Language") %>%
  pull(language)

tree <- read.nexus(here('data', 'mcc.tree'))

tree_df <- tibble(
  original = tree$tip.label
) %>%
  mutate(
    sanitized = str_remove(original, "_[0-9]+$") %>%
      str_replace_all("_", " ") %>%
      str_squish()
  )

# Fuzzy/dialect-level correspondences (spelling variants or a representative
# dialect chosen for a multi-dialect language). Exact-name matches (e.g.
# "Ibaloi", "Pangasinan") need no entry — they fall through the coalesce()
# fallback below.
lookup <- tribble(
  ~tree, ~gram,

  "Agutyanen",                    "Agutaynen",
  "Manobo Ata down-river",        "Ata Manobo",
  "Bontok Guina-ang",             "Central Bontoc",
  "Agta",                         "Central Cagayan Agta",
  "Naga Bikol",                   "Coastal-Naga Bikol",
  "Tagalog",                      "Filipino",
  "Babuyan",                      "Ibatan",
  "Ilongot Kakidugen",            "Ilongot",
  "Isinay Dupax",                 "Isinai",
  "Isneg Dibagat-Kabugao-Isneg",  "Isnag",
  "Kallahan Keleyqiq",            "Keley-i Kallahan",
  "Maguindanaon",                 "Maguindanao",
  "Subanun Sindangan",            "Northern Subanen",
  "Tausug Jolo Dialect",          "Tausug",
  "Tboli Tagabili",               "Tboli",
  "Ifugao Batad",                 "Tuwali Ifugao",
)

tree_df_matched <- tree_df %>%
  left_join(lookup, by = c("sanitized" = "tree")) %>%
  mutate(gram = coalesce(gram, sanitized))

unmatched <- setdiff(Ph_Languages, tree_df_matched$gram)
message(
  "Philippine languages with no tree tip (", length(unmatched), " of ",
  length(Ph_Languages), "):\n  ", paste0(unmatched, collapse = ", ")
)

# --- Pruning -----

tips_keep <- tree_df_matched %>%
  filter(gram %in% Ph_Languages) %>%
  pull(original)

tree_pruned <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_keep)
)

tip_df <- tree_df_matched %>%
  select(original, gram)

Ph_Languages_pruned <- tree_df_matched %>%
  filter(gram %in% Ph_Languages) %>%
  pull(gram)

message(
  "\nKept after pruning: ", length(unique(Ph_Languages_pruned)), " of ",
  length(Ph_Languages), " Philippine languages (", length(tips_keep),
  " tree tips, some languages represented by >1 dialect tip)."
)

# ── Coloured phylogram (tips coloured by Glottolog family subgroup) ──────────
# The MCC/nexus tree carries no clade annotations, so subgroup membership is
# looked up from Glottolog via each language's ISO 639-3 code (lingtypology).
# GRAMBANKdf_full has no iso6393 column (grammar keys on glottocode), so
# bridge Language_ID -> ISO via iso.gltc(), same mechanism used in
# [0]_GRAMBANKdatabase.R to bridge the unrelated control set.
tip_subgroup <- tibble(original = tree_pruned$tip.label) %>%
  left_join(tip_df, by = "original") %>%
  left_join(
    GRAMBANKdf_full %>%
      distinct(gram = language, Language_ID) %>%
      mutate(iso6393 = iso.gltc(Language_ID)),
    by = "gram"
  ) %>%
  left_join(
    lingtypology::glottolog %>% select(iso, affiliation) %>% distinct(),
    by = c("iso6393" = "iso")
  ) %>%
  mutate(subgroup = map_chr(str_split(affiliation, ","), \(x) {
    x <- str_trim(x)
    if (length(x) >= 4) x[4] else tail(x, 1)   # 4th level = finer PH subgroup
  }))

# Palette assigned in order of clade size (largest first).
subgroup_levels <- tip_subgroup %>% count(subgroup, sort = TRUE) %>% pull(subgroup)

soften <- function(cols, s_mult = 0.55, v_mult = 0.95) {
  h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  grDevices::hsv(h["h", ], h["s", ] * s_mult, pmin(h["v", ] * v_mult, 1))
}
.poly <- grDevices::palette.colors(NULL, "Polychrome 36")
.lum  <- colSums(grDevices::col2rgb(.poly) * c(0.299, 0.587, 0.114))  # 0..255
subgroup_pal <- setNames(
  soften(.poly[.lum < 200])[seq_along(subgroup_levels)], subgroup_levels
)

# Export the per-language subgroup -> colour lookup so other scripts colour the
# same languages with the identical palette (e.g. [8]'s tree-vs-network figure
# needs the map points to match these tip colours exactly).
tip_subgroup %>%
  distinct(language = gram, subgroup) %>%
  mutate(colour = unname(subgroup_pal[subgroup])) %>%
  write.csv(here("data", "GRAMMAR_subgroup_lookup.csv"), row.names = FALSE)

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
  geom_text(data = tip_plot_df, aes(x, y, label = gram),
            hjust = 0, nudge_x = x_tip * 0.015, size = 2, colour = "grey15") +
  geom_rect(data = tip_plot_df,
            aes(xmin = bar_x0, xmax = bar_x0 + bar_w,
                ymin = y - 0.5, ymax = y + 0.5, fill = subgroup)) +
  geom_text(data = grp_lab_df, aes(x = grp_x, y = y, label = subgroup),
            hjust = 0, size = 2.3, colour = "grey15") +
  scale_fill_manual(values = subgroup_pal, guide = "none") +
  scale_x_continuous(limits = c(-x_tip * 0.02, x_tip * 1.72), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Philippine study languages (Grammar) - phylogeny (ABVD / King et al. 2024)") +
  theme_void() +
  theme(
    plot.background     = element_rect(fill = "white", colour = NA),
    plot.title.position = "plot",
    plot.title          = element_text(size = 11, hjust = 0.5),
    plot.margin         = margin(6, 6, 6, 6)
  )
print(phylo_tree_plot)

dir.create(here("figures", "shared"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figures", "shared", "grammar_phylogenetic_tree.png"),
       phylo_tree_plot, width = 7.5, height = 9, units = "in", dpi = 300)

# ── Pairwise phylogenetic (patristic) distance matrix ───────────────────────
# Some study languages are represented by more than one tree tip; collapsed by
# averaging every original-tip pair's distance within each gram-to-gram pair,
# matching GRAMMAR_dist_matrix.csv / GRAMMAR_diss_matrix.csv from [5].
phylo_dist_raw <- cophenetic.phylo(tree_pruned)

phylo_dist_long <- as_tibble(phylo_dist_raw, rownames = "original_1") %>%
  pivot_longer(-original_1, names_to = "original_2", values_to = "phylo_dist") %>%
  left_join(tip_df, by = c("original_1" = "original")) %>%
  rename(gram_1 = gram) %>%
  left_join(tip_df, by = c("original_2" = "original")) %>%
  rename(gram_2 = gram) %>%
  filter(gram_1 != gram_2) %>%
  summarise(phylo_dist = mean(phylo_dist), .by = c(gram_1, gram_2))

GRAMMAR_phylo_dist_matrix <- phylo_dist_long %>%
  pivot_wider(names_from = gram_2, values_from = phylo_dist) %>%
  column_to_rownames("gram_1") %>%
  as.matrix()

GRAMMAR_phylo_dist_matrix <- GRAMMAR_phylo_dist_matrix[, rownames(GRAMMAR_phylo_dist_matrix)]

diag(GRAMMAR_phylo_dist_matrix) <- 0

write.csv(GRAMMAR_phylo_dist_matrix,
          file = here("data", "GRAMMAR_phylo_dist_matrix.csv"), row.names = TRUE)
