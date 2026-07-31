# =============================================================================
# [0] Phoneme Analysis — Prune the phylogenetic tree to the study languages
#
# Reads the MCC tree (data/mcc.tree), reconciles its tip labels with the
# Ruhlen/phoneme language names, and prunes it to the Philippine study set.
#
# Run order: requires `Ph_Languages` from PART A of [0]_CREANZA_RUHLENdatabase.R,
# and produces `Ph_Languages_pruned` for its PART B.
#
# Input:   data/mcc.tree
# Outputs: data/PHONEME_phylo_dist_matrix.csv (patristic distances, the
#          phylogenetic predictor in [5]), data/PHONEME_subgroup_lookup.csv
#          (language -> subgroup -> colour, consumed by [8]), and
#          figures/shared/phylogenetic_tree.png.
# In-env:  tree_pruned, tree_df_matched — consumed by [4]_PHONEME_PGLS.R.
# =============================================================================

library(tibble)
library(tidyverse)
library(ape)
library(stringr)
library(ggplot2)
library(lingtypology)   # glottolog affiliation -> family subgroup for tip colours
library(here)

stopifnot(
  "Run PART A of [0]_CREANZA_RUHLENdatabase.R first: `Ph_Languages` is not defined." =
    exists("Ph_Languages")
)

tree <- read.nexus(here('data','mcc.tree'))

tree_df <- tibble(
  original = tree$tip.label
) %>%
  mutate(
    sanitized = str_remove(original, "_[0-9]+$") %>%
      str_replace_all("_", " ") %>%
      str_squish()
  )

lookup <- tribble(
  ~tree, ~ph,
  
  # exact matches
  "Agta", "Agta",
  "Gaddang", "Gaddang",
  "Ibanag", "Ibanag",
  "Ilokano", "Ilokano",
  "Balangaw", "Balangaw",
  "Inibaloi", "Inibaloi",
  "Iraya", "Iraya",
  "Binukid", "Binukid",
  "Mamanwa", "Mamanwa",
  "Hanunoo", "Hanunoo",
  "Buhid", "Buhid",
  "Tagalog", "Tagalog",
  "Kalagan", "Kalagan",
  "Cebuano", "Cebuano",
  "Hiligaynon", "Hiligaynon",
  "Maranao", "Maranao",
  "Tiruray", "Tiruray",
  "Pangasinan", "Pangasinan",
  "Kapampangan", "Kapampangan",
  "Yogad", "Yogad",
  
  # variants
  "Maguindanaon", "Magindanao",
  "Cuyonon", "Kuyunon",
  "Mansaka", "Mansakan",
  "Sangil Saragani Islands", "Sangil",
  "Tausug Jolo Dialect", "Tausug",
  "Samar-Leyte", "Waray",
  
  # Agta / Dumagat
  "Atta Pamplona", "Atta",
  "Dumagat Casiguran", "Casiguran Dumagat",
  
  # Isneg
  "Isneg Dibagat-Kabugao-Isneg", "Isnag",
  
  # Itneg
  "Itneg Binongan", "Itneg",
  
  # Kalinga
  "Kalinga Guinaang Lubuagan Dialect", "North Kalinga",
  "Kalinga Minangali", "Central Kalinga",
  "Kalinga Southern", "South Kalinga",
  
  # Ifugao
  "Ifugao Amganad", "Central Ifugao",
  "Ifugao Batad", "East Ifugao",
  
  # Bontok
  "Bontok Guina-ang", "Central Bontok",
  
  # Kankanaey
  "Kankanay Northern", "North Kankanay",
  
  # Kallahan
  "Kallahan Kayapa Proper", "North Kallahan",
  
  # Ilongot
  "Ilongot Kakidugen", "Ilongot",
  
  # Ivatan
  "Ivatan Basco Dialect", "Ivatanen",
  
  # Sambal
  "Sambal", "Tina",
  
  # Manobo
  "Manobo Ata up-river", "Ata",
  "Manobo Ata down-river", "Ata",
  "Manobo Dibabawon", "Dibabawon",
  "Manobo Ilianen Kibudtungan Dialect", "Ilianen",
  "Manobo Tigwa Iglogsad Dialect", "Tigwa",
  "Manobo Western Bukidnon", "Western Bukidnon",
  "Manobo Sarangani Kayaponga Dialect", "Sarangani Manobo",
  "Manobo Kalamansig Cotabato Paril Dialect", "Kalamansig",
  
  # Bilaan
  "Bilaan Koronadal", "Cotabato Bilaan",
  "Bilaan Sarangani", "Sarangani Bilaan",
  
  # Batak / Tagbanwa
  "Palawan Batak", "Batak",
  "Batak Palawan", "Batak",
  "Tagbanwa Kalamian Coron Island Dialect", "Kalamianen",
  "Tagbanwa Aborlan Dialect", "Tagbanwa",
  
  # Bikol
  "Naga Bikol", "Naga",
  
  # Subanon
  "Subanon Siocon", "Siocon Subanon",
  "Subanun Sindangan", "Sindangan",
  
  # Tboli
  "Tboli Tagabili", "Tboli"
)

Ph_Languages <- setdiff(Ph_Languages, "Ga-dang")

tree_df_matched <- tree_df %>%
  left_join(lookup, by = c("sanitized" = "tree")) %>%
  mutate(ph = coalesce(ph, sanitized))

setdiff(Ph_Languages, tree_df_matched$ph)

# --- Pruning -----

tips_keep <- tree_df_matched %>%
  filter(ph %in% Ph_Languages) %>%
  pull(original)

tree_pruned <- drop.tip(
  tree,
  setdiff(tree$tip.label, tips_keep)
)

tip_df <- tree_df_matched %>%
  select(original, ph)

Ph_Languages_pruned <- tree_df_matched %>%
  filter(ph %in% Ph_Languages) %>%
  pull(ph)


# ── Coloured phylogram (tips coloured by Glottolog family subgroup) ──────────
# The MCC/nexus tree carries no clade annotations, so subgroups come from
# Glottolog via each language's ISO 639-3 code. The affiliation path's 4th level
# is the Philippine subgroup; shorter paths fall back to their deepest level.
# Drawn from ape's own layout coordinates to avoid a ggtree dependency.
tip_subgroup <- tibble(original = tree_pruned$tip.label) %>%
  left_join(tip_df, by = "original") %>%
  left_join(
    read_csv(here("data", "RUHLENdf_PH.csv"), show_col_types = FALSE) %>%
      select(ph = language, iso6393),
    by = "ph"
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

# Palette for the level-4 subgroups: Polychrome 36 minus its lightest colours
# (invisible on white), desaturated, and named by subgroup so the mapping is stable.
soften <- function(cols, s_mult = 0.55, v_mult = 0.95) {
  h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  grDevices::hsv(h["h", ], h["s", ] * s_mult, pmin(h["v", ] * v_mult, 1))
}
.poly <- grDevices::palette.colors(NULL, "Polychrome 36")
.lum  <- colSums(grDevices::col2rgb(.poly) * c(0.299, 0.587, 0.114))  # 0..255
subgroup_pal <- setNames(
  soften(.poly[.lum < 200])[seq_along(subgroup_levels)], subgroup_levels
)

# Export the subgroup -> colour lookup so [8]'s map points match these tip
# colours exactly; distinct() collapses tip_subgroup's dialect duplicates.
tip_subgroup %>%
  distinct(language = ph, subgroup) %>%
  mutate(colour = unname(subgroup_pal[subgroup])) %>%
  write.csv(here("data", "PHONEME_subgroup_lookup.csv"), row.names = FALSE)

# Let ape compute the phylogram layout and read its coordinates back out:
# plot = FALSE fills .PlotPhyloEnv without drawing, absorbed by the null device.
grDevices::pdf(NULL)
plot.phylo(tree_pruned, plot = FALSE)
grDevices::dev.off()
pp     <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
n_tip  <- length(tree_pruned$tip.label)
x_tip  <- max(pp$xx)   # ultrametric tree: every tip sits at this x

# Two segments per edge draw the elbow: a horizontal branch at the child's y and
# a vertical connector at the parent's x spanning the parent->child y gap.
edge_df <- tibble(P = tree_pruned$edge[, 1], C = tree_pruned$edge[, 2]) %>%
  mutate(x0 = pp$xx[P], y0 = pp$yy[P], x1 = pp$xx[C], y1 = pp$yy[C])

tip_plot_df <- tip_subgroup %>%
  mutate(
    x = pp$xx[seq_len(n_tip)], y = pp$yy[seq_len(n_tip)],
    subgroup = factor(subgroup, levels = subgroup_levels)
  )

# Group annotation strip replaces the legend: a colour bar beside the tips,
# labelled once per contiguous block, so non-monophyletic subgroups get one
# label per block rather than a single ambiguous entry.
bar_x0 <- x_tip * 1.37          # strip sits clear of the left-aligned tip labels
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
  geom_text(data = tip_plot_df, aes(x, y, label = ph),
            hjust = 0, nudge_x = x_tip * 0.015, size = 2, colour = "grey15") +
  # colour strip: one tile per tip
  geom_rect(data = tip_plot_df,
            aes(xmin = bar_x0, xmax = bar_x0 + bar_w,
                ymin = y - 0.5, ymax = y + 0.5, fill = subgroup)) +
  # Label text is dark rather than the group colour: it abuts its own tile, so
  # the tie is already spatial and dark type stays legible for pale subgroups.
  geom_text(data = grp_lab_df, aes(x = grp_x, y = y, label = subgroup),
            hjust = 0, size = 2.3, colour = "grey15") +
  scale_fill_manual(values = subgroup_pal, guide = "none") +
  # Right limit crops the panel just past the subgroup labels; clip = "off" lets
  # an overhanging glyph draw into the margin instead of being cut.
  scale_x_continuous(limits = c(-x_tip * 0.02, x_tip * 1.72), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Philippine study languages - phylogeny (ABVD / King et al. 2024)") +
  theme_void() +
  theme(
    # theme_void() leaves the background transparent; set it white so text is
    # legible and the PNG isn't see-through.
    plot.background     = element_rect(fill = "white", colour = NA),
    # centre the title over the whole plot width, not just the panel.
    plot.title.position = "plot",
    plot.title          = element_text(size = 11, hjust = 0.5),
    plot.margin         = margin(6, 6, 6, 6)
  )
print(phylo_tree_plot)

dir.create(here("figures", "shared"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figures", "shared", "phylogenetic_tree.png"),
       phylo_tree_plot, width = 7.5, height = 9, units = "in", dpi = 300)


# ── Pairwise phylogenetic (patristic) distance matrix ───────────────────────
# Tree-only artifact, so it is built here rather than in [4]. Some study
# languages have more than one tree tip; those are collapsed by averaging every
# original-tip pair within each ph-to-ph pair, giving one row/col per language
# to match the [5] matrices.
phylo_dist_raw <- cophenetic.phylo(tree_pruned)

phylo_dist_long <- as_tibble(phylo_dist_raw, rownames = "original_1") %>%
  pivot_longer(-original_1, names_to = "original_2", values_to = "phylo_dist") %>%
  left_join(tip_df, by = c("original_1" = "original")) %>%
  rename(ph_1 = ph) %>%
  left_join(tip_df, by = c("original_2" = "original")) %>%
  rename(ph_2 = ph) %>%
  filter(ph_1 != ph_2) %>%
  summarise(phylo_dist = mean(phylo_dist), .by = c(ph_1, ph_2))

PHONEME_phylo_dist_matrix <- phylo_dist_long %>%
  pivot_wider(names_from = ph_2, values_from = phylo_dist) %>%
  column_to_rownames("ph_1") %>%
  as.matrix()

# pivot_wider orders columns by first appearance, not by row order, so reindex
# to make diag() address true self-pairs.
PHONEME_phylo_dist_matrix <- PHONEME_phylo_dist_matrix[, rownames(PHONEME_phylo_dist_matrix)]

# ph_1 != ph_2 above drops the diagonal (self-pairs) along with same-ph dialect
# pairs; 0 is the conventional value and matches [5]'s dist/diss matrices.
diag(PHONEME_phylo_dist_matrix) <- 0

write.csv(PHONEME_phylo_dist_matrix,
          file = here("data", "PHONEME_phylo_dist_matrix.csv"), row.names = TRUE)

