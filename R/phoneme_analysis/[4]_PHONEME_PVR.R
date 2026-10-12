# =============================================================================
# [4] Phoneme Analysis — geography vs. pure phylogeny
#
# Residualizes the tree's first PCoA axis on a Manila-anchored signed geographic
# axis and on distance from Manila to get `pure_phylogeny`, then regresses
# normalized Spanish phoneme similarity on distance from Manila and on pure phylogeny.
#
# Input:   data/network_distance/PHONEME_final.csv, data/network_distance/PHONEME_dist_matrix.csv;
#          `tree_pruned`, `tree_df_matched` from [0]_Phylogenetic_Tree.R
# Outputs: figures/regression/pvr/phoneme_pvr_{profile,scatter,y_scatter,y_geography,y_pure,varpart}.png
#          `pvr_results_table`, `pvr_checks`, `pvr_geo_lm`, `df` (with `pure_phylogeny`) in the environment
# =============================================================================

library(ape)
library(tidyverse)
library(here)
source(here("R", "shared", "pure_phylogeny.R"))

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

PREDICTOR <- "geodist_H2_span"         # terrain-penalized cost to Manila


# ── Setup: analysis frame ─────────────────────────────────────────────────
# tip_map duplicates multi-tip languages: one row per tip, values repeated.
tip_map <- tree_df_matched |> dplyr::select(original, ph)

df <- read.csv(here("data", "network_distance", "PHONEME_final.csv")) |>
  dplyr::select(language, lat = latitude, y = cossim_span_norm, x_km = all_of(PREDICTOR)) |>
  left_join(tip_map, by = c("language" = "ph")) |>
  filter(!is.na(original), !is.na(y), !is.na(x_km)) |>
  as.data.frame()
df <- df[match(tree_pruned$tip.label, df$original), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the tip set." =
    nrow(df) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order." =
    !anyNA(df$original) && identical(df$original, tree_pruned$tip.label)
)
message(nrow(df), " tips covering ", n_distinct(df$language), " languages.")


G_geo <- as.matrix(read.csv(here("data", "network_distance", "PHONEME_dist_matrix.csv"),
                            row.names = 1, check.names = FALSE))


# ── 1–3. Axes, orientation, residualization, y regressions ────────────────
geo  <- geo_axis_frame(df, "language", G_geo)
pcoa <- pcoa_axis1(cophenetic(tree_pruned)[df$original, df$original])
fit  <- pure_phylo_fit(df, pcoa$axis1, geo, pcoa$var_share,
                       predictor = PREDICTOR, response = "cossim_span_norm")
df         <- fit$frame
pvr_geo_lm <- fit$pvr_geo_lm
pvr_varpart <- fit$vp
pvr_checks <- fit$checks
print(summary(pvr_geo_lm))
print(pvr_varpart)
print(pvr_checks, n = Inf)


# ── 4. Plots ──────────────────────────────────────────────────────────────
dir.create(here("figures", "regression", "pvr"), recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "regression", "pvr", "phoneme_pvr_profile.png"),
       plot_profile(fit, "Phoneme geography and pure phylogeny", "phylogeny"),
       width = 8, height = 6, units = "in", dpi = 300)
ggsave(here("figures", "regression", "pvr", "phoneme_pvr_scatter.png"),
       plot_scatter(fit, "Phoneme phylogeny vs. geography", "phylogeny", "PCoA axis 1"),
       width = 10, height = 4.5, units = "in", dpi = 300)
save_y_scatters(fit, here("figures", "regression", "pvr"), "phoneme_pvr",
                "Spanish phoneme similarity", "phylogeny", "Normalized Spanish phoneme similarity")
ggsave(here("figures", "regression", "pvr", "phoneme_pvr_varpart.png"),
       vp_diagram(pvr_varpart, "Phoneme variance partition",
                  sprintf("r^2 = %.3f", fit$r2_full), "phylogeny"),
       width = 7, height = 3.6, units = "in", dpi = 300)


# ── 5. Results table ──────────────────────────────────────────────────────
pvr_results_table <- fit$rows
print(pvr_results_table)
