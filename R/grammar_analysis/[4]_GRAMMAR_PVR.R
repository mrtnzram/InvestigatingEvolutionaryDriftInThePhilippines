# =============================================================================
# [4] Grammar Analysis — geography vs. pure phylogeny
#
# Residualizes the tree's first PCoA axis on distance from Manila to get
# geography-free `pure_phylogeny`, then regresses normalized Spanish grammar
# similarity on geography and on pure phylogeny.
#
# Input:   data/network_distance/GRAMMAR_final.csv; `tree_pruned`, `tree_df_matched` from [0]_Phylogenetic_Tree.R
# Outputs: figures/regression/grammar_pvr_{profile,scatter,y_scatter,varpart}.png
#          `pvr_results_table`, `pvr_geo_lm`, `df` (with `pure_phylogeny`) in the environment
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
tip_map <- tree_df_matched |> dplyr::select(original, gram)

df <- read.csv(here("data", "network_distance", "GRAMMAR_final.csv")) |>
  dplyr::select(language, y = cossim_span_norm, x_km = all_of(PREDICTOR)) |>
  left_join(tip_map, by = c("language" = "gram")) |>
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


# ── 1–3. Axis 1, orientation, residualization, y regressions ──────────────
pcoa <- pcoa_axis1(cophenetic(tree_pruned)[df$original, df$original])
fit  <- pure_phylo_fit(df, pcoa$axis1, pcoa$var_share,
                       predictor = PREDICTOR, response = "cossim_span_norm")
df         <- fit$frame
pvr_geo_lm <- fit$pvr_geo_lm
pvr_varpart <- fit$vp
print(summary(pvr_geo_lm))
print(pvr_varpart)


# ── 4. Plots ──────────────────────────────────────────────────────────────
dir.create(here("figures", "regression"), recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "regression", "grammar_pvr_profile.png"),
       plot_profile(fit, "Grammar geography and pure phylogeny", "phylogeny"),
       width = 8, height = 6, units = "in", dpi = 300)
ggsave(here("figures", "regression", "grammar_pvr_scatter.png"),
       plot_scatter(fit, "Grammar phylogeny vs. geography", "phylogeny", "PCoA axis 1"),
       width = 10, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "regression", "grammar_pvr_y_scatter.png"),
       plot_y_scatter(fit, "Spanish grammar similarity", "phylogeny",
                      "Normalized Spanish grammar similarity"),
       width = 10, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "regression", "grammar_pvr_varpart.png"),
       vp_diagram(pvr_varpart, "Grammar variance partition",
                  sprintf("r^2 = %.3f", fit$r2_full), "phylogeny"),
       width = 7, height = 3.6, units = "in", dpi = 300)


# ── 5. Results table ──────────────────────────────────────────────────────
pvr_results_table <- fit$rows
print(pvr_results_table)
