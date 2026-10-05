# =============================================================================
# [4] Cognate Analysis — geography vs. pure phylogeny
#
# Residualizes the tree's first PCoA axis on distance from Manila to get
# geography-free `pure_phylogeny`, then regresses the normalized Spanish loanword
# count on geography and on pure phylogeny, with zero-loan languages removed and kept.
#
# Input:   data/network_distance/COGNATE_final.csv; `tree_pruned`, `tip_map` from [0]_Phylogenetic_Tree.R
# Outputs: figures/regression/cognate_pvr_{profile,scatter,y_scatter,varpart}[_zeros_kept].png
#          `pvr_results_table`, `pvr_geo_lm`, `df` (with `pure_phylogeny`) in the environment
# =============================================================================

library(ape)
library(tidyverse)
library(here)
source(here("R", "shared", "pure_phylogeny.R"))

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." = exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tip_map` is not defined."     = exists("tip_map")
)

PREDICTOR <- "geodist_H2_span"         # terrain-penalized cost to Manila


# ── Setup: analysis frames ────────────────────────────────────────────────
# [0] relabelled tips to glottocodes, one per language — no duplication here.
df_full <- read.csv(here("data", "network_distance", "COGNATE_final.csv")) |>
  dplyr::select(glottocode, language, y = loans_norm, x_km = all_of(PREDICTOR)) |>
  semi_join(tip_map, by = "glottocode") |>
  filter(!is.na(y), !is.na(x_km)) |>
  as.data.frame()
df <- df_full[df_full$y > 0, , drop = FALSE]                 # primary: non-zero set

D_phylo <- cophenetic(tree_pruned)
stopifnot("Analysis languages missing from the pruned tree." =
            all(df_full$glottocode %in% rownames(D_phylo)))
message(sum(df_full$y == 0), " zero-loan languages | ", nrow(df), " retained for the primary fit.")


# ── 1–3. Axis 1, orientation, residualization, y regressions ──────────────
# Axis 1 is recomputed within each frame so pure_phylogeny is orthogonal to geography in each fit.
cognate_fit <- function(frame, zeros) {
  pcoa <- pcoa_axis1(D_phylo[frame$glottocode, frame$glottocode])
  pure_phylo_fit(frame, pcoa$axis1, pcoa$var_share,
                 predictor = PREDICTOR, response = "loans_norm", zeros = zeros)
}
removed <- cognate_fit(df,      "removed")
kept    <- cognate_fit(df_full, "kept")

df         <- removed$frame
pvr_geo_lm <- removed$pvr_geo_lm
pvr_varpart <- bind_rows(
  mutate(removed$vp, zeros = "removed", .before = 1),
  mutate(kept$vp,    zeros = "kept",    .before = 1)
)
print(summary(pvr_geo_lm))
print(pvr_varpart)


# ── 4. Plots ──────────────────────────────────────────────────────────────
dir.create(here("figures", "regression"), recursive = TRUE, showWarnings = FALSE)

for (r in list(list(fit = removed, tag = "removed", sfx = ""),
               list(fit = kept,    tag = "kept",    sfx = "_zeros_kept"))) {
  ggsave(here("figures", "regression", paste0("cognate_pvr_profile", r$sfx, ".png")),
         plot_profile(r$fit, sprintf("Cognate geography and pure phylogeny (zeros %s)", r$tag),
                      "phylogeny"),
         width = 8, height = 6, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("cognate_pvr_scatter", r$sfx, ".png")),
         plot_scatter(r$fit, sprintf("Cognate phylogeny vs. geography (zeros %s)", r$tag),
                      "phylogeny", "PCoA axis 1"),
         width = 10, height = 4.5, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("cognate_pvr_y_scatter", r$sfx, ".png")),
         plot_y_scatter(r$fit, sprintf("Spanish loanwords (zeros %s)", r$tag), "phylogeny",
                        "Normalized Spanish loanword count"),
         width = 10, height = 4.5, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("cognate_pvr_varpart", r$sfx, ".png")),
         vp_diagram(r$fit$vp, sprintf("Cognate variance partition (zeros %s)", r$tag),
                    sprintf("r^2 = %.3f", r$fit$r2_full), "phylogeny"),
         width = 7, height = 3.6, units = "in", dpi = 300)
}


# ── 5. Results table ──────────────────────────────────────────────────────
pvr_results_table <- bind_rows(removed$rows, kept$rows)
print(pvr_results_table)
