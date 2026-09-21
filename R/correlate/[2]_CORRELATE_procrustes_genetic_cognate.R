# =============================================================================
# [2] Correlate Analysis — Procrustes: genetic PCA vs. cognate PCA
#
# Where each domain's [9] script fits its PC1/PC2 onto geography, this fits the
# two PCAs onto each other: genetic population structure vs. lexical (cognate)
# structure, on the languages both datasets cover. Pairs by glottocode, bridged
# from genetic populations through GENETIC_subgroup_lookup.csv (the same bridge
# as [0]_CORRELATE_admixture.R). Populations sharing a glottocode are averaged
# to one row, the way genetic [9] averages individuals to populations.
#
# The PCs are reused from the [9] scripts, not recomputed on the intersection:
# the genetic PCA comes from a precomputed PLINK eigenvec, and there is no local
# genotype matrix to rerun it on.
#
# Run order: needs [9]_COGNATE_PROCRUSTES.R and [9]_GENETIC_PROCRUSTES.R to have
# been run (or R/shared/[9]_ALL_PROCRUSTES.R).
#
# Input:   data/procrustes/GENETIC_procrustes_scores.csv
#          data/procrustes/COGNATE_procrustes_scores.csv
#          data/genetic/GENETIC_subgroup_lookup.csv
#          data/shared/subgroup_palette.csv
# Outputs: data/correlate/CORRELATE_procrustes_genetic_cognate.csv
#          data/correlate/CORRELATE_procrustes_genetic_cognate_scores.csv
#          figures/correlate/procrustes_genetic_cognate.png
# =============================================================================

library(adegenet)   # loads ade4 (procuste, procuste.randtest)
library(tidyverse)
library(here)
library(patchwork)

dir.create(here("data", "correlate"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("figures", "correlate"), recursive = TRUE, showWarnings = FALSE)

N_PERM <- 999   # + observed arrangement = 1000 draws

# ── 1. Pair the two PCAs by glottocode ───────────────────────────────────────
GENETIC_lookup <- read_csv(here("data", "genetic", "GENETIC_subgroup_lookup.csv"),
                           show_col_types = FALSE)
GENETIC_scores <- read_csv(here("data", "procrustes", "GENETIC_procrustes_scores.csv"),
                           show_col_types = FALSE)
COGNATE_scores <- read_csv(here("data", "procrustes", "COGNATE_procrustes_scores.csv"),
                           show_col_types = FALSE)

genetic_glc <- GENETIC_scores |>
  left_join(GENETIC_lookup |> select(population, glottocode), by = "population") |>
  filter(glottocode %in% COGNATE_scores$glottocode) |>
  summarise(gen_PC1 = mean(PC1), gen_PC2 = mean(PC2),
            populations = paste(sort(population), collapse = "; "),
            n_populations = n(), .by = glottocode)

subgroup_palette <- read_csv(here("data", "shared", "subgroup_palette.csv"),
                             show_col_types = FALSE)

analysis_df <- COGNATE_scores |>
  select(glottocode, language, subgroup, cog_PC1 = PC1, cog_PC2 = PC2) |>
  left_join(subgroup_palette |> select(subgroup, colour), by = "subgroup") |>
  inner_join(genetic_glc, by = "glottocode") |>
  arrange(glottocode)

stopifnot(
  "Duplicate glottocodes after averaging populations." = !anyDuplicated(analysis_df$glottocode),
  "Some paired rows are missing a subgroup colour." = !anyNA(analysis_df$colour)
)
N <- nrow(analysis_df)
message(N, " glottocodes shared (", sum(analysis_df$n_populations), " genetic populations).")

# ── 2. Procrustes fit + permutation test ─────────────────────────────────────
gen_mat <- analysis_df |> select(PC1 = gen_PC1, PC2 = gen_PC2) |> as.data.frame()
cog_mat <- analysis_df |> select(PC1 = cog_PC1, PC2 = cog_PC2) |> as.data.frame()
rownames(gen_mat) <- rownames(cog_mat) <- analysis_df$glottocode

pr  <- procuste(gen_mat, cog_mat, scale = TRUE, nf = 2)
ss  <- sum((as.matrix(pr$rotX) - as.matrix(pr$tabY))^2)   # Gower's M^2, dimensionless
prt <- procuste.randtest(gen_mat, cog_mat, nrepet = N_PERM)
correlation <- prt$obs
p_value     <- prt$pvalue

# ── 3. Results + scores tables ───────────────────────────────────────────────
CORRELATE_procrustes <- tibble(n = N, ss_gower = ss, correlation = correlation,
                               p_value = p_value, n_perm = N_PERM)
print(CORRELATE_procrustes)
write_csv(CORRELATE_procrustes,
          here("data", "correlate", "CORRELATE_procrustes_genetic_cognate.csv"))

# rotX = genetic rotated onto the cognate configuration; tabY = the cognate
# configuration, centred and scaled — the two sets plotted against each other.
analysis_df <- analysis_df |>
  mutate(gen_rot1 = pr$rotX[, 1], gen_rot2 = pr$rotX[, 2],
         cog_std1 = pr$tabY[, 1], cog_std2 = pr$tabY[, 2],
         residual = sqrt((gen_rot1 - cog_std1)^2 + (gen_rot2 - cog_std2)^2))
write_csv(analysis_df,
          here("data", "correlate", "CORRELATE_procrustes_genetic_cognate_scores.csv"))

# ── 4. Figure: cognate and rotated genetic side by side, shared axes ──────
subgroup_pal <- analysis_df |> distinct(subgroup, colour) |> deframe()

# One set of limits for both panels, so a language's position reads the same
# way left and right.
lims <- range(c(analysis_df$cog_std1, analysis_df$cog_std2,
                analysis_df$gen_rot1, analysis_df$gen_rot2))

panel <- function(x, y, title) {
  ggplot(analysis_df, aes(.data[[x]], .data[[y]], fill = subgroup)) +
    geom_point(shape = 21, colour = "black", size = 3, stroke = 0.4) +
    scale_fill_manual(values = subgroup_pal) +
    coord_fixed(xlim = lims, ylim = lims) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw() +
    theme(legend.position = "none")
}

p_proc <- panel("cog_std1", "cog_std2", "Cognate") +
  panel("gen_rot1", "gen_rot2", "Genetic (Procrustes-rotated)") +
  plot_annotation(title = "Genetic vs cognate Procrustes",
                  subtitle = sprintf("p = %.3f, r^2 = %.3f", p_value, correlation^2))

print(p_proc)
ggsave(here("figures", "correlate", "procrustes_genetic_cognate.png"), p_proc,
       width = 12, height = 6.5, units = "in", dpi = 300, bg = "white")
