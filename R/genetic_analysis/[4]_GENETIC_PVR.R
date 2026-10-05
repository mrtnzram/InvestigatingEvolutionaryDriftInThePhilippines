# =============================================================================
# [4] Genetic Analysis — geography vs. pure population structure
#
# Residualizes population-mean PLINK PC1 on distance from Manila to get
# geography-free `pure_phylogeny`, then regresses normalized Spanish admixture on
# geography and on pure structure, with zero-admixture populations removed and kept.
#
# Input:   data/network_distance/GENETIC_final.csv,
#          data/genetic/Phil_2.24M_pca_results.eigenvec (PLINK --pca, individual-level, headerless)
# Outputs: figures/regression/genetic_pvr_{profile,scatter,y_scatter,varpart}[_zeros_kept].png
#          `pvr_results_table`, `pvr_geo_lm`, `df` (with `pure_phylogeny`) in the environment
# =============================================================================

library(tidyverse)
library(here)
source(here("R", "shared", "pure_phylogeny.R"))

PREDICTOR <- "geodist_H2_span"         # terrain-penalized cost to Manila

# Worst qpAdm-SE bucket to drop; matches the ">= 0.15" tier of [7]_GENETIC_geoplots.R's size legend.
SE_HIGH_CUTOFF <- 0.15


# ── Setup: population PC1 and analysis frames ─────────────────────────────
eigenvec_raw <- read.table(here("data", "genetic", "Phil_2.24M_pca_results.eigenvec"),
                           header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eigenvec_raw) - 2
names(eigenvec_raw) <- c("population", "sample", paste0("PC", seq_len(n_pc)))
# Same population-label reconciliation [5]_GENETIC_MMRR.R applies to .mdist.id.
eigenvec_raw$population <- ifelse(eigenvec_raw$population == "ManoboRajahKabunsuwan",
                                  "ManoboRK", eigenvec_raw$population)
eigenvec_pop <- eigenvec_raw |>
  group_by(population) |>
  summarise(PC1 = mean(PC1), .groups = "drop")

GENETIC_final <- read.csv(here("data", "network_distance", "GENETIC_final.csv"))
stopifnot(
  "eigenvec population set != GENETIC_final population set" =
    setequal(eigenvec_pop$population, GENETIC_final$population)
)

# qpAdm source populations carry no estimate and drop out.
df_full <- GENETIC_final |>
  dplyr::select(population, span_admx, span_se_capped, x_km = all_of(PREDICTOR)) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  as.data.frame()

stopifnot("span_se_capped is NA for a population with a span_admx estimate." =
            !anyNA(df_full$span_se_capped))
n_high_se <- sum(df_full$span_se_capped >= SE_HIGH_CUTOFF)
df_full   <- df_full[df_full$span_se_capped < SE_HIGH_CUTOFF, , drop = FALSE]
message(n_high_se, " high-SE populations dropped (span_se_capped >= ", SE_HIGH_CUTOFF,
        ") | ", nrow(df_full), " retained.")

df_full$y <- df_full$span_admx / max(df_full$span_admx)
df        <- df_full[df_full$span_admx > 0, , drop = FALSE]   # primary: non-zero set
message(sum(df_full$span_admx == 0), " zero-admixture populations | ",
        nrow(df), " retained for the primary fit.")


# ── 1–3. Axis 1, orientation, residualization, y regressions ──────────────
genetic_fit <- function(frame, zeros) {
  pure_phylo_fit(frame, frame$PC1, NA_real_, predictor = PREDICTOR, response = "span_admx_norm",
                 axis_name = "PC1", anc_label = "population structure", zeros = zeros)
}
removed <- genetic_fit(df,      "removed")
kept    <- genetic_fit(df_full, "kept")

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
  ggsave(here("figures", "regression", paste0("genetic_pvr_profile", r$sfx, ".png")),
         plot_profile(r$fit, sprintf("Genetic geography and pure structure (zeros %s)", r$tag),
                      "population structure"),
         width = 8, height = 6, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("genetic_pvr_scatter", r$sfx, ".png")),
         plot_scatter(r$fit, sprintf("Genetic structure vs. geography (zeros %s)", r$tag),
                      "population structure", "PC1"),
         width = 10, height = 4.5, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("genetic_pvr_y_scatter", r$sfx, ".png")),
         plot_y_scatter(r$fit, sprintf("Spanish admixture (zeros %s)", r$tag), "population structure",
                        "Normalized Spanish admixture"),
         width = 10, height = 4.5, units = "in", dpi = 300)
  ggsave(here("figures", "regression", paste0("genetic_pvr_varpart", r$sfx, ".png")),
         vp_diagram(r$fit$vp, sprintf("Genetic variance partition (zeros %s)", r$tag),
                    sprintf("r^2 = %.3f", r$fit$r2_full), "population structure"),
         width = 7, height = 3.6, units = "in", dpi = 300)
}


# ── 5. Results table ──────────────────────────────────────────────────────
pvr_results_table <- bind_rows(removed$rows, kept$rows)
print(pvr_results_table)
