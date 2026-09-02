# =============================================================================
# [4] Genetic Analysis — geography vs. ancestry (population structure)
#
# Two models for normalized Spanish admixture (y):
#   1. Geography — y ~ terrain-penalized migration distance, uncontrolled.
#   2. Structure-scrubbed geography ("pure ancestry") — y ~ population-structure
#      PC eigenvectors (full block) after migration distance is residualised out
#      of each. Genetics has no tree, so PLINK --pca eigenvectors stand in for
#      the other domains' phylogenetic ones (see [0]_GENETIC_subgroups.R). The
#      scrubbed block is orthogonal to distance, so its R^2 is the semipartial
#      (structure-only) share and geography R^2 + structure R^2 reproduces
#      R^2(y ~ E_sel + distance).
#
# B has no single interpretable slope, so the two are compared as predictors:
# RMSE (and leave-one-out RMSE) of actual vs. predicted y. `cor_with_distance`
# carries the directional read — see §4.
#
# Many populations carry zero Spanish admixture (negatives floored to 0 in
# [0]_GENETIC_ADMX_MALDER.R); that spike at 0 dominates the slope. Both models
# are fit twice — zeros removed (primary) and kept — both conditions in the table.
#
# Input:   data/genetic/network_distance/GENETIC_final.csv (from [3]),
#          data/Phil_2.24M_pca_results.eigenvec (PLINK --pca, individual-level,
#          headerless FID/IID/PC1..; same file [9]_GENETIC_PROCRUSTES.R reads)
# Outputs: data/genetic/PVR/GENETIC_geo_vs_ancestry.csv                                (§4 results table)
#          figures/genetic/regression/genetic_regression_geography[_zeros_kept].png          (y vs. distance)
#          figures/genetic/regression/genetic_regression_actual_vs_predicted[_zeros_kept].png (A | B calibration)
# =============================================================================

library(tidyverse)
library(here)

source(here("R", "shared", "select_moran_eigenvectors.R"))

# H1 = penalized cost onto the nearest network node; H2 = cost to Manila.
PREDICTOR       <- "geodist_H1_span"
DISPLAY_UNIT_KM <- 100                 # slopes are reported per this many km
scl             <- DISPLAY_UNIT_KM / 1000


# ── Setup: population-structure eigenvectors, frames, and the geo/ancestry fitter ──
# PLINK --pca is per individual; aggregate to per-population means.
eigenvec_raw <- read.table(here("data", "Phil_2.24M_pca_results.eigenvec"),
                           header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eigenvec_raw) - 2
names(eigenvec_raw) <- c("population", "sample", paste0("PC", seq_len(n_pc)))
# Same population-label reconciliation [5]_GENETIC_MMRR.R applies to .mdist.id.
eigenvec_raw$population <- ifelse(eigenvec_raw$population == "ManoboRajahKabunsuwan",
                                  "ManoboRK", eigenvec_raw$population)
eigenvec_pop <- eigenvec_raw |>
  group_by(population) |>
  summarise(across(starts_with("PC"), mean), .groups = "drop")

GENETIC_final <- read.csv(here("data", "genetic", "network_distance", "GENETIC_final.csv"))
stopifnot(
  "eigenvec population set != GENETIC_final population set" =
    setequal(eigenvec_pop$population, GENETIC_final$population)
)

# qpAdm source populations carry no estimate and drop out. max(span_admx) is
# unchanged by dropping zeros, so y is on the same scale in both frames.
df_full <- GENETIC_final |>
  dplyr::select(population, span_admx, span_w, x_km = all_of(PREDICTOR)) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  as.data.frame()
df_full$y <- df_full$span_admx / max(df_full$span_admx)
df        <- df_full[df_full$span_admx > 0, , drop = FALSE]   # primary: non-zero set
message(sum(df_full$span_admx == 0), " zero-admixture populations | ",
        nrow(df), " retained for the primary fit.")

# Leave-one-out RMSE via the OLS PRESS shortcut — lets A (1 predictor) and B (k)
# be compared on the same footing.
perf <- function(m) {
  r <- residuals(m); h <- hatvalues(m)
  c(rmse = sqrt(mean(r^2)), rmse_loocv = sqrt(mean((r / (1 - h))^2)))
}

# geo_anc_fit(): steps 1 & 2 for one analysis frame.
#   1. Geography — y ~ migration distance, uncontrolled.
#   2. Structure-scrubbed geography — y ~ PC eigenvectors with migration distance
#      residualised out of each. select_moran_eigenvectors() (R/shared/) is
#      genetics' analogue of PVR(method = "moran"), shared with [4.1].
geo_anc_fit <- function(frame, zeros_label) {
  yy <- frame$y
  xx <- frame$x_km / 1000
  Em <- as.matrix(frame[, paste0("PC", seq_len(n_pc)), drop = FALSE])

  m_geo  <- lm(yy ~ xx)
  gcf    <- summary(m_geo)$coefficients["xx", ]
  geo_r2 <- summary(m_geo)$r.squared

  Es <- select_moran_eigenvectors(yy, xx, Em)
  kk <- ncol(Es)
  m_full <- if (kk == 0) m_geo else lm(yy ~ Es + xx)
  if (kk == 0) {
    m_anc <- lm(yy ~ 1); anc_r2 <- 0; anc_F <- NA_real_; anc_p <- NA_real_; anc_dir <- NA_real_
  } else {
    Ep     <- apply(Es, 2, \(col) resid(lm(col ~ xx)))
    m_anc  <- lm(yy ~ Ep)
    anc_r2 <- summary(m_anc)$r.squared
    inc    <- anova(m_geo, m_full)                 # incremental F: structure on top of geography
    anc_F  <- inc[["F"]][2]; anc_p <- inc[["Pr(>F)"]][2]
    anc_dir <- cor(fitted(lm(yy ~ Es)), xx)        # un-scrubbed structure prediction vs distance
    stopifnot(
      "geography R^2 + structure R^2 != R^2(y ~ E_sel + x)." =
        isTRUE(all.equal(geo_r2 + anc_r2, summary(m_full)$r.squared, tolerance = 1e-6))
    )
  }
  gp <- perf(m_geo); ap <- perf(m_anc)
  list(
    frame = frame, m_geo = m_geo, m_anc = m_anc, zeros = zeros_label,
    rows = tibble(
      model      = c("A_geography", "B_ancestry_given_geography"),
      predictor  = c(paste0(PREDICTOR, " (per ", DISPLAY_UNIT_KM, " km)"),
                     "PC eigenvectors, migration distance scrubbed"),
      zeros = zeros_label, n = nrow(frame),
      n_zeros = sum(frame$span_admx == 0), n_evec = c(0L, kk),
      slope             = c(unname(gcf["Estimate"])   * scl, NA_real_),
      slope_se          = c(unname(gcf["Std. Error"]) * scl, NA_real_),
      statistic         = c(unname(gcf["t value"]), anc_F),
      p.value           = c(unname(gcf["Pr(>|t|)"]), anc_p),
      r.squared         = c(geo_r2, anc_r2),
      cor_with_distance = c(cor(yy, xx), anc_dir),
      rmse              = c(unname(gp["rmse"]),       unname(ap["rmse"])),
      rmse_loocv        = c(unname(gp["rmse_loocv"]), unname(ap["rmse_loocv"]))
    )
  )
}


# ── 1 & 2. Geography and structure-scrubbed geography ─────────────────────
removed <- geo_anc_fit(df,      "removed")
kept    <- geo_anc_fit(df_full, "kept")

for (r in list(removed, kept)) {
  rr <- r$rows
  message(sprintf(
    "[zeros %s] n=%d | R^2 geo = %.3f (RMSE %.4f/LOO %.4f) | structure semipartial R^2 = %.3f (RMSE %.4f/LOO %.4f) | cor(distance): admixture %+.3f, structure %+.3f",
    rr$zeros[1], rr$n[1], rr$r.squared[1], rr$rmse[1], rr$rmse_loocv[1],
    rr$r.squared[2], rr$rmse[2], rr$rmse_loocv[2],
    rr$cor_with_distance[1], rr$cor_with_distance[2]))
}


# ── 3. Plots ──────────────────────────────────────────────────────────────
# Point area = span_w (inverse-variance weight), matching the [7] map.
dir.create(here("figures", "genetic", "regression"), recursive = TRUE, showWarnings = FALSE)

mk_geo_plot <- function(r, tag) {
  rg <- r$rows[r$rows$model == "A_geography", ]
  ggplot(r$frame |> mutate(x_disp = x_km / DISPLAY_UNIT_KM), aes(x_disp, y)) +
    geom_point(aes(size = span_w), alpha = 0.5) +
    scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
    theme_bw() +
    labs(title = sprintf("Spanish admixture vs. migration distance (uncontrolled, zeros %s)", tag),
         subtitle = sprintf("slope = %.4f per %d km (p = %.3g), R^2 = %.3f, RMSE = %.4f, n = %d",
                            rg$slope, DISPLAY_UNIT_KM, rg$p.value, rg$r.squared, rg$rmse, rg$n),
         x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
         y = "Normalized Spanish admixture")
}

mk_avp_plot <- function(r, tag) {
  rg <- r$rows[r$rows$model == "A_geography", ]
  ra <- r$rows[r$rows$model == "B_ancestry_given_geography", ]
  d  <- bind_rows(
    tibble(panel = "geography (A)",                 pred = fitted(r$m_geo), y = r$frame$y, span_w = r$frame$span_w),
    tibble(panel = "structure, geography-free (B)", pred = fitted(r$m_anc), y = r$frame$y, span_w = r$frame$span_w)
  ) |> mutate(panel = factor(panel, c("geography (A)", "structure, geography-free (B)")))
  lab <- tibble(
    panel = factor(c("geography (A)", "structure, geography-free (B)"),
                   c("geography (A)", "structure, geography-free (B)")),
    lab = c(sprintf("RMSE %.4f (LOO %.4f)\nR^2 %.3f", rg$rmse, rg$rmse_loocv, rg$r.squared),
            sprintf("RMSE %.4f (LOO %.4f)\nsemipartial R^2 %.3f", ra$rmse, ra$rmse_loocv, ra$r.squared))
  )
  lim <- range(c(d$pred, d$y))
  ggplot(d, aes(pred, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(aes(size = span_w), alpha = 0.5) +
    scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.15,
              size = 3, lineheight = 0.95, inherit.aes = FALSE) +
    facet_wrap(~ panel) +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_bw() +
    labs(title = sprintf("Predicting Spanish admixture: geography vs. geography-free structure (zeros %s)", tag),
         subtitle = "actual vs. model-predicted; dashed line is y = x",
         x = "Predicted normalized admixture", y = "Actual normalized admixture")
}

ggsave(here("figures", "genetic", "regression", "genetic_regression_geography.png"),
       mk_geo_plot(removed, "removed"), width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "genetic", "regression", "genetic_regression_actual_vs_predicted.png"),
       mk_avp_plot(removed, "removed"), width = 9, height = 4.8, units = "in", dpi = 300)
ggsave(here("figures", "genetic", "regression", "genetic_regression_geography_zeros_kept.png"),
       mk_geo_plot(kept, "kept"), width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "genetic", "regression", "genetic_regression_actual_vs_predicted_zeros_kept.png"),
       mk_avp_plot(kept, "kept"), width = 9, height = 4.8, units = "in", dpi = 300)


# ── 4. Results table ──────────────────────────────────────────────────────
# Rows: model x zero-handling condition. slope: A only (B is a k-vector block).
# statistic/p.value: A = t on the distance coefficient, B = incremental F of the
# structure block on top of geography. r.squared: A = R^2, B = semipartial.
# cor_with_distance: A = raw cor(admixture, distance); B = cor(un-scrubbed
# structure prediction, distance) — sign shows whether structure runs with (-)
# or against (+) the geographic gradient.
GENETIC_geo_vs_ancestry <- bind_rows(removed$rows, kept$rows)
print(GENETIC_geo_vs_ancestry)
write.csv(GENETIC_geo_vs_ancestry,
          file = here("data", "genetic", "PVR", "GENETIC_geo_vs_ancestry.csv"),
          row.names = FALSE)
