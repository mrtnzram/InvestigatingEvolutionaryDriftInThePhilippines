# =============================================================================
# [4] Cognate Analysis — geography vs. ancestry
#
# Two models for the normalized Spanish loanword count (y):
#   1. Geography — y ~ terrain-penalized migration distance, uncontrolled.
#   2. Phylogeny-scrubbed geography ("pure ancestry") — y ~ phylogenetic
#      eigenvectors (full block) after migration distance is residualised out of
#      each. The scrubbed block is orthogonal to distance, so its R^2 is the
#      semipartial (ancestry-only) share and geography R^2 + ancestry R^2
#      reproduces R^2(y ~ E_sel + distance).
#
# B has no single interpretable slope, so the two are compared as predictors:
# RMSE (and leave-one-out RMSE) of actual vs. predicted y. `cor_with_distance`
# carries the directional read — see §4.
#
# Most languages have zero confirmed Spanish loans; that spike at 0 dominates
# the slope. Both models are fit twice — with the zero-loan languages removed
# (primary) and kept — and both conditions are in the results table.
#
# Run order: requires `tree_pruned` and `tip_map` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/network_distance/COGNATE_final.csv (from [3]_COGNATE_network_distance.R)
# Outputs: data/pvr/COGNATE_geo_vs_ancestry.csv                                (§4 results table)
#          figures/regression/cognate_regression_geography[_zeros_kept].png          (y vs. distance)
#          figures/regression/cognate_regression_actual_vs_predicted[_zeros_kept].png (A | B calibration)
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." = exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tip_map` is not defined."     = exists("tip_map")
)

# H1 = penalized cost onto the nearest network node; H2 = cost to Manila.
PREDICTOR       <- "geodist_H1_span"
DISPLAY_UNIT_KM <- 100                 # slopes are reported per this many km
scl             <- DISPLAY_UNIT_KM / 1000


# ── Setup: analysis frames, trees, and the geo/ancestry fitter ─────────────
# [0] relabelled tips to glottocodes, one per language — no duplication here.
df_full <- read.csv(here("data", "network_distance", "COGNATE_final.csv")) |>
  dplyr::select(glottocode, language, y = loans_norm, x_km = all_of(PREDICTOR)) |>
  semi_join(tip_map, by = "glottocode") |>
  filter(!is.na(y), !is.na(x_km)) |>
  as.data.frame()

df        <- df_full[df_full$y > 0, , drop = FALSE]          # primary: non-zero set
tree_a    <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, df$glottocode))
tree_full <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, df_full$glottocode))
message(sum(df_full$y == 0), " zero-loan languages | ", nrow(df), " retained for the primary fit.")

# Leave-one-out RMSE via the OLS PRESS shortcut — lets A (1 predictor) and B (k)
# be compared on the same footing.
perf <- function(m) {
  r <- residuals(m); h <- hatvalues(m)
  c(rmse = sqrt(mean(r^2)), rmse_loocv = sqrt(mean((r / (1 - h))^2)))
}

# geo_anc_fit(): steps 1 & 2 for one analysis frame + its pruned tree.
#   1. Geography — y ~ migration distance, uncontrolled.
#   2. Phylogeny-scrubbed geography — y ~ phylogenetic eigenvectors with
#      migration distance residualised out of each ("moran" selection). The
#      block is orthogonal to distance, so its R^2 is the semipartial share.
geo_anc_fit <- function(frame, tree_sub, zeros_label) {
  frame <- frame[match(tree_sub$tip.label, frame$glottocode), ]
  stopifnot(!anyNA(frame$glottocode), identical(frame$glottocode, tree_sub$tip.label))
  yy <- frame$y
  xx <- frame$x_km / 1000

  m_geo  <- lm(yy ~ xx)
  gcf    <- summary(m_geo)$coefficients["xx", ]
  geo_r2 <- summary(m_geo)$r.squared

  pf  <- PVR(PVRdecomp(tree_sub, scale = TRUE), phy = tree_sub,
             trait = yy, envVar = xx, method = "moran")
  Es  <- as.matrix(pf@Selection$Vectors)
  kk  <- ncol(Es)
  m_full <- if (kk == 0) m_geo else lm(yy ~ Es + xx)
  if (kk == 0) {
    m_anc <- lm(yy ~ 1); anc_r2 <- 0; anc_F <- NA_real_; anc_p <- NA_real_; anc_dir <- NA_real_
  } else {
    Ep     <- apply(Es, 2, \(col) resid(lm(col ~ xx)))
    m_anc  <- lm(yy ~ Ep)
    anc_r2 <- summary(m_anc)$r.squared
    inc    <- anova(m_geo, m_full)                 # incremental F: ancestry on top of geography
    anc_F  <- inc[["F"]][2]; anc_p <- inc[["Pr(>F)"]][2]
    anc_dir <- cor(fitted(lm(yy ~ Es)), xx)        # un-scrubbed ancestry prediction vs distance
    stopifnot(
      "geography R^2 + ancestry R^2 != R^2(y ~ E_sel + x)." =
        isTRUE(all.equal(geo_r2 + anc_r2, summary(m_full)$r.squared, tolerance = 1e-6))
    )
  }
  gp <- perf(m_geo); ap <- perf(m_anc)
  list(
    frame = frame, m_geo = m_geo, m_anc = m_anc, zeros = zeros_label,
    rows = tibble(
      model      = c("A_geography", "B_ancestry_given_geography"),
      predictor  = c(paste0(PREDICTOR, " (per ", DISPLAY_UNIT_KM, " km)"),
                     "phylo eigenvectors, migration distance scrubbed"),
      zeros = zeros_label, n = nrow(frame), n_zeros = sum(frame$y == 0), n_evec = c(0L, kk),
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


# ── 1 & 2. Geography and phylogeny-scrubbed geography ─────────────────────
removed <- geo_anc_fit(df,      tree_a,    "removed")
kept    <- geo_anc_fit(df_full, tree_full, "kept")

for (r in list(removed, kept)) {
  rr <- r$rows
  message(sprintf(
    "[zeros %s] n=%d | R^2 geo = %.3f (RMSE %.4f/LOO %.4f) | ancestry semipartial R^2 = %.3f (RMSE %.4f/LOO %.4f) | cor(distance): metric %+.3f, ancestry %+.3f",
    rr$zeros[1], rr$n[1], rr$r.squared[1], rr$rmse[1], rr$rmse_loocv[1],
    rr$r.squared[2], rr$rmse[2], rr$rmse_loocv[2],
    rr$cor_with_distance[1], rr$cor_with_distance[2]))
}


# ── 3. Plots ──────────────────────────────────────────────────────────────
dir.create(here("figures", "regression"), recursive = TRUE, showWarnings = FALSE)

mk_geo_plot <- function(r, tag) {
  rg <- r$rows[r$rows$model == "A_geography", ]
  ggplot(r$frame |> mutate(x_disp = x_km / DISPLAY_UNIT_KM), aes(x_disp, y)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
    theme_bw() +
    labs(title = sprintf("Spanish loanwords vs. migration distance (uncontrolled, zeros %s)", tag),
         subtitle = sprintf("slope = %.4f per %d km (p = %.3g), R^2 = %.3f, RMSE = %.4f, n = %d",
                            rg$slope, DISPLAY_UNIT_KM, rg$p.value, rg$r.squared, rg$rmse, rg$n),
         x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
         y = "Normalized Spanish loanword count")
}

mk_avp_plot <- function(r, tag) {
  rg <- r$rows[r$rows$model == "A_geography", ]
  ra <- r$rows[r$rows$model == "B_ancestry_given_geography", ]
  d  <- bind_rows(
    tibble(panel = "geography (A)",               pred = fitted(r$m_geo), y = r$frame$y),
    tibble(panel = "ancestry, geography-free (B)", pred = fitted(r$m_anc), y = r$frame$y)
  ) |> mutate(panel = factor(panel, c("geography (A)", "ancestry, geography-free (B)")))
  lab <- tibble(
    panel = factor(c("geography (A)", "ancestry, geography-free (B)"),
                   c("geography (A)", "ancestry, geography-free (B)")),
    lab = c(sprintf("RMSE %.4f (LOO %.4f)\nR^2 %.3f", rg$rmse, rg$rmse_loocv, rg$r.squared),
            sprintf("RMSE %.4f (LOO %.4f)\nsemipartial R^2 %.3f", ra$rmse, ra$rmse_loocv, ra$r.squared))
  )
  lim <- range(c(d$pred, d$y))
  ggplot(d, aes(pred, y)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 2, alpha = 0.5) +
    geom_text(data = lab, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.15,
              size = 3, lineheight = 0.95, inherit.aes = FALSE) +
    facet_wrap(~ panel) +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_bw() +
    labs(title = sprintf("Predicting Spanish loanwords: geography vs. geography-free ancestry (zeros %s)", tag),
         subtitle = "actual vs. model-predicted; dashed line is y = x",
         x = "Predicted normalized loanword count", y = "Actual normalized loanword count")
}

ggsave(here("figures", "regression", "cognate_regression_geography.png"),
       mk_geo_plot(removed, "removed"), width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "regression", "cognate_regression_actual_vs_predicted.png"),
       mk_avp_plot(removed, "removed"), width = 9, height = 4.8, units = "in", dpi = 300)
ggsave(here("figures", "regression", "cognate_regression_geography_zeros_kept.png"),
       mk_geo_plot(kept, "kept"), width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "regression", "cognate_regression_actual_vs_predicted_zeros_kept.png"),
       mk_avp_plot(kept, "kept"), width = 9, height = 4.8, units = "in", dpi = 300)


# ── 4. Results table ──────────────────────────────────────────────────────
# Rows: model x zero-handling condition. slope: A only (B is a k-vector block).
# statistic/p.value: A = t on the distance coefficient, B = incremental F of the
# ancestry block on top of geography. r.squared: A = R^2, B = semipartial.
# cor_with_distance: A = raw cor(metric, distance); B = cor(un-scrubbed ancestry
# prediction, distance) — sign shows whether ancestry runs with (-) or against
# (+) the geographic gradient.
COGNATE_geo_vs_ancestry <- bind_rows(removed$rows, kept$rows)
print(COGNATE_geo_vs_ancestry)
write.csv(COGNATE_geo_vs_ancestry,
          file = here("data", "pvr", "COGNATE_geo_vs_ancestry.csv"),
          row.names = FALSE)
