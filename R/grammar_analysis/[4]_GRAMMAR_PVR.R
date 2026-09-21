# =============================================================================
# [4] Grammar Analysis — geography vs. ancestry
#
# Two models for normalized Spanish-grammar cosine similarity (y):
#   1. Geography — y ~ terrain-penalized migration distance, uncontrolled.
#   2. Phylogeny-scrubbed geography ("pure ancestry") — y ~ phylogenetic
#      eigenvectors (full block) after migration distance is residualised out of
#      each. The scrubbed block is orthogonal to distance, so its R^2 is the
#      semipartial (ancestry-only) share and geography R^2 + ancestry R^2
#      reproduces R^2(y ~ E_sel + distance).
#
# B has no single interpretable slope, so the two are compared as predictors:
# RMSE (and leave-one-out RMSE) of actual vs. predicted y. `cor_with_distance`
# carries the directional read — see §5.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/network_distance/GRAMMAR_final.csv (from [3]_GRAMMAR_network_distance.R)
# Outputs: figures/regression/grammar_regression_geography.png           (y vs. distance)
#          figures/regression/grammar_regression_actual_vs_predicted.png (A | B calibration)
#          figures/regression/grammar_pvr_varpart.png                    (a/b/c/d partition)
#
# The §5 results table is left in the environment as `GRAMMAR_geo_vs_ancestry`
# rather than written here. R/shared/[4]_ALL_PVR.R harvests it and writes
# data/pvr/geo_vs_ancestry.csv, consolidated across all four domains.
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

# H1 = penalized cost onto the nearest network node; H2 = cost to Manila.
PREDICTOR       <- "geodist_H1_span"
DISPLAY_UNIT_KM <- 100                 # slopes are reported per this many km


# ── Setup: analysis frame + phylogenetic eigenvectors ──────────────────────
# PVR matches trait to eigenvector rows by POSITION — row i is
# tree_pruned$tip.label[i]. Reorder and assert; a mis-sort fits silently.
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
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(df$original) && identical(df$original, tree_pruned$tip.label)
)

y   <- df$y
x   <- df$x_km / 1000                  # fit on Mm
N   <- nrow(df)
scl <- DISPLAY_UNIT_KM / 1000          # Mm -> display units
# Cosine similarity is strictly positive — no zero spike to drop (cognate/genetic do).
n_zeros <- 0L

message(N, " tips covering ", n_distinct(df$language), " languages ",
        "(multi-dialect languages contribute one row per tip).")

# Phylogenetic eigenvectors: "moran" adds vectors until residual autocorrelation
# on the tree is n.s. (sig.t = 0.05); connectivity left at the package default.
pvr_fit <- PVR(PVRdecomp(tree_pruned, scale = TRUE), phy = tree_pruned,
               trait = y, envVar = x, method = "moran")
E_sel <- as.matrix(pvr_fit@Selection$Vectors)
k     <- ncol(E_sel)
m_full <- if (k == 0) lm(y ~ x) else lm(y ~ E_sel + x)   # reference for the §2 checks

# Leave-one-out RMSE via the OLS PRESS shortcut — lets A (1 predictor) and B (k)
# be compared on the same footing.
perf <- function(m) {
  r <- residuals(m); h <- hatvalues(m)
  c(rmse = sqrt(mean(r^2)), rmse_loocv = sqrt(mean((r / (1 - h))^2)))
}


# ── 1. Geography ──────────────────────────────────────────────────────────
m_geo    <- lm(y ~ x)
geo_cf   <- summary(m_geo)$coefficients["x", ]
geo_r2   <- summary(m_geo)$r.squared
geo_perf <- perf(m_geo)
geo_dir  <- cor(y, x)                  # signed effect direction (magnitude = sqrt(geo_r2))


# ── 2. Phylogeny-scrubbed geography (pure ancestry) ───────────────────────
if (k == 0) {
  m_anc <- lm(y ~ 1); anc_r2 <- 0; anc_F <- NA_real_; anc_p <- NA_real_; anc_dir <- NA_real_
} else {
  E_perp <- apply(E_sel, 2, \(col) resid(lm(col ~ x)))   # distance scrubbed from each
  m_anc  <- lm(y ~ E_perp)
  anc_r2 <- summary(m_anc)$r.squared                     # semipartial (ancestry-only) R^2
  fst    <- summary(m_anc)$fstatistic                    # overall F-test of B: does ancestry itself predict y?
  anc_F  <- unname(fst[1])
  anc_p  <- unname(pf(fst[1], fst[2], fst[3], lower.tail = FALSE))
  # Direction: correlate ancestry's UN-scrubbed predicted metric with distance.
  # Sign says whether the phylogenetic signal, on its own, runs with (-) or
  # against (+) the geographic gradient. After scrubbing it is ~0 by construction.
  anc_dir <- cor(fitted(lm(y ~ E_sel)), x)
  stopifnot(
    "Scrubbed-eigenvector coefficients drift from the reference PVR fit." =
      isTRUE(all.equal(unname(coef(m_anc)[-1]),
                       unname(coef(m_full)[2:(k + 1)]), tolerance = 1e-8)),
    "geography R^2 + ancestry R^2 != R^2(y ~ E_sel + x)." =
      isTRUE(all.equal(geo_r2 + anc_r2, summary(m_full)$r.squared, tolerance = 1e-6))
  )
}
anc_perf <- perf(m_anc)

message(sprintf("R^2  geography = %.3f | ancestry|geography (semipartial) = %.3f | sum = %.3f",
                geo_r2, anc_r2, geo_r2 + anc_r2))
message(sprintf("RMSE geography = %.4f (LOO %.4f) | ancestry = %.4f (LOO %.4f)",
                geo_perf["rmse"], geo_perf["rmse_loocv"], anc_perf["rmse"], anc_perf["rmse_loocv"]))
message(sprintf("cor with migration distance:  metric = %+.3f | raw ancestry prediction = %+.3f",
                geo_dir, anc_dir))


# ── 3. Variance partition ─────────────────────────────────────────────────
# Nested-R^2 split of the full model. b is a difference of R^2, not a variance,
# and goes negative under suppression.
R2_full <- summary(m_full)$r.squared
R2_anc  <- if (k == 0) 0 else summary(lm(y ~ E_sel))$r.squared

GRAMMAR_pvr_varpart <- tibble(
  component = c("a", "b", "c", "d"),
  meaning   = c("geography only", "shared", "ancestry only", "unexplained"),
  fraction  = c(R2_full - R2_anc, geo_r2 + R2_anc - R2_full, anc_r2, 1 - R2_full)
)
print(GRAMMAR_pvr_varpart)
stopifnot("Variance partition does not sum to 1." =
            isTRUE(all.equal(sum(GRAMMAR_pvr_varpart$fraction), 1, tolerance = 1e-8)))


# ── 4. Plots ──────────────────────────────────────────────────────────────
dir.create(here("figures", "regression"), recursive = TRUE, showWarnings = FALSE)

# 4a. Geography: raw metric vs. migration distance.
p_geo <- ggplot(data.frame(x_disp = df$x_km / DISPLAY_UNIT_KM, y = y), aes(x_disp, y)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
  theme_bw() +
  labs(title = "Spanish grammar similarity vs. migration distance (uncontrolled)",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), R^2 = %.3f, RMSE = %.4f, n = %d",
                          unname(geo_cf["Estimate"]) * scl, DISPLAY_UNIT_KM,
                          unname(geo_cf["Pr(>|t|)"]), geo_r2, geo_perf["rmse"], N),
       x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish grammar similarity")
ggsave(here("figures", "regression", "grammar_regression_geography.png"),
       p_geo, width = 7, height = 4.5, units = "in", dpi = 300)

# 4b. Actual vs. predicted for both models, shared square axes, dashed y = x.
# B has no predictor axis, so the models are compared by how tightly each
# prediction tracks the metric (the RMSE printed in each panel).
avp_df <- bind_rows(
  tibble(panel = "geography (A)",               pred = fitted(m_geo), y = y),
  tibble(panel = "ancestry, geography-free (B)", pred = fitted(m_anc), y = y)
) |>
  mutate(panel = factor(panel, c("geography (A)", "ancestry, geography-free (B)")))
avp_lab <- tibble(
  panel = factor(c("geography (A)", "ancestry, geography-free (B)"),
                 c("geography (A)", "ancestry, geography-free (B)")),
  lab   = c(sprintf("RMSE %.4f (LOO %.4f)\nR^2 %.3f",
                    geo_perf["rmse"], geo_perf["rmse_loocv"], geo_r2),
            sprintf("RMSE %.4f (LOO %.4f)\nsemipartial R^2 %.3f",
                    anc_perf["rmse"], anc_perf["rmse_loocv"], anc_r2))
)
lim <- range(c(avp_df$pred, avp_df$y))
p_avp <- ggplot(avp_df, aes(pred, y)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2, alpha = 0.5) +
  geom_text(data = avp_lab, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.15, size = 3, lineheight = 0.95, inherit.aes = FALSE) +
  facet_wrap(~ panel) +
  coord_fixed(xlim = lim, ylim = lim) +
  theme_bw() +
  labs(title = "Predicting Spanish grammar similarity: geography vs. geography-free ancestry",
       subtitle = "actual vs. model-predicted; dashed line is y = x",
       x = "Predicted normalized similarity", y = "Actual normalized similarity")
ggsave(here("figures", "regression", "grammar_regression_actual_vs_predicted.png"),
       p_avp, width = 9, height = 4.8, units = "in", dpi = 300)


# 4c. vp_diagram(): PVR::VarPartplot()'s layout — total variance is the full-width
# line, geography is the bar above it (a + b) and ancestry the bar below (b + c),
# so the two overlap on the shared fraction b. Rebuilt here rather than called
# from PVR so it also serves domains with no PVR fit object.
vp_diagram <- function(vp, title, subtitle, anc_label) {
  f   <- setNames(vp$fraction, vp$component)[c("a", "b", "c", "d")]
  x   <- cumsum(c(0, f[["a"]], f[["b"]], f[["c"]]))
  mid <- x + f / 2

  # Narrow spans cannot sit under their own label, so labels are pushed apart to
  # a minimum gap and tied back to their span by a leader.
  gap <- 0.085
  o <- order(mid); p <- mid[o]
  for (i in seq_along(p)[-1]) p[i] <- max(p[i], p[i - 1] + gap)
  if (max(p) > 1) {
    p <- p - (max(p) - 1)
    for (i in rev(seq_along(p))[-1]) p[i] <- min(p[i], p[i + 1] - gap)
  }
  lx <- numeric(length(mid)); lx[o] <- p

  box <- data.frame(xmin = c(min(x[1], x[3]), min(x[2], x[4])),
                    xmax = c(max(x[1], x[3]), max(x[2], x[4])),
                    ymin = c(0.500, 0.375), ymax = c(0.625, 0.500))
  box <- box[box$xmax > box$xmin, ]   # a zero-width bar would draw as a stray tick

  ggplot() +
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.4) +
    geom_segment(data = data.frame(x = x), aes(x = x, xend = x, y = 0.5, yend = 0.68),
                 linetype = "dashed", linewidth = 0.3) +
    geom_segment(data = data.frame(mid = mid, lx = lx),
                 aes(x = mid, xend = lx, y = 0.68, yend = 0.755),
                 linewidth = 0.25, colour = "grey45") +
    annotate("segment", x = 0, xend = 1, y = 0.5, yend = 0.5, linewidth = 0.9) +
    geom_text(data = data.frame(lx = lx, lab = sprintf("%s\n%.3f", names(f), f)),
              aes(lx, 0.80, label = lab), size = 3.4, lineheight = 0.95, vjust = 0) +
    geom_text(data = data.frame(y = c(0.5625, 0.4375), lab = c("geography", anc_label)),
              aes(-0.012, y, label = lab), hjust = 1, size = 3.4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.34, 0.94), clip = "off") +
    theme_void(base_size = 11) +
    theme(plot.margin = margin(10, 14, 10, 96), plot.title.position = "plot") +
    labs(title = title, subtitle = subtitle)
}

ggsave(here("figures", "regression", "grammar_pvr_varpart.png"),
       vp_diagram(GRAMMAR_pvr_varpart, "Grammar variance partition",
                  sprintf("r^2 = %.3f", R2_full), "ancestry"),
       width = 7, height = 3.6, units = "in", dpi = 300)


# ── 5. Results table ──────────────────────────────────────────────────────
# One row per model. slope: A only (B is a k-vector block). statistic/p.value:
# A = t on the distance coefficient, B = the overall F-test of the ancestry
# model itself (y ~ scrubbed eigenvectors). r.squared: A = R^2, B = semipartial (ancestry-only).
# cor_with_distance: A = raw cor(metric, distance); B = cor(un-scrubbed ancestry
# prediction, distance) — sign shows whether ancestry runs with (-) or against
# (+) the geographic gradient.
GRAMMAR_geo_vs_ancestry <- tibble(
  model      = c("A_geography", "B_ancestry_given_geography"),
  predictor  = c(paste0(PREDICTOR, " (per ", DISPLAY_UNIT_KM, " km)"),
                 "phylo eigenvectors, migration distance scrubbed"),
  n = N, n_zeros = n_zeros, n_evec = c(0L, k),
  slope             = c(unname(geo_cf["Estimate"])   * scl, NA_real_),
  slope_se          = c(unname(geo_cf["Std. Error"]) * scl, NA_real_),
  statistic         = c(unname(geo_cf["t value"]), anc_F),
  p.value           = c(unname(geo_cf["Pr(>|t|)"]), anc_p),
  r.squared         = c(geo_r2, anc_r2),
  cor_with_distance = c(geo_dir, anc_dir),
  rmse              = c(unname(geo_perf["rmse"]),       unname(anc_perf["rmse"])),
  rmse_loocv        = c(unname(geo_perf["rmse_loocv"]), unname(anc_perf["rmse_loocv"]))
)
print(GRAMMAR_geo_vs_ancestry)
