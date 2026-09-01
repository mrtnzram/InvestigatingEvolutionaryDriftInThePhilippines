# =============================================================================
# [4] Phoneme Analysis — Phylogenetic eigenvector regression (PVR)
#
# Regresses normalized Spanish-phoneme cosine similarity on terrain-penalized
# migration distance while controlling for shared ancestry with phylogenetic
# eigenvectors.
#
# §7B additionally splits geography and phylogeny into two comparable regressions
# — y ~ migration distance (uncontrolled) and y ~ phylogenetic eigenvectors with
# geography scrubbed out of the eigenvectors — since the variance partition is
# too collinear to separate the two.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/phoneme/network_distance/PHONEME_final.csv (from [3]_PHONEME_network_distance.R)
# Outputs: data/phoneme/PVR/PHONEME_pvr_results.csv      (coefficients + model summary)
#          data/phoneme/PVR/PHONEME_pvr_varpart.csv      (variance partition a/b/c/d)
#          data/phoneme/PVR/PHONEME_geo_vs_ancestry.csv  (§7B regressions A and B)
#          figures/phoneme/regression/phoneme_pvr_varpart.png
#          figures/phoneme/regression/phoneme_pvr_partial_residuals.png
#          figures/phoneme/regression/phoneme_regression_geography.png
#          figures/phoneme/regression/phoneme_regression_ancestry.png
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
PREDICTOR <- "geodist_H1_span"

# Fit on Mm; DISPLAY_UNIT_KM rescales reported and plotted slopes only.
DISPLAY_UNIT_KM <- 100

PHONEME_final <- read.csv(here("data", "phoneme", "network_distance", "PHONEME_final.csv"))

tip_map <- tree_df_matched |> dplyr::select(original, ph)


# ── 1. Analysis frame, ordered to tip.label ─────────────────────────────────
# PVR matches trait and envVar to eigenvector rows by POSITION — row i is
# tree_pruned$tip.label[i]. Reorder and assert; a mis-sort fits silently.
#
# tip_map duplicates multi-tip languages: one row per tip, values repeated.
df <- PHONEME_final |>
  dplyr::select(language, y = cossim_span_norm, x_km = all_of(PREDICTOR)) |>
  left_join(tip_map, by = c("language" = "ph")) |>
  filter(!is.na(original), !is.na(y), !is.na(x_km)) |>
  as.data.frame()

df <- df[match(tree_pruned$tip.label, df$original), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the tip set." =
    nrow(df) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(df$original) && identical(df$original, tree_pruned$tip.label)
)

y <- df$y
x <- df$x_km / 1000
N <- nrow(df)

# Cosine similarity is strictly positive — no zero spike to drop here (see the
# cognate/genetic [4] scripts, which do filter). Carried for the §7B CSV.
n_zeros_removed <- 0L

message(N, " tips covering ", n_distinct(df$language), " languages ",
        "(multi-dialect languages contribute one row per tip).")


# ── 2. Eigendecomposition + PVR fit ─────────────────────────────────────────
# scale = TRUE rescales eigenvalues only; eigenvectors are unchanged.
pvr_dec <- PVRdecomp(tree_pruned, scale = TRUE)

# "moran" adds eigenvectors until residual autocorrelation is n.s. (sig.t = 0.05);
# connectivity matrix left at the package default.
pvr_fit <- PVR(pvr_dec, phy = tree_pruned, trait = y, envVar = x, method = "moran")

E_sel <- as.matrix(pvr_fit@Selection$Vectors)
k     <- ncol(E_sel)


# ── 3. Explicit refit — source of the reported coefficients ─────────────
# PVR exposes only R2 and residuals, and @PVR$p is the FIRST EIGENVECTOR's
# p-value, not the predictor's. Refit the same design for the coefficient
# table; the checks below confirm it is the same fit.
m_full <- lm(y ~ E_sel + x)

stopifnot(
  "Refit residuals differ from PVR's — check the tip ordering or E_sel." =
    isTRUE(all.equal(unname(resid(m_full)), unname(pvr_fit@PVR$Residuals),
                     tolerance = 1e-8)),
  "Refit R2 differs from PVR's." =
    isTRUE(all.equal(summary(m_full)$r.squared, pvr_fit@PVR$R2, tolerance = 1e-8))
)

message("N = ", N, " tips | ", k, " eigenvectors selected (",
        paste(pvr_fit@Selection$Id$Vectors, collapse = ", "), ") | residual df = ",
        df.residual(m_full))
if (df.residual(m_full) < 10) {
  warning("Only ", df.residual(m_full), " residual df — the eigenvector set is ",
          "consuming most of the sample; treat the p-value with caution.")
}

cf   <- summary(m_full)$coefficients
b_mm <- unname(cf["x", "Estimate"])                 # per Mm, the fitting scale
scl  <- DISPLAY_UNIT_KM / 1000                      # Mm -> display units

print(round(cf[c("(Intercept)", "x"), ], 5))
cat(sprintf("\nslope = %.5f per %d km   (t = %.3f, p = %.4g)\n",
            b_mm * scl, DISPLAY_UNIT_KM,
            cf["x", "t value"], cf["x", "Pr(>|t|)"]))


# ── 4. Partial residuals (Frisch–Waugh–Lovell) ──────────────────────────────
# lm(res_y ~ res_x) reproduces the full model's x coefficient. Asserted so
# §8's scatter cannot drift from the reported effect.
res_y <- resid(lm(y ~ E_sel))
res_x <- resid(lm(x ~ E_sel))
m_av  <- lm(res_y ~ res_x)

stopifnot(
  "AV-plot slope does not match the PVR fit's distance coefficient." =
    isTRUE(all.equal(unname(coef(m_av)["res_x"]), b_mm, tolerance = 1e-8))
)

partial_r <- cor(res_y, res_x)
raw_r     <- cor(y, df$x_km)
cat(sprintf("raw cor(similarity, distance) = %.4f | partial cor given phylogeny = %.4f\n",
            raw_r, partial_r))


# ── 5. Variance partition ───────────────────────────────────────────────────
# PVR's lettering reverses Desdevises et al.: a = geography, c = phylogeny.
# b is a difference of R^2, not a variance, and goes negative under suppression;
# only the sum is constrained.
vp <- pvr_fit@VarPart
PHONEME_pvr_varpart <- tibble(
  component = c("a", "b", "c", "d"),
  meaning   = c("geography only", "shared (geography & phylogeny)",
                "phylogeny only", "unexplained"),
  fraction  = c(vp$a, vp$b, vp$c, vp$d)
)
print(PHONEME_pvr_varpart)
stopifnot("Variance partition does not sum to 1." =
            isTRUE(all.equal(sum(PHONEME_pvr_varpart$fraction), 1, tolerance = 1e-8)))

write.csv(PHONEME_pvr_varpart,
          file = here("data", "phoneme", "PVR", "PHONEME_pvr_varpart.csv"), row.names = FALSE)


# ── 6. Collinearity: geographic vs. phylogenetic distance ───────────────────
# Collinearity of the predictor as fitted (|x_i - x_j|) with patristic distance.
# Runs below the pairwise-matrix r ~ 0.68 in [5]_PHONEME_MMRR.R.
Dgeo   <- as.matrix(dist(df$x_km))
Dphylo <- cophenetic.phylo(tree_pruned)

geo_phylo_cor <- cor(Dgeo[lower.tri(Dgeo)], Dphylo[lower.tri(Dphylo)])
message("cor(geographic distance, phylogenetic distance) = ", round(geo_phylo_cor, 3))


# ── 7. Results table ────────────────────────────────────────────────────────
# Slope on the display scale; intercept is scale-free. n_evec / df_residual
# carried so a thin fit shows in the CSV.
PHONEME_pvr_results <- as_tibble(cf, rownames = "term") |>
  rename(estimate = Estimate, std.error = `Std. Error`,
         statistic = `t value`, p.value = `Pr(>|t|)`) |>
  filter(term %in% c("(Intercept)", "x")) |>
  mutate(term = if_else(term == "x", PREDICTOR, term),
         estimate  = if_else(term == PREDICTOR, estimate * scl, estimate),
         std.error = if_else(term == PREDICTOR, std.error * scl, std.error),
         unit = if_else(term == PREDICTOR, paste0("per ", DISPLAY_UNIT_KM, " km"), NA_character_),
         method = pvr_fit@Selection$Method,
         n = N, n_evec = k, df_residual = df.residual(m_full),
         r.squared = summary(m_full)$r.squared,
         adj.r.squared = summary(m_full)$adj.r.squared,
         aic = AIC(m_full),
         raw_cor = raw_r, partial_cor = partial_r,
         geo_phylo_cor = geo_phylo_cor, .before = 1)
print(PHONEME_pvr_results)
write.csv(PHONEME_pvr_results,
          file = here("data", "phoneme", "PVR", "PHONEME_pvr_results.csv"), row.names = FALSE)


# ── 7B. Geography vs. ancestry — two comparable regressions ─────────────────
# The §5 partition leaves geography and phylogeny too collinear to separate.
# These two regressions split the overlap by design:
#   A  y ~ migration distance, uncontrolled — modern geography, credited with
#      everything it shares with phylogeny.
#   B  y ~ phylogenetic eigenvectors residualized on migration distance (y left
#      untouched). E_perp is orthogonal to x, so B's R^2 is the semipartial
#      (phylogeny-only) share and geo_r2 + anc_r2 reproduces the full model R^2.

# A — modern geography
m_geo  <- lm(y ~ x)
geo_cf <- summary(m_geo)$coefficients["x", ]
geo_r2 <- summary(m_geo)$r.squared

# B — ancestry, geography scrubbed out of the eigenvectors
if (k == 0) {
  m_anc  <- lm(y ~ 1)
  anc_r2 <- 0
  anc_F  <- NA_real_
  anc_p  <- NA_real_
} else {
  E_perp <- apply(E_sel, 2, \(col) resid(lm(col ~ x)))
  m_anc  <- lm(y ~ E_perp)
  anc_r2 <- summary(m_anc)$r.squared
  inc    <- anova(m_geo, m_full)   # incremental F: phylogeny added on top of geography
  anc_F  <- inc[["F"]][2]
  anc_p  <- inc[["Pr(>F)"]][2]
  stopifnot(
    "Scrubbed-eigenvector coefficients drift from the full PVR fit." =
      isTRUE(all.equal(unname(coef(m_anc)[-1]),
                       unname(coef(m_full)[2:(k + 1)]), tolerance = 1e-8)),
    "geo_r2 + anc_r2 != full-model R^2 — E_perp is not orthogonal to x." =
      isTRUE(all.equal(geo_r2 + anc_r2, summary(m_full)$r.squared, tolerance = 1e-6))
  )
}

message(sprintf(
  "R^2 geography (A) = %.3f | semipartial R^2 phylogeny|geography (B) = %.3f | A+B = %.3f (full model R^2 = %.3f)",
  geo_r2, anc_r2, geo_r2 + anc_r2, summary(m_full)$r.squared))

PHONEME_geo_vs_ancestry <- tibble(
  model     = c("A_geography", "B_ancestry_given_geography"),
  predictor = c(paste0(PREDICTOR, " (per ", DISPLAY_UNIT_KM, " km)"),
                "phylo eigenvectors ⊥ migration distance"),
  n = N, n_zeros_removed = n_zeros_removed,
  n_evec      = c(0L, k),
  estimate    = c(unname(geo_cf["Estimate"])   * scl, NA_real_),
  std.error   = c(unname(geo_cf["Std. Error"]) * scl, NA_real_),
  statistic   = c(unname(geo_cf["t value"]), anc_F),
  p.value     = c(unname(geo_cf["Pr(>|t|)"]), anc_p),
  r.squared   = c(geo_r2, anc_r2),
  df.residual = c(df.residual(m_geo), df.residual(m_full))
)
print(PHONEME_geo_vs_ancestry)
write.csv(PHONEME_geo_vs_ancestry,
          file = here("data", "phoneme", "PVR", "PHONEME_geo_vs_ancestry.csv"),
          row.names = FALSE)


# ── 8. Figures ──────────────────────────────────────────────────────────────
dir.create(here("figures", "phoneme", "regression"),
           recursive = TRUE, showWarnings = FALSE)

# VarPartplot is base graphics and draws an unlabelled bar: png() device, and
# the subtitle carries the key.
png(here("figures", "phoneme", "regression", "phoneme_pvr_varpart.png"),
    width = 7, height = 4.5, units = "in", res = 300)
VarPartplot(pvr_fit)
title(main = "Spanish phoneme similarity: variance partition",
      sub = "a = geography | b = shared | c = phylogeny | d = unexplained")
dev.off()

# Added-variable scatter: both axes residualized on E_sel. Line is m_av
# (§4 asserts it carries the full model's slope); ribbon is its 95% CI.
av_df <- data.frame(res_x = res_x / scl, res_y = res_y)
band  <- predict(lm(res_y ~ res_x, data = av_df),
                 newdata = data.frame(res_x = seq(min(av_df$res_x), max(av_df$res_x),
                                                  length.out = 100)),
                 interval = "confidence") |>
  as.data.frame() |>
  mutate(res_x = seq(min(av_df$res_x), max(av_df$res_x), length.out = 100))

p_av <- ggplot() +
  geom_point(data = av_df, aes(x = res_x, y = res_y), size = 2, alpha = 0.5) +
  geom_ribbon(data = band, aes(x = res_x, ymin = lwr, ymax = upr),
              fill = "steelblue", alpha = 0.25) +
  geom_line(data = band, aes(x = res_x, y = fit),
            color = "steelblue", linewidth = 1.2) +
  theme_bw() +
  labs(title = "Spanish phoneme similarity vs. migration distance, phylogeny removed",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), partial r = %.3f, %d eigenvector%s",
                          b_mm * scl, DISPLAY_UNIT_KM, cf["x", "Pr(>|t|)"], partial_r, k,
                          if (k == 1) "" else "s"),
       x = sprintf("Migration distance residual (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish similarity residual")
print(p_av)

ggsave(here("figures", "phoneme", "regression", "phoneme_pvr_partial_residuals.png"),
       p_av, width = 7, height = 4.5, units = "in", dpi = 300)


# §7B regression A: raw metric vs. migration distance, no phylogeny control.
geo_df <- data.frame(x_disp = df$x_km / DISPLAY_UNIT_KM, y = y)
p_geo <- ggplot(geo_df, aes(x_disp, y)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
  theme_bw() +
  labs(title = "Spanish phoneme similarity vs. migration distance (uncontrolled)",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), R^2 = %.3f, n = %d (%d zeros removed)",
                          unname(geo_cf["Estimate"]) * scl, DISPLAY_UNIT_KM,
                          unname(geo_cf["Pr(>|t|)"]), geo_r2, N, n_zeros_removed),
       x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish phoneme similarity")
print(p_geo)
ggsave(here("figures", "phoneme", "regression", "phoneme_regression_geography.png"),
       p_geo, width = 7, height = 4.5, units = "in", dpi = 300)

# §7B regression B: metric vs. the geography-free phylogenetic prediction. The
# fit line is the identity by construction; the scatter is the semipartial R^2.
if (k == 0) {
  p_anc <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0,
             label = "No phylogenetic eigenvector cleared Moran selection (k = 0).")
} else {
  anc_df <- data.frame(fit = fitted(m_anc), y = y)
  p_anc <- ggplot(anc_df, aes(fit, y)) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
    theme_bw() +
    labs(title = "Spanish phoneme similarity vs. ancestry, geography scrubbed out",
         subtitle = sprintf("semipartial R^2 = %.3f, F p = %.3g, %d eigenvector%s",
                            anc_r2, anc_p, k, if (k == 1) "" else "s"),
         x = "Geography-free phylogenetic prediction of similarity",
         y = "Normalized Spanish phoneme similarity")
}
print(p_anc)
ggsave(here("figures", "phoneme", "regression", "phoneme_regression_ancestry.png"),
       p_anc, width = 7, height = 4.5, units = "in", dpi = 300)
