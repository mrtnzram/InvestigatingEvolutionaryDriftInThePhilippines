# =============================================================================
# [4] Genetic Analysis — Population-structure eigenvector regression
#
# Regresses normalized Spanish admixture on terrain-penalized migration
# distance, controlling for population structure via PCA eigenvectors (PLINK
# --pca) in place of the phylogenetic eigenvectors the other three domains use
# — genetics has no tree (see [0]_GENETIC_subgroups.R).
#
# §8B additionally splits geography and population structure into two comparable
# regressions — y ~ migration distance (uncontrolled) and y ~ PC eigenvectors
# with geography scrubbed out of the eigenvectors — since the variance partition
# is too collinear to separate the two. Zero-admixture populations are dropped
# first (the spike at 0 dominates the slope), so this whole script runs on the
# non-zero set; §8C refits with the zeros restored so the slope/p difference is
# on record.
#
# Input:   data/genetic/network_distance/GENETIC_final.csv (from [3]),
#          data/Phil_2.24M_pca_results.eigenvec (PLINK --pca, individual-level,
#          headerless FID/IID/PC1..; same file [9]_GENETIC_PROCRUSTES.R reads)
# Outputs: data/genetic/PVR/GENETIC_pvr_results.csv, GENETIC_pvr_varpart.csv,
#          data/genetic/PVR/GENETIC_geo_vs_ancestry.csv (§8B regressions A and B),
#          data/genetic/PVR/GENETIC_geo_vs_ancestry_sensitivity.csv (§8C zeros kept vs. removed),
#          figures/genetic/regression/genetic_pvr_varpart.png,
#          figures/genetic/regression/genetic_pvr_partial_residuals.png,
#          figures/genetic/regression/genetic_regression_geography.png,
#          figures/genetic/regression/genetic_regression_ancestry.png,
#          figures/genetic/regression/genetic_regression_geography_zeros.png
# =============================================================================

library(tidyverse)
library(here)

source(here("R", "shared", "select_moran_eigenvectors.R"))

# H1 = penalized cost onto the nearest network node; H2 = cost to Manila.
PREDICTOR <- "geodist_H1_span"

# Fit on Mm; DISPLAY_UNIT_KM rescales reported and plotted slopes only.
DISPLAY_UNIT_KM <- 100

GENETIC_final <- read.csv(here("data", "genetic", "network_distance", "GENETIC_final.csv"))


# ── 1. Population-structure eigenvectors ────────────────────────────────────
# PLINK --pca is per individual; aggregate to per-population means (all 115
# populations covered, vs. 70/115 for a tree tip).
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

stopifnot(
  "eigenvec population set != GENETIC_final population set" =
    setequal(eigenvec_pop$population, GENETIC_final$population)
)


# ── 2. Analysis frame ────────────────────────────────────────────────────────
# Matched by name, not tip.label position — no phylo object here. qpAdm
# source populations carry no estimate and drop out.
df_full <- GENETIC_final |>
  dplyr::select(population, span_admx, span_w, x_km = all_of(PREDICTOR)) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  as.data.frame()

# Many populations carry zero Spanish admixture (negatives were floored to 0 in
# [0]_GENETIC_ADMX_MALDER.R); the spike at 0 dominates the slope. The primary
# analysis drops them — every model below, the §6 variance partition included,
# runs on the non-zero set. §8C refits with the zeros restored so the slope/p
# difference is on record. max(span_admx) is unchanged by the drop, so y is on
# the same scale in both frames.
n_zeros_removed <- sum(df_full$span_admx == 0)
df_full$y <- df_full$span_admx / max(df_full$span_admx)
df   <- df_full[df_full$span_admx > 0, , drop = FALSE]

y <- df$y
x <- df$x_km / 1000
E <- as.matrix(df[, paste0("PC", seq_len(n_pc)), drop = FALSE])
N <- nrow(df)

message(N, " populations with a non-zero admixture estimate and PCA coverage | ",
        n_zeros_removed, " zero-admixture populations dropped.")


# ── 3. Eigenvector selection (Moran's-I stopping rule) ──────────────────────
# select_moran_eigenvectors() (R/shared/) is genetics' analogue of PVR(method
# = "moran"). Shared with [4.1]_GENETIC_PVR_sPCA.R so both select the same PCs.
E_sel <- select_moran_eigenvectors(y, x, E)
k     <- ncol(E_sel)

# lm() chokes on a zero-column matrix term, so k = 0 (no PC clears the
# selection bar) drops E_sel from the formula rather than passing an empty one.
fit_full   <- function(y, E, x) if (ncol(E) == 0) lm(y ~ x) else lm(y ~ E + x)
fit_e_only <- function(y, E)    if (ncol(E) == 0) lm(y ~ 1) else lm(y ~ E)


# ── 4. Fit ────────────────────────────────────────────────────────────────
m_full <- fit_full(y, E_sel, x)

message("N = ", N, " populations | ", k, " PC", if (k == 1) "" else "s",
        " selected | residual df = ", df.residual(m_full))
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


# ── 5. Partial residuals (Frisch–Waugh–Lovell) ──────────────────────────────
# lm(res_y ~ res_x) reproduces the full model's x coefficient. Asserted so
# §9's scatter cannot drift from the reported effect.
res_y <- resid(fit_e_only(y, E_sel))
res_x <- resid(fit_e_only(x, E_sel))
m_av  <- lm(res_y ~ res_x)

stopifnot(
  "AV-plot slope does not match the full model's distance coefficient." =
    isTRUE(all.equal(unname(coef(m_av)["res_x"]), b_mm, tolerance = 1e-8))
)

partial_r <- cor(res_y, res_x)
raw_r     <- cor(y, df$x_km)
cat(sprintf("raw cor(admixture, distance) = %.4f | partial cor given population structure = %.4f\n",
            raw_r, partial_r))


# ── 6. Variance partition ───────────────────────────────────────────────────
# Manual Legendre-style nested-R^2 partition — no phylo object for PVR's
# @VarPart to use. a = geography, c = population structure.
r2 <- function(fit) summary(fit)$r.squared
R2_full <- r2(fit_full(y, E_sel, x))
R2_geo  <- r2(lm(y ~ x))
R2_pop  <- r2(fit_e_only(y, E_sel))

vp_a <- R2_full - R2_pop           # geography only
vp_c <- R2_full - R2_geo           # population structure only
vp_b <- R2_geo + R2_pop - R2_full  # shared (can go negative under suppression)
vp_d <- 1 - R2_full                # unexplained

GENETIC_pvr_varpart <- tibble(
  component = c("a", "b", "c", "d"),
  meaning   = c("geography only", "shared (geography & population structure)",
                "population structure only", "unexplained"),
  fraction  = c(vp_a, vp_b, vp_c, vp_d)
)
print(GENETIC_pvr_varpart)
stopifnot("Variance partition does not sum to 1." =
            isTRUE(all.equal(sum(GENETIC_pvr_varpart$fraction), 1, tolerance = 1e-8)))

write.csv(GENETIC_pvr_varpart,
          file = here("data", "genetic", "PVR", "GENETIC_pvr_varpart.csv"), row.names = FALSE)


# ── 7. Collinearity: geographic vs. population-structure distance ───────────
# Same diagnostic as the tree-based domains' cophenetic.phylo() comparison, using PC-space distance.
Dgeo <- as.matrix(dist(df$x_km))
Dpop <- as.matrix(dist(E))

geo_pop_cor <- cor(Dgeo[lower.tri(Dgeo)], Dpop[lower.tri(Dpop)])
message("cor(geographic distance, population-structure distance) = ", round(geo_pop_cor, 3))


# ── 8. Results table ────────────────────────────────────────────────────────
# Same column schema as the other three domains' *_pvr_results.csv.
GENETIC_pvr_results <- as_tibble(cf, rownames = "term") |>
  rename(estimate = Estimate, std.error = `Std. Error`,
         statistic = `t value`, p.value = `Pr(>|t|)`) |>
  filter(term %in% c("(Intercept)", "x")) |>
  mutate(term = if_else(term == "x", PREDICTOR, term),
         estimate  = if_else(term == PREDICTOR, estimate * scl, estimate),
         std.error = if_else(term == PREDICTOR, std.error * scl, std.error),
         unit = if_else(term == PREDICTOR, paste0("per ", DISPLAY_UNIT_KM, " km"), NA_character_),
         method = "moran (population-structure PCA)",
         n = N, n_evec = k, df_residual = df.residual(m_full),
         r.squared = summary(m_full)$r.squared,
         adj.r.squared = summary(m_full)$adj.r.squared,
         aic = AIC(m_full),
         raw_cor = raw_r, partial_cor = partial_r,
         geo_phylo_cor = geo_pop_cor, .before = 1)
print(GENETIC_pvr_results)
write.csv(GENETIC_pvr_results,
          file = here("data", "genetic", "PVR", "GENETIC_pvr_results.csv"), row.names = FALSE)


# ── 8B. Geography vs. ancestry — two comparable regressions ─────────────────
# The §6 partition leaves geography and population structure too collinear to
# separate. These two regressions split the overlap by design (zeros dropped in §2):
#   A  y ~ migration distance, uncontrolled — modern geography, credited with
#      everything it shares with population structure.
#   B  y ~ PC eigenvectors residualized on migration distance (y left untouched).
#      E_perp is orthogonal to x, so B's R^2 is the semipartial (structure-only)
#      share and geo_r2 + anc_r2 reproduces the full model R^2.

# A — modern geography
m_geo  <- lm(y ~ x)
geo_cf <- summary(m_geo)$coefficients["x", ]
geo_r2 <- summary(m_geo)$r.squared

# B — population structure, geography scrubbed out of the eigenvectors
if (k == 0) {
  m_anc  <- lm(y ~ 1)
  anc_r2 <- 0
  anc_F  <- NA_real_
  anc_p  <- NA_real_
} else {
  E_perp <- apply(E_sel, 2, \(col) resid(lm(col ~ x)))
  m_anc  <- lm(y ~ E_perp)
  anc_r2 <- summary(m_anc)$r.squared
  inc    <- anova(m_geo, m_full)   # incremental F: structure added on top of geography
  anc_F  <- inc[["F"]][2]
  anc_p  <- inc[["Pr(>F)"]][2]
  stopifnot(
    "Scrubbed-eigenvector coefficients drift from the full fit." =
      isTRUE(all.equal(unname(coef(m_anc)[-1]),
                       unname(coef(m_full)[2:(k + 1)]), tolerance = 1e-8)),
    "geo_r2 + anc_r2 != full-model R^2 — E_perp is not orthogonal to x." =
      isTRUE(all.equal(geo_r2 + anc_r2, summary(m_full)$r.squared, tolerance = 1e-6))
  )
}

message(sprintf(
  "R^2 geography (A) = %.3f | semipartial R^2 structure|geography (B) = %.3f | A+B = %.3f (full model R^2 = %.3f)",
  geo_r2, anc_r2, geo_r2 + anc_r2, summary(m_full)$r.squared))

GENETIC_geo_vs_ancestry <- tibble(
  model     = c("A_geography", "B_ancestry_given_geography"),
  predictor = c(paste0(PREDICTOR, " (per ", DISPLAY_UNIT_KM, " km)"),
                "PC eigenvectors ⊥ migration distance"),
  n = N, n_zeros_removed = n_zeros_removed,
  n_evec      = c(0L, k),
  estimate    = c(unname(geo_cf["Estimate"])   * scl, NA_real_),
  std.error   = c(unname(geo_cf["Std. Error"]) * scl, NA_real_),
  statistic   = c(unname(geo_cf["t value"]), anc_F),
  p.value     = c(unname(geo_cf["Pr(>|t|)"]), anc_p),
  r.squared   = c(geo_r2, anc_r2),
  df.residual = c(df.residual(m_geo), df.residual(m_full))
)
print(GENETIC_geo_vs_ancestry)
write.csv(GENETIC_geo_vs_ancestry,
          file = here("data", "genetic", "PVR", "GENETIC_geo_vs_ancestry.csv"),
          row.names = FALSE)


# ── 8C. Zeros-in sensitivity ───────────────────────────────────────────────
# Refit A and B with the zero-admixture populations restored so the migration-
# distance slope and its p-value can be read off with and without them. The PC
# selection is re-run on each set — its eigenvector count can differ.
geo_anc_fit <- function(frame) {
  yy  <- frame$y
  xx  <- frame$x_km / 1000
  Em  <- as.matrix(frame[, paste0("PC", seq_len(n_pc)), drop = FALSE])
  mg  <- lm(yy ~ xx)
  gcf <- summary(mg)$coefficients["xx", ]
  gr2 <- summary(mg)$r.squared
  Es  <- select_moran_eigenvectors(yy, xx, Em)
  kk  <- ncol(Es)
  mf  <- if (kk == 0) lm(yy ~ xx) else lm(yy ~ Es + xx)
  tibble(
    n = nrow(frame), n_evec = kk,
    geo_slope     = unname(gcf["Estimate"])   * scl,
    geo_std_error = unname(gcf["Std. Error"]) * scl,
    geo_t         = unname(gcf["t value"]),
    geo_p         = unname(gcf["Pr(>|t|)"]),
    geo_r2        = gr2,
    anc_semipartial_r2 = if (kk == 0) 0 else summary(mf)$r.squared - gr2,
    anc_F_p            = if (kk == 0) NA_real_ else anova(mg, mf)[["Pr(>F)"]][2]
  )
}

GENETIC_geo_vs_ancestry_sensitivity <- bind_rows(
  tibble(zeros = "removed") |> bind_cols(geo_anc_fit(df)),
  tibble(zeros = "kept")    |> bind_cols(geo_anc_fit(df_full))
)
stopifnot(
  "§8C 'removed' refit disagrees with the §8B primary fit." =
    isTRUE(all.equal(GENETIC_geo_vs_ancestry_sensitivity$geo_r2[1], geo_r2, tolerance = 1e-8)) &&
    isTRUE(all.equal(GENETIC_geo_vs_ancestry_sensitivity$anc_semipartial_r2[1], anc_r2, tolerance = 1e-6))
)
print(GENETIC_geo_vs_ancestry_sensitivity)
write.csv(GENETIC_geo_vs_ancestry_sensitivity,
          file = here("data", "genetic", "PVR", "GENETIC_geo_vs_ancestry_sensitivity.csv"),
          row.names = FALSE)


# ── 9. Figures ──────────────────────────────────────────────────────────────
dir.create(here("figures", "genetic", "regression"),
           recursive = TRUE, showWarnings = FALSE)

# ggplot bar chart, not PVR's VarPartplot() — no PVR fit object here.
p_vp <- ggplot(GENETIC_pvr_varpart, aes(x = component, y = fraction, fill = component)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", fraction)), vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = c(a = "#4C72B0", b = "grey70", c = "#DD8452", d = "grey90"),
                     guide = "none") +
  theme_bw() +
  labs(title = "Spanish admixture: variance partition",
       subtitle = "a = geography | b = shared | c = population structure | d = unexplained",
       x = NULL, y = expression(R^2 ~ fraction))
print(p_vp)
ggsave(here("figures", "genetic", "regression", "genetic_pvr_varpart.png"),
       p_vp, width = 7, height = 4.5, units = "in", dpi = 300)

# Added-variable scatter, axes residualized on E_sel; line is m_av (§5 asserts its slope).
# Point area = span_w, matching the [7] map.
av_df <- data.frame(res_x = res_x / scl, res_y = res_y, span_w = df$span_w)
band  <- predict(lm(res_y ~ res_x, data = av_df),
                 newdata = data.frame(res_x = seq(min(av_df$res_x), max(av_df$res_x),
                                                  length.out = 100)),
                 interval = "confidence") |>
  as.data.frame() |>
  mutate(res_x = seq(min(av_df$res_x), max(av_df$res_x), length.out = 100))

p_av <- ggplot() +
  geom_point(data = av_df, aes(x = res_x, y = res_y, size = span_w), alpha = 0.5) +
  scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
  geom_ribbon(data = band, aes(x = res_x, ymin = lwr, ymax = upr),
              fill = "steelblue", alpha = 0.25) +
  geom_line(data = band, aes(x = res_x, y = fit),
            color = "steelblue", linewidth = 1.2) +
  theme_bw() +
  labs(title = "Spanish admixture vs. migration distance, population structure removed",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), partial r = %.3f, %d PC%s",
                          b_mm * scl, DISPLAY_UNIT_KM, cf["x", "Pr(>|t|)"], partial_r, k,
                          if (k == 1) "" else "s"),
       x = sprintf("Migration distance residual (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish admixture residual")
print(p_av)

ggsave(here("figures", "genetic", "regression", "genetic_pvr_partial_residuals.png"),
       p_av, width = 7, height = 4.5, units = "in", dpi = 300)


# §8B regression A: raw admixture vs. migration distance, no structure control.
# Point area = span_w, matching the [7] map.
geo_df <- data.frame(x_disp = df$x_km / DISPLAY_UNIT_KM, y = y, span_w = df$span_w)
p_geo <- ggplot(geo_df, aes(x_disp, y)) +
  geom_point(aes(size = span_w), alpha = 0.5) +
  scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
  theme_bw() +
  labs(title = "Spanish admixture vs. migration distance (uncontrolled, zeros removed)",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), R^2 = %.3f, n = %d (%d zeros removed)",
                          unname(geo_cf["Estimate"]) * scl, DISPLAY_UNIT_KM,
                          unname(geo_cf["Pr(>|t|)"]), geo_r2, N, n_zeros_removed),
       x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish admixture")
print(p_geo)
ggsave(here("figures", "genetic", "regression", "genetic_regression_geography.png"),
       p_geo, width = 7, height = 4.5, units = "in", dpi = 300)

# §8B regression B: admixture vs. the geography-free population-structure
# prediction. The fit line is the identity by construction; the scatter is the
# semipartial R^2.
if (k == 0) {
  p_anc <- ggplot() + theme_void() +
    annotate("text", x = 0, y = 0,
             label = "No PC eigenvector cleared Moran selection (k = 0).")
} else {
  anc_df <- data.frame(fit = fitted(m_anc), y = y, span_w = df$span_w)
  p_anc <- ggplot(anc_df, aes(fit, y)) +
    geom_point(aes(size = span_w), alpha = 0.5) +
    scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
    theme_bw() +
    labs(title = "Spanish admixture vs. population structure, geography scrubbed out",
         subtitle = sprintf("semipartial R^2 = %.3f, F p = %.3g, %d PC%s",
                            anc_r2, anc_p, k, if (k == 1) "" else "s"),
         x = "Geography-free population-structure prediction of admixture",
         y = "Normalized Spanish admixture")
}
print(p_anc)
ggsave(here("figures", "genetic", "regression", "genetic_regression_ancestry.png"),
       p_anc, width = 7, height = 4.5, units = "in", dpi = 300)

# §8C: regression A side by side, zero-admixture populations kept vs. removed.
s <- GENETIC_geo_vs_ancestry_sensitivity
zc_df <- bind_rows(
  df_full |> mutate(zeros = "kept",    x_disp = x_km / DISPLAY_UNIT_KM),
  df      |> mutate(zeros = "removed", x_disp = x_km / DISPLAY_UNIT_KM)
)
p_zc <- ggplot(zc_df, aes(x_disp, y)) +
  geom_point(aes(size = span_w), alpha = 0.5) +
  scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.1) +
  facet_wrap(~ factor(zeros, c("kept", "removed"))) +
  theme_bw() +
  labs(title = "Spanish admixture vs. migration distance: zero-admixture populations kept vs. removed",
       subtitle = sprintf(
         "kept: slope %.4f/%dkm, p = %.3g, n = %d   |   removed: slope %.4f/%dkm, p = %.3g, n = %d",
         s$geo_slope[s$zeros == "kept"],    DISPLAY_UNIT_KM, s$geo_p[s$zeros == "kept"],    s$n[s$zeros == "kept"],
         s$geo_slope[s$zeros == "removed"], DISPLAY_UNIT_KM, s$geo_p[s$zeros == "removed"], s$n[s$zeros == "removed"]),
       x = sprintf("Migration distance (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish admixture")
print(p_zc)
ggsave(here("figures", "genetic", "regression", "genetic_regression_geography_zeros.png"),
       p_zc, width = 9, height = 4.5, units = "in", dpi = 300)
