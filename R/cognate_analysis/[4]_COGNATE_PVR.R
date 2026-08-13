# =============================================================================
# [4] Cognate Analysis — Phylogenetic eigenvector regression (PVR)
#
# Regresses the normalized Spanish loanword count on terrain-penalized migration
# distance while controlling for shared ancestry with phylogenetic eigenvectors.
#
# Replaces [4]_COGNATE_PGLS.R: phylogeny moves from the error structure into the
# mean structure, so the geography/ancestry confound is partitioned rather than
# absorbed — a = geography only, b = shared, c = phylogeny only, d = unexplained.
#
# Run order: requires `tree_pruned` and `tip_map` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/cognate/COGNATE_final.csv (from [3]_COGNATE_network_distance.R)
# Outputs: data/cognate/COGNATE_pvr_results.csv  (coefficients + model summary)
#          data/cognate/COGNATE_pvr_varpart.csv  (variance partition a/b/c/d)
#          figures/cognate/regression/cognate_pvr_varpart.png
#          figures/cognate/regression/cognate_pvr_partial_residuals.png
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tip_map` is not defined." =
    exists("tip_map")
)

# H1 = penalized cost onto the nearest network node; H2 = cost to Manila.
PREDICTOR <- "geodist_H1_span"

# Fit on Mm; DISPLAY_UNIT_KM rescales reported and plotted slopes only.
DISPLAY_UNIT_KM <- 100

COGNATE_final <- read.csv(here("data", "cognate", "COGNATE_final.csv"))


# ── 1. Analysis frame, ordered to tip.label ─────────────────────────────────
# PVR matches trait and envVar to eigenvector rows by POSITION — row i is
# tree_pruned$tip.label[i]. Reorder and assert; a mis-sort fits silently.
#
# [0] relabelled tips to glottocodes, one per language — no duplication here.
# semi_join is redundant against [3]'s output but guards a regenerated file.
df <- COGNATE_final |>
  dplyr::select(glottocode, language, y = loans_norm, x_km = all_of(PREDICTOR)) |>
  semi_join(tip_map, by = "glottocode") |>
  filter(!is.na(y), !is.na(x_km)) |>
  as.data.frame()

df <- df[match(tree_pruned$tip.label, df$glottocode), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the language set." =
    nrow(df) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(df$glottocode) && identical(df$glottocode, tree_pruned$tip.label)
)

y <- df$y
x <- df$x_km / 1000
N <- nrow(df)


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

message("N = ", N, " languages | ", k, " eigenvectors selected (",
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
cat(sprintf("raw cor(loans, distance) = %.4f | partial cor given phylogeny = %.4f\n",
            raw_r, partial_r))


# ── 5. Variance partition ───────────────────────────────────────────────────
# PVR's lettering reverses Desdevises et al.: a = geography, c = phylogeny.
# b is a difference of R^2, not a variance, and goes negative under suppression;
# only the sum is constrained.
vp <- pvr_fit@VarPart
COGNATE_pvr_varpart <- tibble(
  component = c("a", "b", "c", "d"),
  meaning   = c("geography only", "shared (geography & phylogeny)",
                "phylogeny only", "unexplained"),
  fraction  = c(vp$a, vp$b, vp$c, vp$d)
)
print(COGNATE_pvr_varpart)
stopifnot("Variance partition does not sum to 1." =
            isTRUE(all.equal(sum(COGNATE_pvr_varpart$fraction), 1, tolerance = 1e-8)))

write.csv(COGNATE_pvr_varpart,
          file = here("data", "cognate", "COGNATE_pvr_varpart.csv"), row.names = FALSE)


# ── 6. Collinearity: geographic vs. phylogenetic distance ───────────────────
# Collinearity of the predictor as fitted (|x_i - x_j|) with patristic distance.
# Runs below the pairwise-matrix r ~ 0.7 the phoneme/grammar MMRR scripts report.
Dgeo   <- as.matrix(dist(df$x_km))
Dphylo <- cophenetic.phylo(tree_pruned)

geo_phylo_cor <- cor(Dgeo[lower.tri(Dgeo)], Dphylo[lower.tri(Dphylo)])
message("cor(geographic distance, phylogenetic distance) = ", round(geo_phylo_cor, 3))


# ── 7. Results table ────────────────────────────────────────────────────────
# Slope on the display scale; intercept is scale-free. n_evec / df_residual
# carried so a thin fit shows in the CSV.
COGNATE_pvr_results <- as_tibble(cf, rownames = "term") |>
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
print(COGNATE_pvr_results)
write.csv(COGNATE_pvr_results,
          file = here("data", "cognate", "COGNATE_pvr_results.csv"), row.names = FALSE)


# ── 8. Figures ──────────────────────────────────────────────────────────────
dir.create(here("figures", "cognate", "regression"),
           recursive = TRUE, showWarnings = FALSE)

# VarPartplot is base graphics and draws an unlabelled bar: png() device, and
# the subtitle carries the key.
png(here("figures", "cognate", "regression", "cognate_pvr_varpart.png"),
    width = 7, height = 4.5, units = "in", res = 300)
VarPartplot(pvr_fit)
title(main = "Spanish loanwords: variance partition",
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
  labs(title = "Spanish loanwords vs. migration distance, phylogeny removed",
       subtitle = sprintf("slope = %.4f per %d km (p = %.3g), partial r = %.3f, %d eigenvector%s",
                          b_mm * scl, DISPLAY_UNIT_KM, cf["x", "Pr(>|t|)"], partial_r, k,
                          if (k == 1) "" else "s"),
       x = sprintf("Migration distance residual (%d km)", DISPLAY_UNIT_KM),
       y = "Normalized Spanish loanword count residual")
print(p_av)

ggsave(here("figures", "cognate", "regression", "cognate_pvr_partial_residuals.png"),
       p_av, width = 7, height = 4.5, units = "in", dpi = 300)
