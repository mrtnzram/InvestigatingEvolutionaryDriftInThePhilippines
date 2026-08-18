# =============================================================================
# [4] Genetic Analysis — Population-structure eigenvector regression
#
# Regresses normalized Spanish admixture on terrain-penalized migration
# distance, controlling for population structure via PCA eigenvectors (PLINK
# --pca) in place of the phylogenetic eigenvectors the other three domains use
# — genetics has no tree (see [0]_GENETIC_subgroups.R).
#
# Input:   data/GENETIC_final.csv (from [3]), data/pca_results_phil_only.eigenvec
#          (PLINK --pca, individual-level, headerless FID/IID/PC1..PC10)
# Outputs: data/GENETIC_pvr_results.csv, data/GENETIC_pvr_varpart.csv,
#          figures/genetic/regression/genetic_pvr_varpart.png,
#          figures/genetic/regression/genetic_pvr_partial_residuals.png
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
eigenvec_raw <- read.table(here("data", "genetic", "PVR", "pca_results_phil_only.eigenvec"),
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
df <- GENETIC_final |>
  dplyr::select(population, span_admx, span_w, x_km = all_of(PREDICTOR)) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  mutate(y = span_admx / max(span_admx)) |>
  as.data.frame()

y <- df$y
x <- df$x_km / 1000
E <- as.matrix(df[, paste0("PC", seq_len(n_pc)), drop = FALSE])
N <- nrow(df)

message(N, " populations with both an admixture estimate and PCA coverage.")


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
