# =============================================================================
# [4] Genetic Analysis — Bayesian inverse-variance weighted regression
#
# Fits linear, exponential-decay and cubic-spline models of normalized Spanish
# admixture vs. terrain-penalized migration distance with ulam(), each under an
# independent-normal likelihood whose per-population residual sd is scaled by the
# qpAdm weight from [3], then compares them by WAIC/PSIS-LOO.
#
# Input:   data/GENETIC_final.csv (from [3]_GENETIC_network_distance.R)
# Outputs: data/GENETIC_ols_results.csv   (OLS / WLS reference coefficient table)
#          data/GENETIC_waic_compare.csv  (WAIC model comparison + nRMSE)
#          figures/genetic/regression/genetic_{linear,exponential,spline}_model.png
# =============================================================================

library(tidyverse)
library(rethinking)
library(splines)
library(loo)
library(here)

# rethinking pulls in `posterior`, which masks sd() and filter(); nothing here
# uses rvars or time-series filtering, so the base/dplyr versions are the ones needed.
if (requireNamespace("conflicted", quietly = TRUE)) {
  conflicted::conflicts_prefer(stats::sd, .quiet = TRUE)
  conflicted::conflicts_prefer(dplyr::filter, .quiet = TRUE)
}

# geodist_H1_span matches [4]_PHONEME_PGLS.R; switch to "geodist_H2_span" to
# regress on distance to Manila instead of distance onto the network.
PREDICTOR <- "geodist_H1_span"

GENETIC_final <- read.csv(here("data", "GENETIC_final.csv"))


# ── 1. Analysis frame + WLS reference fit ───────────────────────────────────
# The three qpAdm source populations carry no estimate and drop out here.
# span_admx is min-max normalized anchored at 0, matching [1]'s cosine
# convention, so a slope reads as fraction of the observed range per unit distance.
df <- GENETIC_final |>
  dplyr::select(population, span_admx, span_w, dist = all_of(PREDICTOR)) |>
  filter(!is.na(span_admx), !is.na(span_w)) |>
  mutate(span_admx_norm = span_admx / max(span_admx)) |>
  as.data.frame()

message("N = ", nrow(df), " populations with a Spanish-admixture estimate.")

# Weights normalized to mean 1 so sigma stays on the y scale and the dexp prior
# stays calibrated.
w <- df$span_w / mean(df$span_w)

ols_fit <- lm(span_admx_norm ~ dist, data = df)
wls_fit <- lm(span_admx_norm ~ dist, data = df, weights = w)

tidy_lm <- function(fit, model_label) {
  as_tibble(summary(fit)$coefficients, rownames = "term") |>
    rename(estimate = Estimate, std.error = `Std. Error`,
           statistic = `t value`, p.value = `Pr(>|t|)`) |>
    mutate(model = model_label,
           r.squared = summary(fit)$r.squared,
           adj.r.squared = summary(fit)$adj.r.squared,
           aic = AIC(fit), n = nobs(fit), .before = 1)
}

GENETIC_ols_results <- bind_rows(
  tidy_lm(ols_fit, "OLS (unweighted)"),
  tidy_lm(wls_fit, "WLS (inverse variance)")
)
print(GENETIC_ols_results, n = Inf)
write.csv(GENETIC_ols_results,
          file = here("data", "GENETIC_ols_results.csv"), row.names = FALSE)


# ── 2. Model data ───────────────────────────────────────────────────────────
# Distance rescaled km -> Mm (thousands of km) so the decay/slope priors sit on
# an O(1) scale; all plotting converts back to km.
y      <- df$span_admx_norm
x_km   <- df$dist
x      <- x_km / 1000
N      <- length(y)
inv_sw <- 1 / sqrt(w)   # sigma_i = sigma * inv_sw[i]

dat <- list(N = N, y = y, x = x, inv_sw = inv_sw)


# ── 3. Model 1: linear ──────────────────────────────────────────────────────
# Priors are [4]_PHONEME_PGLS.R's generic weakly-informative ones on the
# normalized [0,1] scale. log_lik = FALSE: WAIC is built in §6 instead.
m_lin <- ulam(
  alist(
    y ~ normal(mu, s),
    vector[N]:mu <- a + b * x,
    vector[N]:s  <- sigma * inv_sw,
    a ~ dnorm(0.5, 0.3),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_lin, pars = c("a", "b", "sigma"))


# ── 4. Model 2: exponential decay ───────────────────────────────────────────
# b ~ dexp(1) is positive-only, so the model can express decay but never growth.
m_exp <- ulam(
  alist(
    y ~ normal(mu, s),
    vector[N]:mu <- a * exp(-b * x),
    vector[N]:s  <- sigma * inv_sw,
    a ~ dnorm(0.5, 0.3),
    b ~ dexp(1),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_exp, pars = c("a", "b", "sigma"))


# ── 5. Model 3: cubic spline ────────────────────────────────────────────────
# B is the bs() design matrix (5 columns) passed as data; w are the basis weights.
Bmat <- bs(x, df = 5)
dat_sp <- c(dat, list(B = matrix(as.numeric(Bmat), nrow = N)))

m_spline <- ulam(
  alist(
    y ~ normal(mu, s),
    vector[N]:mu <- a + B %*% w,
    vector[N]:s  <- sigma * inv_sw,
    a ~ dnorm(0.5, 0.3),
    vector[5]:w ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat_sp, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_spline, pars = c("a", "sigma"))


# ── 6. Model comparison (WAIC + PSIS-LOO via loo) ───────────────────────────
# Pointwise log-likelihood under the per-population sd sigma * inv_sw.
weighted_loglik <- function(post, mu_of) {
  S  <- length(post$sigma)
  ll <- matrix(NA_real_, nrow = S, ncol = N)
  for (s in seq_len(S)) {
    ll[s, ] <- dnorm(y, mu_of(post, s), post$sigma[s] * inv_sw, log = TRUE)
  }
  ll
}

post_lin <- extract.samples(m_lin)
post_exp <- extract.samples(m_exp)
post_sp  <- extract.samples(m_spline)

ll_lin <- weighted_loglik(post_lin, \(p, s) p$a[s] + p$b[s] * x)
ll_exp <- weighted_loglik(post_exp, \(p, s) p$a[s] * exp(-p$b[s] * x))
ll_sp  <- weighted_loglik(post_sp,  \(p, s) p$a[s] + as.numeric(Bmat %*% p$w[s, ]))

waic_list <- list(m_lin = loo::waic(ll_lin),
                  m_exp = loo::waic(ll_exp),
                  m_spline = loo::waic(ll_sp))
loo_list  <- list(m_lin = loo::loo(ll_lin),
                  m_exp = loo::loo(ll_exp),
                  m_spline = loo::loo(ll_sp))
cat("\n--- WAIC comparison ---\n");     print(loo::loo_compare(waic_list))
cat("\n--- PSIS-LOO comparison ---\n"); print(loo::loo_compare(loo_list))


# ── 7. Posterior predictions + 95% credible intervals ───────────────────────
# The ribbon is the 95% quantile interval of mu across posterior draws. Fitting
# stays on §2's Mm scale; DISPLAY_UNIT_KM only rescales what is plotted.
DISPLAY_UNIT_KM <- 100
xseq_km <- seq(min(x_km), max(x_km), length.out = 100)
xseq    <- xseq_km / 1000

mu_lin <- sapply(xseq, function(xx) post_lin$a + post_lin$b * xx)
mu_exp <- sapply(xseq, function(xx) post_exp$a * exp(-post_exp$b * xx))
# B_grid: bs() basis evaluated on the prediction grid, using the fitted knots.
B_grid <- predict(Bmat, newx = xseq)
# as.numeric(post_sp$a): extract.samples returns the scalar intercept as a
# 1-column matrix, which would clash with the (samples x grid) spline term.
mu_sp  <- t(apply(post_sp$w, 1, function(w_s) as.numeric(B_grid %*% w_s))) + as.numeric(post_sp$a)

summarise_mu <- function(mu_draws, xseq_disp) {
  data.frame(
    dist  = xseq_disp,
    mean  = apply(mu_draws, 2, mean),
    lower = apply(mu_draws, 2, quantile, probs = 0.025),
    upper = apply(mu_draws, 2, quantile, probs = 0.975)
  )
}
xseq_disp <- xseq_km / DISPLAY_UNIT_KM
ribbon_lin <- summarise_mu(mu_lin, xseq_disp)
ribbon_exp <- summarise_mu(mu_exp, xseq_disp)
ribbon_sp  <- summarise_mu(mu_sp,  xseq_disp)

# Point area carries span_w, matching the [7] map.
scatter_df <- data.frame(dist = x_km / DISPLAY_UNIT_KM, span_admx_norm = y,
                         span_w = df$span_w)

# nRMSE at the observed points, computed once here so the plot subtitles and
# §8's comparison table report the same numbers; weighted to match the likelihood.
mu_lin_obs <- sapply(x, function(xx) post_lin$a + post_lin$b * xx)
mu_exp_obs <- sapply(x, function(xx) post_exp$a * exp(-post_exp$b * xx))
mu_sp_obs  <- as.numeric(post_sp$a) + post_sp$w %*% t(matrix(as.numeric(Bmat), nrow = N))

y_range <- max(y) - min(y)
nrmse <- function(mu_draws) {
  sqrt(sum(w * (y - apply(mu_draws, 2, mean))^2) / sum(w)) / y_range
}
nrmse_lin <- nrmse(mu_lin_obs)
nrmse_exp <- nrmse(mu_exp_obs)
nrmse_sp  <- nrmse(mu_sp_obs)

# Posterior-mean equations for the plot subtitles, with slopes/decay rates
# converted from the Mm fitting scale to the DISPLAY_UNIT_KM plotting scale.
disp_per_mm <- DISPLAY_UNIT_KM / 1000
a_lin <- mean(post_lin$a); b_lin_disp <- mean(post_lin$b) * disp_per_mm
eq_lin <- sprintf("y = %.4f %s %.4f*x", a_lin, ifelse(b_lin_disp >= 0, "+", "-"), abs(b_lin_disp))

a_exp <- mean(post_exp$a); b_exp_disp <- mean(post_exp$b) * disp_per_mm
eq_exp <- sprintf("y = %.4f * exp(-%.4f*x)", a_exp, b_exp_disp)

a_sp <- mean(post_sp$a); w_sp <- colMeans(post_sp$w)
eq_sp <- paste0(
  "y = ", sprintf("%.4f", a_sp),
  paste0(sprintf(" %s %.4f*B%d(x)", ifelse(w_sp >= 0, "+", "-"), abs(w_sp), seq_along(w_sp)),
         collapse = "")
)

# subtitle wraps so the 6-term spline equation does not clip in the saved PNG.
plot_fit <- function(ribbon_df, title, eq, nrmse_val) {
  ggplot() +
    geom_point(data = scatter_df,
               aes(x = dist, y = span_admx_norm, size = span_w), alpha = 0.5) +
    geom_ribbon(data = ribbon_df, aes(x = dist, ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.25) +
    geom_line(data = ribbon_df, aes(x = dist, y = mean),
              color = "steelblue", linewidth = 1.2) +
    scale_size(trans = "log10", range = c(0.7, 3.5), guide = "none") +
    theme_bw() +
    labs(title = title,
         subtitle = str_wrap(sprintf("%s   |   weighted nRMSE = %.4f", eq, nrmse_val), width = 70),
         x = "Relative Migration Distance (100 km)", y = "Normalized Spanish Admixture")
}

p_lin <- plot_fit(ribbon_lin, "Linear (inverse-variance weighted)", eq_lin, nrmse_lin)
p_exp <- plot_fit(ribbon_exp, "Exponential decay (inverse-variance weighted)", eq_exp, nrmse_exp)
p_sp  <- plot_fit(ribbon_sp,  "Cubic spline (inverse-variance weighted)", eq_sp, nrmse_sp)
print(p_lin); print(p_exp); print(p_sp)

dir.create(here("figures", "genetic", "regression"),
           recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "genetic", "regression", "genetic_linear_model.png"),
       p_lin, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "genetic", "regression", "genetic_exponential_model.png"),
       p_exp, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "genetic", "regression", "genetic_spline_model.png"),
       p_sp,  width = 7, height = 4.5, units = "in", dpi = 300)


# ── 8. Comparison report (WAIC + normalized RMSE) ───────────────────────────
# nRMSE values (nrmse_lin/exp/sp) computed in §7 alongside the plot subtitles.
report <- tibble(
  model  = c("m_lin", "m_exp", "m_spline"),
  WAIC   = sapply(waic_list, \(w) w$estimates["waic", "Estimate"]),
  WAIC_SE = sapply(waic_list, \(w) w$estimates["waic", "SE"]),
  pWAIC  = sapply(waic_list, \(w) w$estimates["p_waic", "Estimate"]),
  LOOIC  = sapply(loo_list,  \(l) l$estimates["looic", "Estimate"]),
  nRMSE  = c(nrmse_lin, nrmse_exp, nrmse_sp)
) |>
  mutate(dWAIC = WAIC - min(WAIC), .after = WAIC) |>
  arrange(WAIC)
print(report)
write.csv(report, file = here("data", "GENETIC_waic_compare.csv"), row.names = FALSE)
