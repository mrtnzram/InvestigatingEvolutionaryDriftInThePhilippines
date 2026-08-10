# =============================================================================
# [4] Cognate Analysis — Bayesian distance-decay regression
#
# Fits linear, exponential-decay and cubic-spline models of the normalized
# Spanish loanword count vs. terrain-penalized migration distance with ulam(),
# each under an independent-normal likelihood with a single residual sd, then
# compares them by WAIC/PSIS-LOO.
#
# Input:   data/cognate/COGNATE_final.csv (from [3]_COGNATE_network_distance.R)
# Outputs: data/cognate/COGNATE_ols_results.csv   (OLS reference coefficient table)
#          data/cognate/COGNATE_waic_compare.csv  (WAIC model comparison + nRMSE)
#          figures/cognate/regression/cognate_{linear,exponential,spline}_model.png
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

# geodist_H1_span matches [4]_GENETIC_OLS.R; switch to "geodist_H2_span" to
# regress on distance to Manila instead of distance onto the network.
PREDICTOR <- "geodist_H1_span"

COGNATE_final <- read.csv(here("data", "cognate", "COGNATE_final.csv"))


# ── 1. Analysis frame + OLS reference fit ───────────────────────────────────
# loans_norm is the min-max normalized loan count from [0]; there is no standard
# error on this side, so the fit is unweighted throughout.
df <- COGNATE_final |>
  dplyr::select(language, loans_norm, dist = all_of(PREDICTOR)) |>
  filter(!is.na(loans_norm), !is.na(dist)) |>
  as.data.frame()

message("N = ", nrow(df), " languages with a loanword count and a distance.")

ols_fit <- lm(loans_norm ~ dist, data = df)

tidy_lm <- function(fit, model_label) {
  as_tibble(summary(fit)$coefficients, rownames = "term") |>
    rename(estimate = Estimate, std.error = `Std. Error`,
           statistic = `t value`, p.value = `Pr(>|t|)`) |>
    mutate(model = model_label,
           r.squared = summary(fit)$r.squared,
           adj.r.squared = summary(fit)$adj.r.squared,
           aic = AIC(fit), n = nobs(fit), .before = 1)
}

COGNATE_ols_results <- tidy_lm(ols_fit, "OLS (unweighted)")
print(COGNATE_ols_results, n = Inf)
write.csv(COGNATE_ols_results,
          file = here("data", "cognate", "COGNATE_ols_results.csv"), row.names = FALSE)


# ── 2. Model data ───────────────────────────────────────────────────────────
# Distance rescaled km -> Mm (thousands of km) so the decay/slope priors sit on
# an O(1) scale; all plotting converts back to km.
y    <- df$loans_norm
x_km <- df$dist
x    <- x_km / 1000
N    <- length(y)

dat <- list(N = N, y = y, x = x)

# Fitting stays on the Mm scale above; DISPLAY_UNIT_KM only rescales what is
# plotted and printed, so coefficients read in "per 100 km".
DISPLAY_UNIT_KM <- 100
DISP_PER_MM     <- DISPLAY_UNIT_KM / 1000  # e.g. 0.1 when DISPLAY_UNIT_KM = 100

# Posterior-mean equation per model form, printed under that model's precis() so
# the console shows the fitted curve next to the parameters it came from.
eq_linear <- function(post) {
  b_disp <- mean(post$b) * DISP_PER_MM
  sprintf("y = %.4f %s %.4f*x", mean(post$a),
          ifelse(b_disp >= 0, "+", "-"), abs(b_disp))
}
eq_exponential <- function(post) {
  sprintf("y = %.4f * exp(-%.4f*x)", mean(post$a), mean(post$b) * DISP_PER_MM)
}
eq_spline <- function(post) {
  w_mean <- colMeans(post$w)
  paste0("y = ", sprintf("%.4f", mean(post$a)),
         paste0(sprintf(" %s %.4f*B%d(x)", ifelse(w_mean >= 0, "+", "-"),
                        abs(w_mean), seq_along(w_mean)),
                collapse = ""))
}


# ── 3. Model 1: linear ──────────────────────────────────────────────────────
# Priors are [4]_GENETIC_OLS.R's generic weakly-informative ones on the
# normalized [0,1] scale. log_lik = FALSE: WAIC is built in §6 instead.
m_lin <- ulam(
  alist(
    y ~ normal(mu, sigma),
    vector[N]:mu <- a + b * x,
    a ~ dnorm(0.5, 0.3),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
post_lin <- extract.samples(m_lin)
precis(m_lin, pars = c("a", "b", "sigma"))
cat("equation:", eq_linear(post_lin), "\n")


# ── 4. Model 2: exponential decay ───────────────────────────────────────────
# b ~ dexp(1) is positive-only, so the model can express decay but never growth.
m_exp <- ulam(
  alist(
    y ~ normal(mu, sigma),
    vector[N]:mu <- a * exp(-b * x),
    a ~ dnorm(0.5, 0.3),
    b ~ dexp(1),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
post_exp <- extract.samples(m_exp)
precis(m_exp, pars = c("a", "b", "sigma"))
cat("equation:", eq_exponential(post_exp), "\n")


# ── 5. Model 3: cubic spline ────────────────────────────────────────────────
# B is the bs() design matrix (5 columns) passed as data; w are the basis weights.
Bmat <- bs(x, df = 5)
dat_sp <- c(dat, list(B = matrix(as.numeric(Bmat), nrow = N)))

m_spline <- ulam(
  alist(
    y ~ normal(mu, sigma),
    vector[N]:mu <- a + B %*% w,
    a ~ dnorm(0.5, 0.3),
    vector[5]:w ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat_sp, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
post_sp <- extract.samples(m_spline)
precis(m_spline, pars = c("a", "sigma"))
cat("equation:", eq_spline(post_sp), "\n")


# ── 6. Model comparison (WAIC + PSIS-LOO via loo) ───────────────────────────
pointwise_loglik <- function(post, mu_of) {
  S  <- length(post$sigma)
  ll <- matrix(NA_real_, nrow = S, ncol = N)
  for (s in seq_len(S)) {
    ll[s, ] <- dnorm(y, mu_of(post, s), post$sigma[s], log = TRUE)
  }
  ll
}

# post_lin / post_exp / post_sp were extracted alongside each fit above.
ll_lin <- pointwise_loglik(post_lin, \(p, s) p$a[s] + p$b[s] * x)
ll_exp <- pointwise_loglik(post_exp, \(p, s) p$a[s] * exp(-p$b[s] * x))
ll_sp  <- pointwise_loglik(post_sp,  \(p, s) p$a[s] + as.numeric(Bmat %*% p$w[s, ]))

waic_list <- list(m_lin = loo::waic(ll_lin),
                  m_exp = loo::waic(ll_exp),
                  m_spline = loo::waic(ll_sp))
loo_list  <- list(m_lin = loo::loo(ll_lin),
                  m_exp = loo::loo(ll_exp),
                  m_spline = loo::loo(ll_sp))
cat("\n--- WAIC comparison ---\n");     print(loo::loo_compare(waic_list))
cat("\n--- PSIS-LOO comparison ---\n"); print(loo::loo_compare(loo_list))


# ── 7. Posterior predictions + 95% credible intervals ───────────────────────
# The ribbon is the 95% quantile interval of mu across posterior draws;
# DISPLAY_UNIT_KM (set in §2) rescales the x axis to hundreds of km.
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

# Fixed point size, matching the [7] map: with no standard error there is
# nothing for point area to encode.
scatter_df <- data.frame(dist = x_km / DISPLAY_UNIT_KM, loans_norm = y)

# nRMSE at the observed points, for §8\'s comparison table.
mu_lin_obs <- sapply(x, function(xx) post_lin$a + post_lin$b * xx)
mu_exp_obs <- sapply(x, function(xx) post_exp$a * exp(-post_exp$b * xx))
mu_sp_obs  <- as.numeric(post_sp$a) + post_sp$w %*% t(matrix(as.numeric(Bmat), nrow = N))

y_range <- max(y) - min(y)
nrmse <- function(mu_draws) {
  sqrt(mean((y - apply(mu_draws, 2, mean))^2)) / y_range
}
nrmse_lin <- nrmse(mu_lin_obs)
nrmse_exp <- nrmse(mu_exp_obs)
nrmse_sp  <- nrmse(mu_sp_obs)

plot_fit <- function(ribbon_df, title) {
  ggplot() +
    geom_point(data = scatter_df,
               aes(x = dist, y = loans_norm), size = 2, alpha = 0.5) +
    geom_ribbon(data = ribbon_df, aes(x = dist, ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.25) +
    geom_line(data = ribbon_df, aes(x = dist, y = mean),
              color = "steelblue", linewidth = 1.2) +
    theme_bw() +
    labs(title = title,
         x = "Relative Migration Distance (100 km)",
         y = "Normalized Spanish Loanword Count")
}

p_lin <- plot_fit(ribbon_lin, "Linear")
p_exp <- plot_fit(ribbon_exp, "Exponential decay")
p_sp  <- plot_fit(ribbon_sp,  "Cubic spline")
print(p_lin); print(p_exp); print(p_sp)

dir.create(here("figures", "cognate", "regression"),
           recursive = TRUE, showWarnings = FALSE)

ggsave(here("figures", "cognate", "regression", "cognate_linear_model.png"),
       p_lin, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "cognate", "regression", "cognate_exponential_model.png"),
       p_exp, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "cognate", "regression", "cognate_spline_model.png"),
       p_sp,  width = 7, height = 4.5, units = "in", dpi = 300)


# ── 8. Comparison report (WAIC + normalized RMSE) ───────────────────────────
# nRMSE values (nrmse_lin/exp/sp) computed in §7.
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
write.csv(report, file = here("data", "cognate", "COGNATE_waic_compare.csv"),
          row.names = FALSE)
