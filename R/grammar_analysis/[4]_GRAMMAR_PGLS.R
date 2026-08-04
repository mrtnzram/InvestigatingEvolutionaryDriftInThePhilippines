# =============================================================================
# [4] Grammar Analysis — Bayesian phylogenetic regression
#
# Fits linear, exponential-decay and cubic-spline models of normalized
# Spanish-grammar cosine similarity (cossim_span_norm, min-max normalized 0 to
# max(cossim_span) in [1] — see there) vs. terrain-penalized migration distance
# with ulam(), each under a multi_normal likelihood carrying the same
# Pagel's-lambda-scaled phylogenetic covariance, then compares them by
# WAIC/PSIS-LOO. Lambda is fixed at the PGLS ML estimate and shared across
# models, so WAIC differences reflect only the mean structure. Fitting on the
# normalized scale means a slope reads as "fraction of the observed Spanish-
# similarity range per unit distance" rather than a raw cosine-similarity delta.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/GRAMMAR_final.csv (from [3]_GRAMMAR_network_distance.R)
# Outputs: data/GRAMMAR_pgls_results.csv   (PGLS lambda / slope table)
#          data/GRAMMAR_waic_compare.csv   (WAIC model comparison + nRMSE)
#          figures/grammar/regression/grammar_{linear,exponential,spline}_model.png
# =============================================================================

library(caper)
library(ape)
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

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

GRAMMAR_final <- read.csv(here("data", "GRAMMAR_final.csv"))

# dplyr:: qualified throughout: caper loads MASS, whose select() masks
# dplyr::select() once attached, regardless of library() order.
tip_map <- tree_df_matched |> dplyr::select(original, gram)


# ── 1. Comparative dataset + PGLS lambda estimate ───────────────────────────
# comparative.data(vcv = TRUE) returns the VCV already ordered to match its
# $data, and that ordering is reused below so y, x and R_lambda stay in sync.
# The tip_map join duplicates multi-tip languages, giving one point per tip.
df <- GRAMMAR_final |>
  dplyr::select(language, cossim = cossim_span_norm, geodist_H1_span) |>
  left_join(tip_map, by = c("language" = "gram")) |>
  filter(!is.na(original)) |>
  as.data.frame()

comp_data <- comparative.data(
  phy = tree_pruned, data = df, names.col = "original", vcv = TRUE
)

# PGLS is the lambda estimator. Fixed-lambda~0 fit is retained only for the
# reporting table (near-OLS reference); pgls_fit (ML) supplies lambda_hat.
pgls_ols <- pgls(cossim ~ geodist_H1_span, data = comp_data, lambda = 1e-6)
pgls_fit <- pgls(cossim ~ geodist_H1_span, data = comp_data, lambda = "ML")
lambda_hat <- unname(pgls_fit$param["lambda"])
message("PGLS ML lambda_hat = ", round(lambda_hat, 4))

tidy_pgls <- function(fit, model_label) {
  as_tibble(summary(fit)$coefficients, rownames = "term") |>
    rename(estimate = Estimate, std.error = `Std. Error`,
           statistic = `t value`, p.value = `Pr(>|t|)`) |>
    mutate(model = model_label,
           lambda = unname(fit$param["lambda"]),
           r.squared = summary(fit)$r.squared,
           adj.r.squared = summary(fit)$adj.r.squared,
           aic = as.numeric(fit$aic), n = fit$n, .before = 1)
}

GRAMMAR_pgls_results <- bind_rows(
  tidy_pgls(pgls_ols, "OLS (lambda fixed ~0)"),
  tidy_pgls(pgls_fit, "ML (lambda estimated)")
)
print(GRAMMAR_pgls_results, n = Inf)
write.csv(GRAMMAR_pgls_results,
          file = here("data", "GRAMMAR_pgls_results.csv"), row.names = FALSE)


# ── 2. Build the fixed phylogenetic correlation R_lambda ────────────────────
# R_lambda = lambda*R + (1-lambda)*I on the correlation-scaled VCV; sigma^2 in
# the models scales it up to the residual covariance. The matrix must be stripped
# of its caper "VCV.array" class and dimnames, which break ulam's dimension registration.
ord  <- rownames(comp_data$vcv)
cd   <- comp_data$data[ord, , drop = FALSE]
y    <- cd$cossim
x_km <- cd$geodist_H1_span

V_phylo <- unclass(comp_data$vcv)[ord, ord]  # raw phylogenetic VCV, kept for the collinearity check in §9
R <- cov2cor(V_phylo)
Rlambda <- lambda_hat * R + (1 - lambda_hat) * diag(nrow(R))
Rlambda <- matrix(as.numeric(Rlambda), nrow = nrow(Rlambda))  # plain matrix, no attrs
stopifnot(length(y) == nrow(Rlambda))

# Distance rescaled km -> Mm (thousands of km) so the decay/slope priors sit on
# an O(1) scale; all plotting converts back to km.
x <- x_km / 1000
N <- length(y)

dat <- list(N = N, y = y, x = x, Rlambda = Rlambda)


# ── 3. Model 1: linear ──────────────────────────────────────────────────────
# Priors are generic weakly-informative ones on the normalized [0,1] scale,
# wide enough to stay agnostic about where on the range the intercept and
# effect sizes fall. log_lik = FALSE because WAIC is built from the whitened
# likelihood in §6 instead.
m_lin <- ulam(
  alist(
    y ~ multi_normal(mu, K),
    vector[N]:mu  <- a + b * x,
    matrix[N,N]:K <- square(sigma) * Rlambda,
    a ~ dnorm(0.5, 0.3),
    b ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_lin, pars = c("a", "b", "sigma"))


# ── 4. Model 2: exponential decay ───────────────────────────────────────────
# a is on the normalized [0,1] scale, so it gets the same generic wide prior as
# the linear model's intercept. b ~ dexp(1) is a positive-only prior on the
# distance-decay rate (an x-scale parameter, unaffected by normalizing y), so
# the model can only express decay (mu falling with distance), never growth.
m_exp <- ulam(
  alist(
    y ~ multi_normal(mu, K),
    vector[N]:mu  <- a * exp(-b * x),
    matrix[N,N]:K <- square(sigma) * Rlambda,
    a ~ dnorm(0.5, 0.3),
    b ~ dexp(1),
    sigma ~ dexp(2)
  ),
  data = dat, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_exp, pars = c("a", "b", "sigma"))


# ── 5. Model 3: cubic spline ────────────────────────────────────────────────
# B is the bs() design matrix (5 columns) passed as data; w are the basis
# weights, given the same generic wide prior as the other models' effect terms.
Bmat <- bs(x, df = 5)
dat_sp <- c(dat, list(B = matrix(as.numeric(Bmat), nrow = N)))

m_spline <- ulam(
  alist(
    y ~ multi_normal(mu, K),
    vector[N]:mu  <- a + B %*% w,
    matrix[N,N]:K <- square(sigma) * Rlambda,
    a ~ dnorm(0.5, 0.3),
    vector[5]:w ~ dnorm(0, 0.3),
    sigma ~ dexp(2)
  ),
  data = dat_sp, chains = 4, cores = 4, cmdstan = TRUE, log_lik = FALSE
)
precis(m_spline, pars = c("a", "sigma"))


# ── 6. Model comparison (WAIC + PSIS-LOO via loo) ───────────────────────────
# Pointwise log-likelihood by Cholesky whitening: Lc is the lower factor of the
# fixed R_lambda, and forwardsolve() applies Lc^{-1} to the residuals.
Lc     <- t(chol(Rlambda))        # lower-triangular: Rlambda = Lc %*% t(Lc)
Lc_diag <- diag(Lc)

phylo_loglik <- function(post, mu_of) {
  S  <- length(post$sigma)
  ll <- matrix(NA_real_, nrow = S, ncol = N)
  for (s in seq_len(S)) {
    z <- forwardsolve(Lc, y - mu_of(post, s)) / post$sigma[s]
    ll[s, ] <- dnorm(z, 0, 1, log = TRUE) - log(post$sigma[s] * Lc_diag)
  }
  ll
}

post_lin <- extract.samples(m_lin)
post_exp <- extract.samples(m_exp)
post_sp  <- extract.samples(m_spline)

ll_lin <- phylo_loglik(post_lin, \(p, s) p$a[s] + p$b[s] * x)
ll_exp <- phylo_loglik(post_exp, \(p, s) p$a[s] * exp(-p$b[s] * x))
ll_sp  <- phylo_loglik(post_sp,  \(p, s) p$a[s] + as.numeric(Bmat %*% p$w[s, ]))

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
# stays on the Mm scale of §2 (the priors are tuned to it); DISPLAY_UNIT_KM only
# rescales what is plotted and printed, so coefficients read in "per 100 km".
DISPLAY_UNIT_KM <- 100
xseq_km <- seq(min(x_km), max(x_km), length.out = 100)
xseq    <- xseq_km / 1000

# post_lin / post_exp / post_sp already extracted in §6.
# mu draws: rows = posterior samples, cols = grid points.
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

scatter_df <- data.frame(geodist_H1_span = x_km / DISPLAY_UNIT_KM, cossim_span_norm = y)

# nRMSE at the observed points, computed once here so the plot subtitles and
# §8's comparison table report the same numbers.
mu_lin_obs <- sapply(x, function(xx) post_lin$a + post_lin$b * xx)
mu_exp_obs <- sapply(x, function(xx) post_exp$a * exp(-post_exp$b * xx))
mu_sp_obs  <- as.numeric(post_sp$a) + post_sp$w %*% t(matrix(as.numeric(Bmat), nrow = N))

y_range <- max(y) - min(y)
nrmse <- function(mu_draws) sqrt(mean((y - apply(mu_draws, 2, mean))^2)) / y_range
nrmse_lin <- nrmse(mu_lin_obs)
nrmse_exp <- nrmse(mu_exp_obs)
nrmse_sp  <- nrmse(mu_sp_obs)

# Posterior-mean equations for the plot subtitles, with slopes/decay rates
# converted from the Mm fitting scale to the DISPLAY_UNIT_KM plotting scale.
disp_per_mm <- DISPLAY_UNIT_KM / 1000  # e.g. 0.1 when DISPLAY_UNIT_KM = 100
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

# subtitle wraps long lines (the spline equation has 6 terms) so it stays
# readable in the saved PNG rather than clipping.
plot_fit <- function(ribbon_df, title, eq, nrmse_val) {
  ggplot() +
    geom_point(data = scatter_df,
               aes(x = geodist_H1_span, y = cossim_span_norm), size = 2) +
    geom_ribbon(data = ribbon_df, aes(x = dist, ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.25) +
    geom_line(data = ribbon_df, aes(x = dist, y = mean),
              color = "steelblue", linewidth = 1.2) +
    theme_bw() +
    labs(title = title,
         subtitle = str_wrap(sprintf("%s   |   nRMSE = %.4f", eq, nrmse_val), width = 70),
         x = "Relative Migration Distance (100 km)", y = "Normalized Spanish Similarity")
}

p_lin <- plot_fit(ribbon_lin, "Linear (phylogenetic)", eq_lin, nrmse_lin)
p_exp <- plot_fit(ribbon_exp, "Exponential decay (phylogenetic)", eq_exp, nrmse_exp)
p_sp  <- plot_fit(ribbon_sp,  "Cubic spline (phylogenetic)", eq_sp, nrmse_sp)
print(p_lin); print(p_exp); print(p_sp)

ggsave(here("figures", "grammar", "regression", "grammar_linear_model.png"),
       p_lin, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "grammar", "regression", "grammar_exponential_model.png"),
       p_exp, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "grammar", "regression", "grammar_spline_model.png"),
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
write.csv(report, file = here("data", "GRAMMAR_waic_compare.csv"), row.names = FALSE)


# ── 9. Collinearity: geographic vs. phylogenetic distance ───────────────────
# Geographic distance from the model predictor, phylogenetic distance from the
# raw VCV of §2 via D_ij = V_ii + V_jj - 2*V_ij.
stopifnot(
  "V_phylo is not square/NxN — rerun from §2 (comp_data build) in a clean session." =
    is.matrix(V_phylo) && nrow(V_phylo) == N && ncol(V_phylo) == N
)
Dgeo   <- as.matrix(dist(x_km))
Dphylo <- outer(diag(V_phylo), diag(V_phylo), "+") - 2 * V_phylo

geo_phylo_cor <- cor(Dgeo[lower.tri(Dgeo)], Dphylo[lower.tri(Dphylo)])
message("cor(geographic distance, phylogenetic distance) = ", round(geo_phylo_cor, 3))
