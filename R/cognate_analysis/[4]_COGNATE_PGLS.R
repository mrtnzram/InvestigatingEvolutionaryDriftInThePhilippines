# =============================================================================
# [4] Cognate Analysis — Bayesian phylogenetic regression
#
# Fits linear, exponential-decay and cubic-spline models of the normalized
# Spanish loanword count vs. terrain-penalized migration distance with ulam(),
# each under a multi_normal likelihood carrying the same Pagel's-lambda-scaled
# phylogenetic covariance, then compares them by WAIC/PSIS-LOO. Lambda is fixed
# at the PGLS ML estimate and shared across models, so WAIC differences reflect
# only the mean structure. Fitting on the normalized scale means a slope reads
# as "fraction of the observed borrowing range per unit distance" rather than a
# raw loan-count delta.
#
# Replaces the independent-normal [4]_COGNATE_OLS.R: shared ancestry is now
# modelled, and the study set is the 94 languages that have a tree tip.
#
# Run order: requires `tree_pruned` and `tip_map` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/cognate/COGNATE_final.csv (from [3]_COGNATE_network_distance.R)
# Outputs: data/cognate/COGNATE_pgls_results.csv  (PGLS lambda / slope table)
#          data/cognate/COGNATE_waic_compare.csv  (WAIC model comparison + nRMSE)
#          figures/cognate/regression/cognate_{linear,exponential,spline}_model.png
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
  "Run [0]_Phylogenetic_Tree.R first: `tip_map` is not defined." =
    exists("tip_map")
)

# geodist_H1_span matches [4]_GENETIC_OLS.R; switch to "geodist_H2_span" to
# regress on distance to Manila instead of distance onto the network.
PREDICTOR <- "geodist_H1_span"

COGNATE_final <- read.csv(here("data", "cognate", "COGNATE_final.csv"))


# ── 1. Comparative dataset + PGLS lambda estimate ───────────────────────────
# comparative.data(vcv = TRUE) returns the VCV already ordered to match its
# $data, and that ordering is reused below so y, x and R_lambda stay in sync.
#
# names.col is the glottocode, not a tree tip label as in the phoneme/grammar
# versions: [0] reduced the tree to one tip per language, so this is one row per
# language with no duplication. [3] already restricted COGNATE_final.csv to the
# tree-matched set, so the semi_join is a no-op there — it is kept so the script
# stays correct against a COGNATE_final.csv regenerated without that filter.
#
# dplyr:: qualified throughout: caper loads MASS, whose select() masks
# dplyr::select() once attached, regardless of library() order.
df <- COGNATE_final |>
  dplyr::select(glottocode, language, loans_norm, geodist_H1_span = all_of(PREDICTOR)) |>
  semi_join(tip_map, by = "glottocode") |>
  filter(!is.na(loans_norm), !is.na(geodist_H1_span)) |>
  as.data.frame()

message("N = ", nrow(df), " languages of ", nrow(COGNATE_final),
        " in COGNATE_final.csv (tree-matched, with a loanword count and a distance).")
stopifnot(
  "Analysis frame and pruned tree disagree on the language set." =
    nrow(df) == length(tree_pruned$tip.label)
)

comp_data <- comparative.data(
  phy = tree_pruned, data = df, names.col = "glottocode", vcv = TRUE
)

# PGLS is the lambda estimator. Fixed-lambda~0 fit is retained only for the
# reporting table (near-OLS reference); pgls_fit (ML) supplies lambda_hat.
pgls_ols <- pgls(loans_norm ~ geodist_H1_span, data = comp_data, lambda = 1e-6)
pgls_fit <- pgls(loans_norm ~ geodist_H1_span, data = comp_data, lambda = "ML")
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

COGNATE_pgls_results <- bind_rows(
  tidy_pgls(pgls_ols, "OLS (lambda fixed ~0)"),
  tidy_pgls(pgls_fit, "ML (lambda estimated)")
)
print(COGNATE_pgls_results, n = Inf)
write.csv(COGNATE_pgls_results,
          file = here("data", "cognate", "COGNATE_pgls_results.csv"), row.names = FALSE)


# ── 2. Build the fixed phylogenetic correlation R_lambda ────────────────────
# R_lambda = lambda*R + (1-lambda)*I on the correlation-scaled VCV; sigma^2 in
# the models scales it up to the residual covariance. The matrix must be stripped
# of its caper "VCV.array" class and dimnames, which break ulam's dimension registration.
ord  <- rownames(comp_data$vcv)
cd   <- comp_data$data[ord, , drop = FALSE]
y    <- cd$loans_norm
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
post_lin <- extract.samples(m_lin)
precis(m_lin, pars = c("a", "b", "sigma"))
cat("equation:", eq_linear(post_lin), "\n")


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
post_exp <- extract.samples(m_exp)
precis(m_exp, pars = c("a", "b", "sigma"))
cat("equation:", eq_exponential(post_exp), "\n")


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
post_sp <- extract.samples(m_spline)
precis(m_spline, pars = c("a", "sigma"))
cat("equation:", eq_spline(post_sp), "\n")


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

# post_lin / post_exp / post_sp were extracted alongside each fit above.
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
# The ribbon is the 95% quantile interval of mu across posterior draws;
# DISPLAY_UNIT_KM (set in §2) rescales the x axis to hundreds of km.
xseq_km <- seq(min(x_km), max(x_km), length.out = 100)
xseq    <- xseq_km / 1000

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

# Fixed point size, matching the [7] map: with no standard error there is
# nothing for point area to encode.
scatter_df <- data.frame(dist = x_km / DISPLAY_UNIT_KM, loans_norm = y)

# nRMSE at the observed points, for §8's comparison table.
mu_lin_obs <- sapply(x, function(xx) post_lin$a + post_lin$b * xx)
mu_exp_obs <- sapply(x, function(xx) post_exp$a * exp(-post_exp$b * xx))
mu_sp_obs  <- as.numeric(post_sp$a) + post_sp$w %*% t(matrix(as.numeric(Bmat), nrow = N))

y_range <- max(y) - min(y)
nrmse <- function(mu_draws) sqrt(mean((y - apply(mu_draws, 2, mean))^2)) / y_range
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

p_lin <- plot_fit(ribbon_lin, "Linear (phylogenetic)")
p_exp <- plot_fit(ribbon_exp, "Exponential decay (phylogenetic)")
p_sp  <- plot_fit(ribbon_sp,  "Cubic spline (phylogenetic)")
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
