# =============================================================================
# [4] Phoneme Analysis — Distance-decay regression
# Fits linear, exponential, and spline models of Spanish phoneme cosine
# similarity against relative migration distance, then compares their fit
# (RMSE / p-value / posterior probability).
# Input: data/PHONEME_cossim_dist.csv  (from [3]_PHONEME_network_distance.R)
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(rethinking)
library(splines)
library(here)

PHONEME_cossim_dist <- read.csv(here("data", "PHONEME_cossim_dist.csv"))

# ----- linear model ----------------------------------------------------------

lm_span <- lm(cossim_span ~ geodist_H1_span, data = PHONEME_cossim_dist)
summary_lm_span <- summary(lm_span)

slope <- round(coef(lm_span)[2], 5)
coefficient <- round(coef(lm_span)[1], 5)
r_squared <- round(summary_lm_span$r.squared, 3)
linear_model_equation <- paste0("Y = ", slope, "x + ", coefficient)
print(linear_model_equation)

ggplot(data = PHONEME_cossim_dist, aes(x = geodist_H1_span, y = delta_span)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  theme_bw() +
  labs(title = 'Linear Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity')

# ---- exponential model ------------------------------------------------------

# Standardize predictors
PHONEME_cossim_dist$dist_std <- rethinking::standardize(PHONEME_cossim_dist$geodist_H1_span)

# Model: exponential decay
m_exp <- rethinking::map(
  alist(
    cossim_span ~ dnorm(mu, sigma),
    mu <- a * exp(-b * geodist_H1_span),
    a ~ dnorm(0.5, 0.2),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = PHONEME_cossim_dist
)

# Summarize
precis(m_exp)

# Generate predictions
dist_seq <- seq(from = min(PHONEME_cossim_filter_test$geodist_H1_span), to = max(PHONEME_cossim_filter_test$geodist_H1_span), length.out = 100)
preds <- link(m_exp, data = data.frame(geodist_H1_span = dist_seq))
mu_mean <- apply(preds, 2, mean)
mu_PI <- apply(preds, 2, PI)
mu_PI_t <- t(mu_PI)
# Plot
scatter_df <- PHONEME_cossim_filter_test

# Line + ribbon data
ribbon_df <- data.frame(
  dist = dist_seq,
  mean = mu_mean,
  lower = mu_PI_t[, 1],
  upper = mu_PI_t[, 2]
)

ggplot() +
  geom_point(data = scatter_df, aes(x = geodist_H1_span, y = delta_spanish),
             color = "black", size = 2) +
  geom_line(data = ribbon_df, aes(x = dist, y = mean),
            color = "blue", linewidth = 1.2) +
  labs(title = 'Exponential Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity') +
  theme_bw()

# Extract coefficients
coefs <- coef(m_exp)
a <- round(coefs["a"], 5)
b <- round(coefs["b"], 5)

# Write out the equation
exp_eq <- paste0("μ = ", a, " * exp(-", b, " * x)")
print(exp_eq)

# ----- spline regression -----------------------------------------------------

spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)
x_vals <- PHONEME_cossim$geodist_H1_span
y_spline <- predict(spline_fit, newdata = data.frame(geodist_H1_span = x_vals))

ggplot(PHONEME_cossim, aes(x = geodist_H1_span, y = cossim_span)) +
  geom_point(color = "gray") +
  geom_line(aes(y = y_spline), color = "blue", linewidth = 1.2) +
  theme_bw() +
  labs(title = 'Spline Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity')

paste0("y(x) = ", paste0("β", 1:5, "·B", 1:5, "(x)", collapse = " + "))
coefs <- coef(spline_fit)

paste0(
  "y(x) = ", round(coefs[1], 5), " + ",
  paste0(round(coefs[-1], 5), "·B", 1:5, "(x)", collapse = " + ")
)

# 1. Extract the spline basis object from the model frame
spline_basis_object <- model.frame(spline_fit)$'bs(geodist_H1_span, df = 5)'

# 2. Extract the 'knots' attribute from that object
knots <- attr(spline_basis_object, "knots")

# Print the resulting knot vector
print(knots)

# ----- comparing errors ------------------------------------------------------

spline_fit <- lm(cossim_span ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)

# Linear model predictions
PHONEME_cossim$pred_linear <- predict(lm_span)

# Exponential model predictions using link()
exp_preds <- link(m_exp)
PHONEME_cossim$pred_exp <- apply(exp_preds, 2, mean)

# RMSE calculation
rmse_linear <- sqrt(mean((PHONEME_cossim$cossim_span - PHONEME_cossim$pred_linear)^2))
rmse_exp <- sqrt(mean((PHONEME_cossim$cossim_span - PHONEME_cossim$pred_exp)^2))
rmse_spline <- sqrt(mean((PHONEME_cossim$cossim_span - y_spline)^2))

p_linear <- summary(lm_span)$coefficients[2, 4]  # Slope p-value
p_spline <- pf(summary(spline_fit)$fstatistic[1],
               summary(spline_fit)$fstatistic[2],
               summary(spline_fit)$fstatistic[3],
               lower.tail = FALSE)

exp_post <- extract.samples(m_exp)
pr_b_gt_0 <- mean(exp_post$b > 0)  # Replace 'b' with your actual slope parameter name

# Compare

cossim_span_range <- max(PHONEME_cossim$cossim_span) - min(PHONEME_cossim$cossim_span)

rmse_df <- data.frame(
  Model = c("Linear", "Exponential", "Spline"),
  RMSE = round(c(rmse_linear, rmse_exp, rmse_spline), 5),
  normalized_rmse = round(c(rmse_linear, rmse_exp, rmse_spline) / cossim_span_range, 5),
  p_value = c(round(p_linear, 5), NA, round(p_spline, 5)),
  posterior_prob = c(NA, round(pr_b_gt_0, 3), NA)
)

print(rmse_df)
