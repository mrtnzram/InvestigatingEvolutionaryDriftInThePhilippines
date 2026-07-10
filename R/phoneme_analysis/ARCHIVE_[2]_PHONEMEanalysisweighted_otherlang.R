# =============================================================================
# [2] Phoneme Analysis — Regression & distance matrices for a comparison language
# Continuation of [1] (uses its in-memory objects: PHONEME_cossim, cosine_matrix,
# PHOIBLEdf, nodes, graph, land_sf, find_nearest_node, shortest_path_trace).
# Fits linear/exponential/spline distance-decay models to Japanese similarity and
# builds the phoneme distance/dissimilarity matrices. Reads/writes no files.
# Run after sourcing [1]_PHONEMEanalysisweighted_span.R in the same session.
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(geosphere)
library(rethinking)
library(lingtypology)
library(infotheo)
library(reshape2)
library(ggplot2)
library(purrr)
library(patchwork)
library(igraph)
library(sfheaders)
library(splines)

# ---- 1. Linear model: Japanese similarity vs. migration distance ----
lm_jap <- lm(cossim_jap ~ geodist_H1_span, data = PHONEME_cossim)
summary_lm_jap <- summary(lm_jap)

slope <- round(coef(lm_jap)[2],5)
coefficient <- round(coef(lm_jap)[1],5)
r_squared <- round(summary_lm_jap$r.squared,3)
linear_model_equation <- paste0("Y = ", slope, "x + ", coefficient)

ggplot(data = PHONEME_cossim, aes(x = geodist_H1_span, y = cossim_jap)) +
  geom_point() +
  geom_smooth(method = 'lm',se = FALSE) +
  theme_bw() +
  labs(title = 'Linear Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity')

# ---- 2. Exponential decay model (rethinking::map) with prediction ribbon ----
PHONEME_cossim$dist_std <- standardize(PHONEME_cossim$geodist_H1_span)

m_exp <- rethinking::map(
  alist(
    cossim_jap ~ dnorm(mu, sigma),
    mu <- a * exp(-b * dist_std),
    a ~ dnorm(0.5, 0.2),
    b ~ dnorm(0, 1),
    sigma ~ dexp(1)
  ),
  data = PHONEME_cossim
)

precis(m_exp)

dist_seq <- seq(from = min(PHONEME_cossim$dist_std), to = max(PHONEME_cossim$dist_std), length.out = 100)
preds <- link(m_exp, data = data.frame(dist_std = dist_seq))
mu_mean <- apply(preds, 2, mean)
mu_PI <- apply(preds, 2, PI)
mu_PI_t <- t(mu_PI)
scatter_df <- PHONEME_cossim

ribbon_df <- data.frame(
  dist = dist_seq,
  mean = mu_mean,
  lower = mu_PI_t[,1],
  upper = mu_PI_t[,2]
)

ggplot() +
  geom_point(data = scatter_df, aes(x = dist_std, y = cossim_jap),
             color = "black", size = 2) +
  geom_line(data = ribbon_df, aes(x = dist, y = mean),
            color = "blue", linewidth = 1.2) +
  labs(title = 'Exponential Model',
       x = 'Relative Migration Distance (km)',
       y = 'Cosine Similarity') +
  theme_bw()

coefs <- coef(m_exp)
a <- round(coefs["a"], 5)
b <- round(coefs["b"], 5)

exp_eq <- paste0("μ = ", a, " * exp(-", b, " * x)")
print(exp_eq)

# ---- 3. Spline regression and its knot vector ----
spline_fit <- lm(cossim_jap ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)
x_vals <- PHONEME_cossim$geodist_H1_span
y_spline <- predict(spline_fit, newdata = data.frame(geodist_H1_span = x_vals))

ggplot(PHONEME_cossim, aes(x = geodist_H1_span, y = cossim_jap)) +
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

# Extract the spline basis object and its knots
spline_basis_object <- model.frame(spline_fit)$'bs(geodist_H1_span, df = 5)'
knots <- attr(spline_basis_object, "knots")
print(knots)

# ---- 4. Compare models by RMSE / significance ----
spline_fit <- lm(cossim_jap ~ bs(geodist_H1_span, df = 5), data = PHONEME_cossim)

PHONEME_cossim$pred_linear <- predict(lm_jap)

exp_preds <- link(m_exp)
PHONEME_cossim$pred_exp <- apply(exp_preds, 2, mean)

rmse_linear <- sqrt(mean((PHONEME_cossim$cossim_jap - PHONEME_cossim$pred_linear)^2))
rmse_exp <- sqrt(mean((PHONEME_cossim$cossim_jap - PHONEME_cossim$pred_exp)^2))
rmse_spline <- sqrt(mean((PHONEME_cossim$cossim_jap - y_spline)^2))

p_linear <- summary(lm_jap)$coefficients[2, 4]
p_spline <- pf(summary(spline_fit)$fstatistic[1],
               summary(spline_fit)$fstatistic[2],
               summary(spline_fit)$fstatistic[3],
               lower.tail = FALSE)

exp_post <- extract.samples(m_exp)
pr_b_gt_0 <- mean(exp_post$b > 0)

cossim_jap_range <- max(PHONEME_cossim$cossim_jap) - min(PHONEME_cossim$cossim_jap)

rmse_df_jap <- data.frame(
  Model = c("Linear", "Exponential", "Spline"),
  RMSE = round(c(rmse_linear, rmse_exp, rmse_spline), 5),
  normalized_rmse = round(c(rmse_linear, rmse_exp, rmse_spline) / cossim_jap_range, 5),
  p_value = c(round(p_linear, 5), NA, round(p_spline, 5)),
  posterior_prob = c(NA, round(pr_b_gt_0, 3), NA)
)

print(rmse_df_jap)

# ---- 5. Phoneme dissimilarity matrix (1 - cosine) for Philippine languages ----
cosine_matrix_phil <- cosine_matrix[ph_lang, ph_lang]
lang_order <- rownames(cosine_matrix_phil)

cosine_matrix_phil <- 1-cosine_matrix_phil

PHONEME_diss_matrix <- cosine_matrix_phil

# ---- 6. Phoneme distance matrix: land-penalized geodesic network distances ----
ph_lang <- PHOIBLEdf %>%
  filter(Language_type == 'Philippine Language') %>%
  pull(language)

phil_df <- PHOIBLEdf %>%
  filter(Language_type == "Philippine Language") %>%
  mutate(
    start_coords = map2(longitude, latitude, ~ c(.x, .y)),
    nearest_node = map_chr(start_coords, find_nearest_node)
  )

land_penalty <- 4.44

connector_df <- phil_df %>%
  mutate(
    connector_geom = map2(start_coords, nearest_node, ~ st_linestring(rbind(
      .x,
      c(nodes$longitude[nodes$id == .y], nodes$latitude[nodes$id == .y])
    )))
  )

connector_df <- connector_df %>%
  mutate(connector_geom_sfc = st_sfc(connector_geom, crs = 4326))

connector_df <- connector_df %>%
  rowwise() %>%
  mutate(
    land_part = list(st_intersection(connector_geom_sfc, land_sf)),
    sea_part  = list(st_difference(connector_geom_sfc, land_sf)),

    land_len = as.numeric(if (!is.null(land_part) && length(land_part) > 0) st_length(land_part) else 0),
    sea_len  = as.numeric(if (!is.null(sea_part)  && length(sea_part)  > 0) st_length(sea_part)  else 0),

    connector_penalty = land_len * land_penalty + sea_len
  ) %>%
  ungroup()

phil_pairs <- expand_grid(lang1 = phil_df$language, lang2 = phil_df$language) %>%
  filter(lang1 != lang2) %>%
  left_join(phil_df %>% select(language, node1 = nearest_node), by = c("lang1" = "language")) %>%
  left_join(phil_df %>% select(language, node2 = nearest_node), by = c("lang2" = "language")) %>%
  left_join(connector_df %>% select(language, penalty1 = connector_penalty), by = c("lang1" = "language")) %>%
  left_join(connector_df %>% select(language, penalty2 = connector_penalty), by = c("lang2" = "language"))

phil_pairs <- phil_pairs %>%
  rowwise() %>%
  mutate(
    trace = list(shortest_path_trace(node1, node2, graph)),
    tree_dist = trace$distance,
    geodist_H1_span = if (is.na(tree_dist)) NA_real_ else
      (penalty1 + tree_dist + penalty2) / 1000
  ) %>%
  ungroup()

dist_matrix <- phil_pairs %>%
  select(lang1, lang2, geodist_H1_span) %>%
  pivot_wider(names_from = lang2, values_from = geodist_H1_span) %>%
  column_to_rownames("lang1") %>%
  as.matrix()

PHONEME_dist_matrix <- dist_matrix[ph_lang, ph_lang]
PHONEME_dist_matrix[is.na(PHONEME_dist_matrix)] <- 0

# ---- 7. Plot the distance and dissimilarity matrices side by side ----
melt_phoneme_dist_matrix <- melt(PHONEME_dist_matrix)
melt_phoneme_diss_matrix <- melt(PHONEME_diss_matrix)

dist_matrix_p <- ggplot(melt_phoneme_dist_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

diss_matrix_p  <- ggplot(melt_phoneme_diss_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_fixed()

dist_matrix_p + diss_matrix_p
