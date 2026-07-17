# =============================================================================
# [1.5] Phoneme Analysis — Cosine-similarity distribution analysis
# Visualizes the distribution of each Philippine language's phoneme cosine
# similarity to Spanish / English / Japanese / unrelated controls (ridge +
# density plots), then tests whether the baselines differ (Friedman + linear
# mixed model).
# Input: data/PHONEME_cossim.csv   (from [1]_PHONEME_cosine_similarity.R)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggridges)
library(patchwork)
library(lme4)
library(lmerTest)
library(here)

# Loading Dataframe ------
PHONEME_cossim <- read_csv(here('data','PHONEME_cossim.csv'))
PHONEME_cossim <- PHONEME_cossim[, -1]

# --- Distribution plots ------------------------------------------------------
# Long form keyed by baseline, using friendly labels for the legend/axes.
combined_scores <- PHONEME_cossim |>
  select(Spanish   = cossim_span,
         Japanese  = cossim_jap,
         English   = cossim_eng,
         Unrelated = cossim_unr) |>
  pivot_longer(cols = everything(),
               names_to = "Language",
               values_to = "Similarity_Score")

combined_scores_summary <- combined_scores |>
  group_by(Language) |>
  summarize(mean_score = mean(Similarity_Score), .groups = "drop")

cossim_phoneme_density_ridge <- ggplot(combined_scores, aes(x = Similarity_Score, y = Language, fill = Language)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2, color = "black") +
  # geom_segment keeps each mean line contained within its own ridge
  geom_segment(
    data = combined_scores_summary,
    aes(
      x = mean_score,
      xend = mean_score,
      y = as.numeric(factor(Language)),
      yend = as.numeric(factor(Language)) + 0.9,
      color = Language
    ),
    linetype = "dashed",
    linewidth = 1.2,
    inherit.aes = FALSE
  ) +
  labs(
    title = "Phoneme Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Language"
  ) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  theme(legend.position = "none")

cossim_phoneme_density_ridge

ggsave(
  filename = here("figures", "phoneme", "distributions", "phoneme_ridgeplot.png"),
  plot = cossim_phoneme_density_ridge,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300
)

# Individual baseline-vs-unrelated density comparisons
density_compare <- function(langs) {
  ggplot(combined_scores |> filter(Language %in% langs),
         aes(x = Similarity_Score, fill = Language)) +
    geom_density(alpha = 0.5) +
    geom_vline(
      data = combined_scores_summary |> filter(Language %in% langs),
      aes(xintercept = mean_score, color = Language),
      linetype = "dashed",
      linewidth = 1.2
    ) +
    labs(title = "Cosine Similarity Distribution",
         x = "Similarity Score",
         y = "Density") +
    theme_bw() +
    scale_x_continuous(breaks = seq(0, 0.4, by = 0.02))
}

phoneme_cos_s <- density_compare(c("Unrelated", "Spanish"))
phoneme_cos_e <- density_compare(c("Unrelated", "English"))
phoneme_cos_j <- density_compare(c("Unrelated", "Japanese"))

phoneme_cos_s + phoneme_cos_e + phoneme_cos_j

# --- Frieddman Test --------
sim_matrix <- PHONEME_cossim |>
  select(cossim_span, cossim_eng, cossim_jap, cossim_unr) |>
  as.matrix()

friedman_result <- friedman.test(sim_matrix)
W_kendall       <- friedman_result$statistic / (nrow(sim_matrix) * (ncol(sim_matrix) - 1))

cat(sprintf("Friedman: chi2(3) = %.3f, p = %.4f, W = %.3f\n",
            friedman_result$statistic, friedman_result$p.value, W_kendall))

# no siginificant difference between the distirbutions
# no need for Wilcoxon signed rank test

# --- Linear Mixed Model ------

phoneme_long <- PHONEME_cossim |> 
  pivot_longer(
    cols = c(cossim_span,cossim_eng,cossim_jap,cossim_unr),
    names_to = 'baseline',
    values_to = 'similarity'
  ) |>
  mutate(
    baseline = factor(baseline,
                      levels = c("cossim_unr", "cossim_span", "cossim_eng", "cossim_jap"),
                      labels = c("Unrelated", "Spanish", "English", "Japanese"))
  )


phoneme_long

m_lmm <- lmer(similarity ~ baseline + (1 | language),
              data = phoneme_long, REML = TRUE)

summary(m_lmm)
confint(m_lmm, method = "profile")

vc      <- as.data.frame(VarCorr(m_lmm))
icc     <- vc$vcov[1] / sum(vc$vcov)
sigma_r <- sigma(m_lmm)
d_vec   <- fixef(m_lmm)[-1] / sigma_r

cat(sprintf("ICC = %.3f\n", icc))
cat("Cohen's d per contrast:\n"); print(round(d_vec, 3))

# BH correction for p-values

lmm_coefs <- summary(m_lmm)$coefficients

# Rows 2-4 are the three contrasts (row 1 is the intercept)
p_raw <- lmm_coefs[-1, "Pr(>|t|)"]
p_bh  <- p.adjust(p_raw, method = "BH")

tibble(
  contrast   = names(p_raw),
  beta       = lmm_coefs[-1, "Estimate"],
  cohens_d   = c(0.780, -0.072, 0.522),
  p_raw      = p_raw,
  p_bh       = p_bh,
  reject_bh  = p_bh < 0.05
)

