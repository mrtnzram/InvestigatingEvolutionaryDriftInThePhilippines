library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(here)

# Loading Dataframe ------
PHONEME_cossim <- read_csv(here('data','PHONEME_cossim.csv'))
PHONEME_cossim <- PHONEME_cossim[, -1]

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

