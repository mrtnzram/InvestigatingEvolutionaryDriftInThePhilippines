library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(here)

# Loading Dataframe ------
GRAMMAR_cossim <- read_csv(here('data','GRAMMAR_cossim.csv'))
GRAMMAR_cossim <- GRAMMAR_cossim[, -1]

# --- Freidman Test --------
sim_matrix <- GRAMMAR_cossim |>
  select(cossim_span, cossim_eng, cossim_jap, cossim_unr) |>
  as.matrix()

friedman_result <- friedman.test(sim_matrix)
W_kendall       <- friedman_result$statistic / (nrow(sim_matrix) * (ncol(sim_matrix) - 1))

cat(sprintf("Friedman: chi2(3) = %.3f, p = %.4f, W = %.3f\n",
            friedman_result$statistic, friedman_result$p.value, W_kendall))

# results are significant at X2 = 96.446 and p-value < 0.00001

# ---- Wilcoxon post-hoc -------
wt_sp_un <- wilcox.test(GRAMMAR_cossim$cossim_span, GRAMMAR_cossim$cossim_unr,
                        paired = TRUE, alternative = "greater", exact = TRUE)
wt_en_un <- wilcox.test(GRAMMAR_cossim$cossim_eng,  GRAMMAR_cossim$cossim_unr,
                        paired = TRUE, alternative = "greater", exact = TRUE)
wt_sp_ja <- wilcox.test(GRAMMAR_cossim$cossim_span, GRAMMAR_cossim$cossim_jap,
                        paired = TRUE, alternative = "greater", exact = TRUE)

rank_biserial <- function(W_plus, n) {
  W_minus <- n * (n + 1) / 2 - W_plus
  1 - (4 * W_minus) / (n * (n + 1))
}

n_pairs   <- nrow(GRAMMAR_cossim)
pvals_raw <- c(wt_sp_un$p.value, wt_en_un$p.value, wt_sp_ja$p.value)
pvals_bh  <- p.adjust(pvals_raw, method = "BH")

results_wilcox <- tibble(
  contrast   = c("Spanish vs Unrelated", "English vs Unrelated", "Spanish vs Japanese"),
  W_plus     = c(wt_sp_un$statistic, wt_en_un$statistic, wt_sp_ja$statistic),
  p_raw      = pvals_raw,
  p_bh       = pvals_bh,
  r_biserial = c(rank_biserial(wt_sp_un$statistic, n_pairs),
                 rank_biserial(wt_en_un$statistic, n_pairs),
                 rank_biserial(wt_sp_ja$statistic, n_pairs)),
  reject_bh  = pvals_bh < 0.05
)

print(results_wilcox)

# --- Linear Mixed Model ------

grammar_long <- GRAMMAR_cossim |> 
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


grammar_long

m_lmm <- lmer(similarity ~ baseline + (1 | language),
              data = grammar_long, REML = TRUE)

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

