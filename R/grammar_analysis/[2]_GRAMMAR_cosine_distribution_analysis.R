# =============================================================================
# [2] Grammar Analysis — Cosine-similarity distribution analysis
# Visualizes the distribution of each Philippine language's grammar cosine
# similarity to Spanish / English / Japanese / unrelated controls (ridge +
# density plots), then tests whether the baselines stand out from the unrelated
# controls at two levels:
#
#   Population level — Friedman test, Wilcoxon signed-rank test (+ box plot),
#                      and a BIC-selected Gaussian mixture (bimodality) test that
#                      marks each language as <baseline>_influenced.
#   Individual level — Shapiro-Wilk normality of each language's per-language
#                      unrelated null (+ QQ plots and a null-ridge with the
#                      observed baselines marked), and an empirical percentile
#                      of each language against its own null.
#
# Input:   data/GRAMMAR_cossim.csv         (from [1]_GRAMMAR_cosine_similarity.R)
#          data/GRAMMAR_cosine_matrix.csv  (per-language null: cols = unrelated)
#          data/GRAMBANKdf_full.csv        (language-type lookup for ph/unr sets)
# Outputs: data/GRAMMAR_cossim_marked.csv  (GRAMMAR_cossim + *_influenced / sig_*
#                                            columns for [3] to carry downstream)
# Next:    [3]_GRAMMAR_network_distance.R
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggridges)
library(patchwork)
library(here)
library(broom)
library(mclust)
library(lme4)
library(lmerTest)

# --- Loading data ------------------------------------------------------------
# read_csv names an unnamed leading index column "...1"; drop it if present.
# Also drop any *_influenced / sig_* columns that a previous run may have left in
# the file, so the classification joins below stay idempotent (no .x/.y dupes).
GRAMMAR_cossim <- read_csv(here("data", "GRAMMAR_cossim.csv")) |>
  dplyr::select(-any_of(c("...1", "span_influenced", "jap_influenced", "eng_influenced",
                   "sig_span", "sig_jap", "sig_eng")))

# Philippine / unrelated language sets (derived locally so [2] is self-contained
# and does not depend on [1] having been sourced in the same session).
GRAMBANKdf_PH <- read_csv(here("data", "GRAMBANKdf_full.csv"))
ph_lang  <- GRAMBANKdf_PH |> filter(Language_Type == "Philippine Language") |> pull(language)
unr_lang <- GRAMBANKdf_PH |> filter(Language_Type == "Unrelated Language")  |> pull(language)

# Per-language unrelated "null" matrix: rows = Philippine languages, cols = the
# individual unrelated controls (each language's null distribution of cosine
# similarities). Used by the bimodality and individual-level analyses below.
load_null_matrix <- function(ph_lang, unr_lang) {
  read.csv(here("data", "GRAMMAR_cosine_matrix.csv"),
           row.names = 1, check.names = FALSE) |>
    as.matrix() |>
    (\(m) m[ph_lang, unr_lang, drop = FALSE])()
}
null_mat <- load_null_matrix(ph_lang, unr_lang)

# Per-language excess similarity over the (mean) unrelated baseline. cossim_unr
# is already rowMeans(null_mat) by construction (see [1]), so reuse it directly.
delta_df <- GRAMMAR_cossim |>
  transmute(language,
            delta_span = cossim_span - cossim_unr,
            delta_jap  = cossim_jap  - cossim_unr,
            delta_eng  = cossim_eng  - cossim_unr)


# --- Distribution overview plots ---------------------------------------------
# Long form keyed by baseline, using friendly labels for the legend/axes.
combined_scores <- GRAMMAR_cossim |>
  select(Spanish   = cossim_span,
         Japanese  = cossim_jap,
         English   = cossim_eng,
         Unrelated = cossim_unr) |>
  pivot_longer(cols = everything(),
               names_to = "Language",
               values_to = "Similarity_Score")

combined_scores_summary <- combined_scores |>
  group_by(Language) |>
  summarize(mean_score = mean(Similarity_Score),
            median_score = median(Similarity_Score),
            .groups = "drop")

cossim_grammar_density_ridge <- ggplot(combined_scores, aes(x = Similarity_Score, y = Language, fill = Language)) +
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
    title = "Grammar Cosine Similarity Distribution",
    x = "Similarity Score",
    y = "Language"
  ) +
  theme_minimal() +
  # Extend the axis to the full similarity range (grammar English similarity
  # reaches ~0.78, so a 0..0.5 axis clipped the labelled ticks). Labelled major
  # ticks every 0.1; minor ticks every 0.05. ceiling() rounds the top up to the
  # next 0.05 so the last ridge's tail still sits under a tick.
  (\(xmax) scale_x_continuous(
    breaks       = seq(0, xmax, by = 0.10),
    minor_breaks = seq(0, xmax, by = 0.05),
    guide        = guide_axis(minor.ticks = TRUE)
  ))(ceiling(max(combined_scores$Similarity_Score) / 0.05) * 0.05) +
  scale_y_discrete(expand = c(0.01, 0)) +
  # theme_minimal() blanks axis ticks; re-enable so the minor ticks render.
  theme(legend.position = "none",
        axis.ticks.x = element_line(linewidth = 0.3, colour = "grey40"),
        axis.minor.ticks.x.bottom = element_line(linewidth = 0.3, colour = "grey70"))

cossim_grammar_density_ridge

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_ridgeplot.png"),
  plot = cossim_grammar_density_ridge,
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

grammar_cos_s <- density_compare(c("Unrelated", "Spanish"))
grammar_cos_e <- density_compare(c("Unrelated", "English"))
grammar_cos_j <- density_compare(c("Unrelated", "Japanese"))

grammar_cos_s + grammar_cos_e + grammar_cos_j


# =============================================================================
# Population Level
# =============================================================================

# --- Friedman test -----------------------------------------------------------
# Do the four baselines differ across languages (repeated-measures, one row per
# language)? Accompanied by Kendall's W as an effect size.
sim_matrix <- GRAMMAR_cossim |>
  select(cossim_span, cossim_eng, cossim_jap, cossim_unr) |>
  as.matrix()

friedman_result <- friedman.test(sim_matrix)
W_kendall       <- friedman_result$statistic / (nrow(sim_matrix) * (ncol(sim_matrix) - 1))

friedman_tbl <- tibble(
  test      = "Friedman",
  chi_sq    = unname(friedman_result$statistic),
  df        = unname(friedman_result$parameter),
  p_value   = friedman_result$p.value,
  W_kendall = unname(W_kendall)
)

print(friedman_tbl)

# --- Wilcoxon signed-rank test -----------------------------------------------
# One-sample signed-rank on the per-language delta (colonial − unrelated), one
# test per baseline, BH-corrected across the three. r_rb = rank-biserial effect.
delta_long <- delta_df |>
  pivot_longer(
    cols      = c(delta_span, delta_jap, delta_eng),
    names_to  = "baseline",
    values_to = "delta"
  ) |>
  mutate(baseline = recode(baseline,
                           delta_span = "Spanish",
                           delta_jap  = "Japanese",
                           delta_eng  = "English"))

wilcox_results <- delta_long |>
  summarise(
    test      = list(wilcox.test(delta, mu = 0, alternative = "greater")),
    n_nonzero = sum(delta != 0),
    .by = baseline
  ) |>
  mutate(
    V       = map_dbl(test, \(t) unname(t$statistic)),
    p_value = map_dbl(test, \(t) t$p.value),
    r_rb    = (V / (n_nonzero * (n_nonzero + 1) / 2)) * 2 - 1,
    p_adj   = p.adjust(p_value, method = "BH")   # across the 3 baselines
  ) |>
  select(baseline, n_nonzero, V, r_rb, p_value, p_adj) |>
  arrange(p_value)

print(wilcox_results)

# --- Wilcoxon box plot -------------------------------------------------------
# order ridges by median excess (largest at top)
baseline_order <- delta_long |>
  summarise(m = median(delta), .by = baseline) |>
  arrange(m) |>
  pull(baseline)

delta_long <- delta_long |>
  mutate(baseline = factor(baseline, levels = baseline_order))

# significance labels from the Wilcoxon results
sig_labels <- wilcox_results |>
  mutate(
    baseline = factor(baseline, levels = baseline_order),
    label = paste0("V = ", V,
                   ", p ", ifelse(p_adj < 0.001, "< 0.001",
                                  paste0("= ", round(p_adj, 3))),
                   ", r = ", round(r_rb, 2))
  )

wilcox_boxplot <- ggplot(delta_long, aes(x = delta, y = baseline, color = baseline)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_boxplot(aes(fill = baseline), alpha = 0.25, outlier.shape = NA,
               width = 0.5, color = "grey40") +
  geom_jitter(height = 0.12, alpha = 0.6, size = 1.4) +
  geom_text(data = sig_labels, aes(x = Inf, y = baseline, label = label),
            hjust = 1.1, vjust = -2.5, size = 4, color = "grey20",
            inherit.aes = FALSE) +
  scale_color_manual(values = c(Spanish = "#2ca6a4", Japanese = "#8fb339",
                                English = "#e07a5f")) +
  scale_fill_manual(values = c(Spanish = "#2ca6a4", Japanese = "#8fb339",
                               English = "#e07a5f")) +
  labs(title = "Delta similarity over the unrelated baseline",
       x = expression(delta ~ "= colonial similarity − unrelated baseline"),
       y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

wilcox_boxplot

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_wilcoxon_boxplot.png"),
  plot = wilcox_boxplot,
  width = 8, height = 3.5, units = "in", dpi = 300
)

# --- Linear Mixed Model ------------------------------------------------------
# Baseline as fixed effect (Unrelated = reference), language as random intercept.
# lmerTest supplies Kenward-Roger df for the t-tests.
grammar_long <- GRAMMAR_cossim |>
  select(language, cossim_span, cossim_eng, cossim_jap, cossim_unr) |>
  pivot_longer(
    cols      = c(cossim_span, cossim_eng, cossim_jap, cossim_unr),
    names_to  = "baseline",
    values_to = "similarity"
  ) |>
  mutate(
    baseline = factor(baseline,
                      levels = c("cossim_unr", "cossim_span", "cossim_eng", "cossim_jap"),
                      labels = c("Unrelated", "Spanish", "English", "Japanese"))
  )

m_lmm <- lmer(similarity ~ baseline + (1 | language),
              data = grammar_long, REML = TRUE)

# Variance components and ICC
vc       <- as.data.frame(VarCorr(m_lmm))
sigma_u2 <- vc$vcov[1]
sigma2   <- vc$vcov[2]

lmm_icc_tbl <- tibble(
  sigma_u2     = sigma_u2,                          # between-unit-of-observation variance
  sigma2       = sigma2,                             # residual variance
  icc          = sigma_u2 / (sigma_u2 + sigma2),    # proportion explained by language identity
  sigma_resid  = sigma(m_lmm)                       # residual SD; denominator for Cohen's d
)

print(lmm_icc_tbl)

# Profile CIs (asymmetric, more accurate than Wald at n = 14; rows matched to contrasts)
ci_mat <- confint(m_lmm, method = "profile", quiet = TRUE) |>
  as_tibble(rownames = "term") |>
  filter(str_detect(term, "^baseline")) |>
  rename(ci_lo = `2.5 %`, ci_hi = `97.5 %`)

# Fixed effects table: contrasts only (drop intercept)
coef_mat <- summary(m_lmm)$coefficients |>
  as_tibble(rownames = "term") |>
  filter(term != "(Intercept)")

p_raw_lmm <- coef_mat$`Pr(>|t|)`

lmm_results <- tibble(
  baseline  = str_remove(coef_mat$term, "^baseline"),
  estimate  = coef_mat$Estimate,
  se        = coef_mat$`Std. Error`,
  ci_lo     = ci_mat$ci_lo,
  ci_hi     = ci_mat$ci_hi,
  df        = coef_mat$df,
  t         = coef_mat$`t value`,
  p_raw     = p_raw_lmm,
  p_bh      = p.adjust(p_raw_lmm, method = "BH"),
  cohens_d  = coef_mat$Estimate / sigma(m_lmm),
  reject_bh = p.adjust(p_raw_lmm, method = "BH") < 0.05
) |>
  arrange(p_raw)

print(lmm_results)


ranef_tbl <- ranef(m_lmm)$language |>
  as_tibble(rownames = "language") |>
  rename(u_i = `(Intercept)`) |>
  arrange(u_i) |>
  mutate(language = factor(language, levels = language),
         direction = if_else(u_i >= 0, "above", "below"))

lmm_ranef_plot <- ggplot(ranef_tbl, aes(x = u_i, y = language, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_segment(aes(x = 0, xend = u_i,
                   y = language, yend = language),
               linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = c(above = "#2ca6a4", below = "#e07a5f")) +
  labs(
    title    = "LMM random intercepts per language",
    subtitle = paste0("û_i = language's deviation from population mean across all baselines",
                      "  |  σ_u = ", round(sqrt(lmm_icc_tbl$sigma_u2), 4)),
    x        = "Random intercept (û_i)",
    y        = NULL
  ) +
  theme_bw() +
  theme(legend.position = "none")

lmm_ranef_plot

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_lmm_ranef.png"),
  plot = lmm_ranef_plot,
  width = 10, height = 6, units = "in", dpi = 300
)


# --- Bimodality test (BIC-selected Gaussian mixture) -------------------------
# Fit a 1–3 component Gaussian mixture to each baseline's deltas; BIC picks the
# component count. A k > 1 fit implies a distinct "influenced" sub-population
# (the highest-mean component), which we use to classify each language.
deltas <- list(
  Spanish  = setNames(delta_df$delta_span, delta_df$language),
  Japanese = setNames(delta_df$delta_jap,  delta_df$language),
  English  = setNames(delta_df$delta_eng,  delta_df$language)
)

# ----- fit GMM per baseline, return both a summary row and per-language class -
fit_gmm <- function(delta_vec, baseline_name) {
  mix <- Mclust(delta_vec, G = 1:3, verbose = FALSE)
  k   <- mix$G
  inf_comp <- which.max(mix$parameters$mean)   # highest-δ component = influenced

  # per-language classification (only meaningful if k > 1)
  per_lang <- tibble(
    baseline     = baseline_name,
    language     = names(delta_vec),
    delta        = as.numeric(delta_vec),
    component    = mix$classification,
    p_influenced = if (k > 1) mix$z[, inf_comp] else NA_real_,
    influenced   = if (k > 1) mix$classification == inf_comp else FALSE
  )

  # one-row fit summary
  summary_row <- tibble(
    baseline    = baseline_name,
    k_selected  = k,
    bic         = max(mix$BIC, na.rm = TRUE),
    n_influenced = sum(per_lang$influenced),
    prop_influenced = mean(per_lang$influenced),
    mean_uninfluenced = min(mix$parameters$mean),
    mean_influenced   = max(mix$parameters$mean)
  )

  list(summary = summary_row, per_lang = per_lang)
}

gmm_fits <- imap(deltas, fit_gmm)

# ----- the two tibbles ------------------------------------------------------
# purrr::map explicitly namespaced: library(mclust) is loaded after tidyverse and
# masks purrr::map with mclust::map (classification error), so a bare map() here
# fails with a cryptic "invalid 'length' argument".
gmm_summary        <- purrr::map(gmm_fits, "summary") |> list_rbind()
gmm_classification <- purrr::map(gmm_fits, "per_lang") |> list_rbind()

print(gmm_summary)

# --- Mark languages as <baseline>_influenced ---------------------------------
influenced_wide <- gmm_classification |>
  select(baseline, language, influenced) |>
  pivot_wider(names_from = baseline, values_from = influenced) |>
  rename(span_influenced = Spanish,
         jap_influenced  = Japanese,
         eng_influenced  = English)

GRAMMAR_cossim <- GRAMMAR_cossim |>
  left_join(influenced_wide, by = "language")


# =============================================================================
# Individual Level
# =============================================================================

# --- Shapiro-Wilk normality of each language's null --------------------------
# Test each language's per-language unrelated null (its row of null_mat, ~211
# controls) for normality, then split languages into p-value tertiles used to
# sample representatives for the QQ / ridge diagnostics below.
shapiro_tbl <- tibble(language = rownames(null_mat)) |>
  mutate(
    test    = purrr::map(language, \(lang) shapiro.test(null_mat[lang, ])),
    W       = map_dbl(test, \(t) unname(t$statistic)),
    p_value = map_dbl(test, \(t) t$p.value),
    p_adj   = p.adjust(p_value, method = "BH")
  ) |>
  select(-test) |>
  mutate(
    level = case_when(
      ntile(p_value, 3) == 3 ~ "most normal",
      ntile(p_value, 3) == 2 ~ "mid",
      TRUE                   ~ "least normal"
    ) |> factor(levels = c("most normal", "mid", "least normal"))
  ) |>
  arrange(p_value)

print(shapiro_tbl)
cat(sprintf("Shapiro-Wilk: %d / %d languages reject normality (p_adj < 0.05)\n",
            sum(shapiro_tbl$p_adj < 0.05), nrow(shapiro_tbl)))

# Three representative languages per tertile: spread across the p-value range
# within each level (low / mid / high) so the grid shows variation within tiers.
rep_langs <- shapiro_tbl |>
  arrange(level, p_value) |>
  group_by(level) |>
  slice(round(quantile(seq_len(n()), probs = c(0.25, 0.5, 0.75)))) |>
  ungroup()

sel_langs <- rep_langs$language

# Long null draws for the selected languages (rows of null_mat -> long).
null_sel_long <- as_tibble(null_mat[sel_langs, , drop = FALSE], rownames = "language") |>
  pivot_longer(-language, names_to = "unr", values_to = "similarity") |>
  left_join(rep_langs |> select(language, level), by = "language") |>
  mutate(language = factor(language, levels = sel_langs))

# --- QQ plots: 3 x 3 grid (3 languages per Shapiro tertile) ------------------
qq_plot <- ggplot(null_sel_long, aes(sample = similarity)) +
  stat_qq(size = 0.8, alpha = 0.6) +
  stat_qq_line(color = "#e07a5f") +
  facet_wrap(~ level + language, scales = "free", ncol = 3,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  labs(title = "Grammar individual null distribution QQ plots",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_bw()

qq_plot

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_null_qqplots.png"),
  plot = qq_plot,
  width = 8, height = 3.5, units = "in", dpi = 300
)

# --- Null-distribution ridge for influence-selected languages ----------------
# Three languages influenced by each interest language (Spanish / Japanese /
# English), each shown as its own unrelated null with dashed markers for where
# its observed Spanish / Japanese / English similarities fall — so we can see how
# far into (or beyond) its own null each baseline sits. Grouped on the y-axis by
# which interest language influenced it.
pick_influenced <- function(inf_col, obs_col, influence_name) {
  GRAMMAR_cossim |>
    filter(.data[[inf_col]]) |>
    transmute(language, influence = influence_name, observed = .data[[obs_col]]) |>
    arrange(observed) |>
    # low/mid/high quantile indices; unique() guards small n (e.g. n = 3), where
    # round()'s banker's rounding can otherwise collapse two indices together
    # and duplicate a row (breaking the unique per-language `label` downstream).
    slice(unique(round(quantile(seq_len(n()), probs = c(0.25, 0.5, 0.75)))))
}

null_ridge_langs <- bind_rows(
  pick_influenced("span_influenced", "cossim_span", "Spanish"),
  pick_influenced("jap_influenced",  "cossim_jap",  "Japanese"),
  pick_influenced("eng_influenced",  "cossim_eng",  "English")
) |>
  mutate(influence = factor(influence, levels = c("Spanish", "Japanese", "English")),
         # unique label so a language influenced by >1 baseline still gets its own ridge
         label     = paste0(language, " (", influence, ")")) |>
  arrange(influence, observed)

# insert a blank spacer level between influence groups so their boxes (and the
# ridges, which extend ~1 unit up from each baseline) never overlap
spacers <- c(Spanish = " ", Japanese = "  ")   # unique empty labels
lab_levels <- null_ridge_langs |>
  group_split(influence) |>
  purrr::imap(\(g, i) if (i < 3) c(g$label, spacers[[i]]) else g$label) |>
  unlist()

# per-language null draws, keyed by the unique label
null_ridge_long <- null_ridge_langs |>
  mutate(similarity = purrr::map(language, \(l) as.numeric(null_mat[l, ]))) |>
  select(label, similarity) |>
  unnest_longer(similarity) |>
  mutate(label = factor(label, levels = lab_levels))

# observed Spanish / Japanese / English similarity for each selected language
marker_df <- null_ridge_langs |>
  left_join(GRAMMAR_cossim |> select(language, Spanish = cossim_span,
                                     Japanese = cossim_jap, English = cossim_eng),
            by = "language") |>
  select(label, Spanish, Japanese, English) |>
  pivot_longer(-label, names_to = "baseline", values_to = "observed") |>
  mutate(label = factor(label, levels = lab_levels))

# background band per influence group so the three are visually boxed together,
# using the (spacer-adjusted) y positions of each group's languages
band_palette <- c(Spanish = "#2ca6a4", Japanese = "#8fb339", English = "#e07a5f")
label_pos <- tibble(label = lab_levels, ypos = seq_along(lab_levels))
group_bands <- null_ridge_langs |>
  left_join(label_pos, by = "label") |>
  summarise(ymin = min(ypos) - 0.4, ymax = max(ypos) + 0.95, .by = influence)

null_ridge <- ggplot(null_ridge_long, aes(x = similarity, y = label)) +
  geom_rect(data = group_bands, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = influence),
            alpha = 0.12) +
  geom_density_ridges(fill = "grey80", color = "black", alpha = 0.6, scale = 0.95) +
  # markers constrained to each language's own ridge (same technique as the mean
  # lines in the overview ridge above)
  geom_segment(
    data = marker_df,
    aes(x = observed, xend = observed,
        y = as.numeric(label), yend = as.numeric(label) + 0.9,
        color = baseline),
    linewidth = 0.8, inherit.aes = FALSE, linetype = 'dashed'
  ) +
  # group name at the (empty) right of each band
  geom_text(data = group_bands, inherit.aes = FALSE,
            aes(x = Inf, y = ymax, label = paste0(influence, "-influenced")),
            hjust = 1.05, vjust = 1.4, color = "grey30", fontface = "bold", size = 3.5) +
  scale_color_manual(values = band_palette) +
  scale_fill_manual(values = band_palette, guide = "none") +
  # y labels show just the language (the group is shown by the band instead);
  # drop = FALSE keeps the empty spacer rows that separate the groups
  scale_y_discrete(expand = c(0.01, 0), drop = FALSE,
                   labels = \(x) sub(" \\(.*\\)$", "", x)) +
  labs(title = "Influenced languages' nulls with observed baselines",
       x = "Cosine similarity", y = NULL, color = "Baseline") +
  theme_minimal()

null_ridge

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_null_ridge_markers.png"),
  plot = null_ridge,
  width = 7, height = 4.5, units = "in", dpi = 300
)

# --- Empirical percentile against per-language null --------------------------
# Non-parametric: per-language nulls are frequently non-normal (Shapiro-Wilk
# above), so we use an empirical rank p-value rather than a parametric test.
test_baseline <- function(null_mat, observed_vec, baseline_name, alpha = 0.05) {
  languages <- rownames(null_mat)
  n_unr     <- ncol(null_mat)

  p_raw <- map_dbl(languages, \(lang) {
    u <- null_mat[lang, ]
    y <- observed_vec[[lang]]
    (1 + sum(u >= y)) / (n_unr + 1)         # one-sided upper tail, +1 corrected
  })

  tibble(
    baseline    = baseline_name,
    language    = languages,
    observed    = observed_vec[languages],
    null_mean   = rowMeans(null_mat),
    z_score     = (observed_vec[languages] - rowMeans(null_mat)) /
      apply(null_mat, 1, sd),        # descriptive only
    p_raw       = p_raw,
    p_adj       = p.adjust(p_raw, method = "BH"),
    significant = p.adjust(p_raw, method = "BH") < alpha
  )
}

# ----- run all three baselines ----------------------------------------------
observed <- list(
  Spanish  = setNames(GRAMMAR_cossim$cossim_span, GRAMMAR_cossim$language),
  Japanese = setNames(GRAMMAR_cossim$cossim_jap,  GRAMMAR_cossim$language),
  English  = setNames(GRAMMAR_cossim$cossim_eng,  GRAMMAR_cossim$language)
)

influence_results <- imap(
  observed,
  \(vec, name) test_baseline(null_mat, vec, name)
) |>
  list_rbind()

print(influence_results)

# Percentile ridge (percentile = 1 - empirical upper-tail p)
percentile_df <- influence_results |>
  mutate(percentile = 1 - p_raw)

pct_summary <- percentile_df |>
  summarise(mean_pct = mean(percentile),
            median_pct = median(percentile), .by = baseline) |>
  arrange(mean_pct) |>
  mutate(baseline = factor(baseline, levels = baseline))

percentile_df <- percentile_df |>
  mutate(baseline = factor(baseline, levels = levels(pct_summary$baseline)))

percentile_ridge <- ggplot(percentile_df, aes(x = percentile, y = baseline, fill = baseline)) +
  ggridges::geom_density_ridges(
    scale = 0.9, alpha = 0.5,
    quantile_lines = TRUE, quantiles = 2
  ) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey40", linewidth = 1) +
  scale_fill_manual(values = c(Spanish = "#2ca6a4", Japanese = "#8fb339",
                               English = "#e07a5f")) +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  labs(title = "Percentile against per-language unrelated null",
       x = "Percentile within per-language null", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

percentile_ridge

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_percentile_ridge.png"),
  plot = percentile_ridge,
  width = 7, height = 4.5, units = "in", dpi = 300
)

# --- Influenced languages against their nulls --------------------------------
# One ridge per <baseline>-influenced language: its unrelated null distribution,
# with a dashed line at its observed similarity to that baseline. n = 10 / 20 /
# 10 languages for Spanish / Japanese / English respectively.
plot_influenced_nulls <- function(inf_col, obs_col, baseline_name, fill_color) {
  # influenced languages, ordered by observed similarity (largest at top)
  obs_df <- GRAMMAR_cossim |>
    filter(.data[[inf_col]]) |>
    transmute(language, observed = .data[[obs_col]]) |>
    arrange(observed)

  null_long <- as_tibble(null_mat[obs_df$language, , drop = FALSE], rownames = "language") |>
    pivot_longer(-language, names_to = "unr", values_to = "similarity") |>
    mutate(language = factor(language, levels = obs_df$language))

  obs_df <- obs_df |> mutate(language = factor(language, levels = obs_df$language))

  obs_label <- paste0("Observed ", baseline_name, " similarity")

  ggplot(null_long, aes(x = similarity, y = language)) +
    geom_density_ridges(fill = fill_color, color = "black", alpha = 0.5, scale = 1.1) +
    geom_segment(
      data = obs_df,
      aes(x = observed, xend = observed,
          y = as.numeric(language), yend = as.numeric(language) + 0.9,
          linetype = obs_label),
      color = "black", linewidth = 0.8
    ) +
    scale_linetype_manual(values = setNames("dashed", obs_label), name = NULL) +
    scale_y_discrete(expand = c(0.01, 0)) +
    labs(title = paste0(baseline_name, "-influenced languages against their unrelated nulls"),
         x = "Cosine similarity", y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}

null_influenced_span <- plot_influenced_nulls("span_influenced", "cossim_span", "Spanish",  "#2ca6a4")
null_influenced_jap  <- plot_influenced_nulls("jap_influenced",  "cossim_jap",  "Japanese", "#8fb339")
null_influenced_eng  <- plot_influenced_nulls("eng_influenced",  "cossim_eng",  "English",  "#e07a5f")

null_influenced_span
null_influenced_jap
null_influenced_eng

ggsave(here("figures", "grammar", "distributions", "grammar_influenced_null_spanish.png"),
       null_influenced_span, width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "grammar", "distributions", "grammar_influenced_null_japanese.png"),
       null_influenced_jap,  width = 7, height = 7,   units = "in", dpi = 300)
ggsave(here("figures", "grammar", "distributions", "grammar_influenced_null_english.png"),
       null_influenced_eng,  width = 7, height = 4.5, units = "in", dpi = 300)

# --- Influenced baseline vs. the unrelated baseline (one plot per baseline) ---
# Two population distributions per plot, over the SAME set of languages that the
# baseline influenced: their observed similarity to the baseline vs. their own
# unrelated baseline. Both groups are therefore n = 10 / 20 / 10 for Spanish /
# Japanese / English. Same ridge style as the overview plot.
plot_influenced_vs_unrelated <- function(inf_col, obs_col, baseline_name, fill_color) {
  sub <- GRAMMAR_cossim |> filter(.data[[inf_col]])
  scores <- bind_rows(
    sub |> transmute(Language = "Unrelated", Similarity_Score = cossim_unr),
    sub |> transmute(Language = baseline_name, Similarity_Score = .data[[obs_col]])
  ) |>
    mutate(Language = factor(Language, levels = c("Unrelated", baseline_name)))

  # axis labels carry the sample size per group
  labels <- scores |>
    dplyr::count(Language) |>
    mutate(label = paste0(Language, " (n=", n, ")")) |>
    (\(d) setNames(d$label, d$Language))()

  score_means <- scores |>
    group_by(Language) |>
    summarize(mean_score = mean(Similarity_Score), .groups = "drop")

  ggplot(scores, aes(x = Similarity_Score, y = Language, fill = Language)) +
    geom_density_ridges(alpha = 0.5, scale = 1.1, color = "black") +
    geom_segment(
      data = score_means,
      aes(x = mean_score, xend = mean_score,
          y = as.numeric(factor(Language, levels = levels(scores$Language))),
          yend = as.numeric(factor(Language, levels = levels(scores$Language))) + 0.9,
          color = Language),
      linetype = "dashed", linewidth = 1.2, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = setNames(c("grey70", fill_color), c("Unrelated", baseline_name))) +
    scale_color_manual(values = setNames(c("grey40", fill_color), c("Unrelated", baseline_name))) +
    scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
    scale_y_discrete(expand = c(0.01, 0), labels = labels) +
    labs(title = paste0(baseline_name, "-influenced languages vs. the unrelated baseline"),
         x = "Cosine similarity", y = NULL) +
    theme_minimal() +
    theme(legend.position = "none")
}

influenced_vs_unr_span <- plot_influenced_vs_unrelated("span_influenced", "cossim_span", "Spanish",  "#2ca6a4")
influenced_vs_unr_jap  <- plot_influenced_vs_unrelated("jap_influenced",  "cossim_jap",  "Japanese", "#8fb339")
influenced_vs_unr_eng  <- plot_influenced_vs_unrelated("eng_influenced",  "cossim_eng",  "English",  "#e07a5f")

influenced_vs_unr_span
influenced_vs_unr_jap
influenced_vs_unr_eng

ggsave(here("figures", "grammar", "distributions", "grammar_influenced_vs_unrelated_spanish.png"),
       influenced_vs_unr_span, width = 7, height = 3.5, units = "in", dpi = 300)
ggsave(here("figures", "grammar", "distributions", "grammar_influenced_vs_unrelated_japanese.png"),
       influenced_vs_unr_jap,  width = 7, height = 3.5, units = "in", dpi = 300)
ggsave(here("figures", "grammar", "distributions", "grammar_influenced_vs_unrelated_english.png"),
       influenced_vs_unr_eng,  width = 7, height = 3.5, units = "in", dpi = 300)

# --- Venn diagram of the influence classification ----------------------------
# How the span/jap/eng-influenced sets overlap (a language can be influenced by
# more than one interest language). Drawn by hand (no Venn package dependency):
# three circles, each region listing the languages that fall in it.
venn_pal <- c(Spanish = "#2ca6a4", Japanese = "#8fb339", English = "#e07a5f")

# assign each influenced language to exactly one of the seven regions
region_langs <- GRAMMAR_cossim |>
  filter(span_influenced | jap_influenced | eng_influenced) |>
  mutate(region = case_when(
    span_influenced & jap_influenced & eng_influenced ~ "SJE",
    span_influenced & jap_influenced                  ~ "SJ",
    span_influenced & eng_influenced                  ~ "SE",
    jap_influenced  & eng_influenced                  ~ "JE",
    span_influenced                                   ~ "S_only",
    jap_influenced                                    ~ "J_only",
    eng_influenced                                    ~ "E_only"
  )) |>
  summarise(langs = paste(sort(language), collapse = "\n"), .by = region)

# region text anchors
region_coords <- tibble(
  region = c("S_only", "J_only", "E_only", "SJ",   "SE",  "JE",   "SJE"),
  x      = c( 0.00,    -0.95,     0.95,    -0.55,   0.55,  0.00,   0.00),
  y      = c( 0.95,    -0.55,    -0.55,     0.30,   0.30, -0.62,  -0.05)
)
venn_regions <- region_coords |> inner_join(region_langs, by = "region")

# three circle outlines centered at 90 / 210 / 330 degrees
circle_pts <- function(cx, cy, r, set, n = 200) {
  t <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(t), y = cy + r * sin(t), set = set)
}
venn_circles <- bind_rows(
  circle_pts( 0.000,  0.500, 0.9, "Spanish"),
  circle_pts(-0.433, -0.250, 0.9, "Japanese"),
  circle_pts( 0.433, -0.250, 0.9, "English")
) |>
  mutate(set = factor(set, levels = c("Spanish", "Japanese", "English")))

venn_setlabs <- tibble(
  x   = c(0.00, -1.20, 1.20),
  y   = c(1.62, -0.95, -0.95),
  lab = c(paste0("Spanish (n=",  sum(GRAMMAR_cossim$span_influenced), ")"),
          paste0("Japanese (n=", sum(GRAMMAR_cossim$jap_influenced),  ")"),
          paste0("English (n=",  sum(GRAMMAR_cossim$eng_influenced),  ")"))
)

venn_plot <- ggplot() +
  geom_polygon(data = venn_circles, aes(x, y, fill = set, group = set),
               alpha = 0.35, color = "grey30") +
  geom_text(data = venn_regions, aes(x, y, label = langs),
            size = 2.6, lineheight = 0.9) +
  geom_text(data = venn_setlabs, aes(x, y, label = lab), fontface = "bold", size = 4) +
  scale_fill_manual(values = venn_pal) +
  coord_equal(clip = "off") +
  labs(title = "Languages influenced by each interest language") +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

venn_plot

ggsave(
  filename = here("figures", "grammar", "distributions", "grammar_influence_venn.png"),
  plot = venn_plot,
  width = 6, height = 6, units = "in", dpi = 300
)

# Fold per-language significance into GRAMMAR_cossim (sig_<baseline>).
influence_sig_wide <- influence_results |>
  select(baseline, language, significant) |>
  pivot_wider(names_from = baseline, values_from = significant) |>
  rename(sig_span = Spanish,
         sig_jap  = Japanese,
         sig_eng  = English)

GRAMMAR_cossim_marked <- GRAMMAR_cossim |>
  left_join(influence_sig_wide, by = "language")


# =============================================================================
# =============================================================================
write.csv(GRAMMAR_cossim_marked, file = here("data", "GRAMMAR_cossim_marked.csv"), row.names = FALSE)
