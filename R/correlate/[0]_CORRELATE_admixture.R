# =============================================================================
# [0] Correlate Analysis — Linguistic influence vs. Spanish genetic admixture
#
# Tests whether Spanish-colonial linguistic influence — loanword borrowing
# (cognate), phoneme-inventory similarity to Spanish (phoneme), or grammatical
# similarity to Spanish (grammar) — predicts the presence of Spanish genetic
# admixture in the population(s) speaking that language. One logistic
# regression per linguistic domain (admixture_exists ~ <domain predictor>),
# fit with glm(family = binomial), sharing one genetic-admixture prep.
#
# Predictors are continuous: number_of_loans for cognate, and for phoneme and
# grammar span_delta = cossim_span - cossim_unr — the same quantity [2]'s
# Gaussian mixture thresholded into span_influenced, used here undichotomized
# so the ordering is not discarded.
#
# NOTE: Chavacano is absent from all three models and cannot be added. It is the
# strongest genetic signal in the dataset (span_admx 0.053 at 4.4x its SE, MALDER
# date inside the colonial window) and, as a Spanish-lexified creole, the natural
# high anchor for every predictor here — but Glottolog files it as Indo-European
# (chav1241, Latino-Faliscan) and it has no PHOIBLE inventory, no Grambank
# features, and no ABVD wordlist. The regressions therefore span a range that
# excludes the clearest case of the phenomenon.
#
# admixture_exists is kept as the response in all three models: for a single
# binary predictor and binary outcome the fitted log-odds-ratio is the same
# either direction (equivalent to the 2x2 table), but the intercept — and so
# the predicted-probability curve — is anchored to whichever variable is the
# response, and the ask here is a probability score *for admixture*.
#
# The response is coded span_admx > span_se ("detectably nonzero"), not
# span_admx > 0 — see ADMX_DETECTABLE below. Event counts are small under this
# rule (4-5 per domain), so every coefficient here is imprecise; treat the
# direction as suggestive and the magnitudes as unresolved.
#
# Only cognate carries glottocode already; phoneme and grammar are bridged to
# it per domain: grammar's GRAMBANKdf_full$Language_ID *is* the glottocode;
# phoneme's RUHLENdf_PH$iso6393 resolves to glottocode via
# lingtypology::gltc.iso() (the reverse of iso.gltc(), already used in
# [0]_CREANZA_RUHLENdatabase.R).
#
# Inputs:  data/network_distance/GENETIC_final.csv
#          data/genetic/GENETIC_subgroup_lookup.csv
#          data/network_distance/COGNATE_final.csv
#          data/cosine_distribution/GRAMMAR_cossim_marked.csv, data/grammar/GRAMBANKdf_full.csv
#          data/cosine_distribution/PHONEME_cossim_marked.csv, data/phoneme/RUHLENdf_PH.csv
# Outputs: data/correlate/CORRELATE_<domain>_admixture.csv      (merged analysis tables)
#          data/correlate/CORRELATE_<domain>_admixture_glm.rds  (fitted models)
#          data/correlate/CORRELATE_admixture_prob_by_domain.csv (combined probabilities)
#          figures/correlate/<domain>_admixture_glm.png
# =============================================================================

library(here)
library(tidyverse)
library(ggplot2)
library(lingtypology)

dir.create(here("data", "correlate"), showWarnings = FALSE)
dir.create(here("figures", "correlate"), showWarnings = FALSE)

# SPAN_SE_MAX: a population whose span_se exceeds the largest observed span_admx
# (0.083) cannot test positive at any true value — it is unmeasurable, not
# negative, and coding it 0 is a misclassification. One absolute bar across all
# domains, replacing an earlier per-domain quantile that held cognate, grammar
# and phoneme to three different standards (0.0935 / 0.017 / 0.168).
SPAN_SE_MAX <- 0.083

# ADMX_DETECTABLE: require the qpAdm estimate to exceed its own standard error
# before calling admixture present, rather than merely being > 0. Of the 42
# populations with span_admx > 0, only 9 clear their own SE and 2 clear 2*SE —
# negatives are already floored to 0 upstream in [0]_GENETIC_ADMX_MALDER.R, so
# a "> 0" rule largely asks which side of zero a noisy point estimate landed on,
# and it reverses the sign of the fitted coefficient in all three domains.
# Set FALSE to reproduce the original "> 0" coding.
ADMX_DETECTABLE <- TRUE


# ── 1. Genetic prep (shared across all three domains) ───────────────────────
GENETIC_final  <- read_csv(here("data", "network_distance", "GENETIC_final.csv"),
                            show_col_types = FALSE)
GENETIC_lookup <- read_csv(here("data", "genetic", "GENETIC_subgroup_lookup.csv"),
                            show_col_types = FALSE)

GENETIC_prepped <- GENETIC_final |>
  left_join(GENETIC_lookup |> select(population, glottocode), by = "population")

n_na_admx <- sum(is.na(GENETIC_prepped$span_admx))
GENETIC_prepped <- GENETIC_prepped |>
  filter(!is.na(span_admx)) |>
  mutate(admixture_exists = if (ADMX_DETECTABLE) {
           as.integer(span_admx > span_se)
         } else {
           as.integer(span_admx > 0)
         })

message(n_na_admx, " populations dropped (span_admx is NA — qpAdm source populations) | ",
        nrow(GENETIC_prepped), " retained.")
message("admixture_exists rule: span_admx > ",
        if (ADMX_DETECTABLE) "span_se" else "0", " | ",
        sum(GENETIC_prepped$admixture_exists), " of ", nrow(GENETIC_prepped),
        " populations positive.")

# Under the detectable rule a population whose span_se exceeds the largest
# observed span_admx can never test positive — it is unmeasurable, not negative,
# and coding it 0 is a misclassification the SE filter below only partly offsets.
if (ADMX_DETECTABLE) {
  n_censored <- sum(GENETIC_prepped$span_se > max(GENETIC_prepped$span_admx))
  message(n_censored, " populations are structurally unable to test positive ",
          "(span_se > max span_admx = ",
          round(max(GENETIC_prepped$span_admx), 3), ") — counted as negative.")
}


# ── 2. Per-domain join + SE filter ───────────────────────────────────────────
# Shared step once a domain's frame carries `glottocode` and its continuous
# predictor column: join to genetic, drop the top SE quartile, report counts.
# many-to-many: legitimate on both sides — a glottocode can have >1 genetic
# population, and (for phoneme/grammar) >1 dialect-level language name.
join_and_filter <- function(pred_df, predictor, domain_label) {
  merged <- pred_df |>
    inner_join(
      GENETIC_prepped |> select(population, glottocode, span_admx, span_se, admixture_exists),
      by = "glottocode",
      relationship = "many-to-many"
    )

  n_before <- nrow(merged)
  merged   <- merged |> filter(span_se < SPAN_SE_MAX)

  message(domain_label, ": ", n_before - nrow(merged), " populations dropped (span_se >= ",
          SPAN_SE_MAX, ", unmeasurable) | ",
          nrow(merged), " retained, ", sum(merged$admixture_exists), " positive.")

  stopifnot(
    "predictor column has NA after the join"  = !anyNA(merged[[predictor]]),
    "admixture_exists has NA after the join"   = !anyNA(merged$admixture_exists),
    "predictor column has no variation"        = length(unique(merged[[predictor]])) > 1,
    "admixture_exists is constant"             = length(unique(merged$admixture_exists)) == 2
  )

  merged
}


# ── 3. Cognate: loan count ───────────────────────────────────────────────────
# number_of_loans, not the >0 indicator: binarizing collapses 1 loanword and 13
# into the same value, discarding the ordering the small sample can least afford
# to lose.
COGNATE_final <- read_csv(here("data", "network_distance", "COGNATE_final.csv"),
                           show_col_types = FALSE)

COGNATE_pred <- COGNATE_final |>
  select(glottocode, language, number_of_loans)

cognate_merged <- join_and_filter(COGNATE_pred, "number_of_loans", "Cognate")

write.csv(cognate_merged, here("data", "correlate", "CORRELATE_cognate_admixture.csv"),
          row.names = FALSE)


# ── 4. Grammar: Spanish-similarity delta ─────────────────────────────────────
# span_delta = cossim_span - cossim_unr, the same quantity [2]'s Gaussian
# mixture thresholded into span_influenced — used here undichotomized.
# GRAMMAR_cossim_marked carries no glottocode; GRAMBANKdf_full$Language_ID
# already *is* the glottocode (renamed from it in [0]_GRAMBANKdatabase.R).
GRAMMAR_cossim  <- read_csv(here("data", "cosine_distribution", "GRAMMAR_cossim_marked.csv"),
                             show_col_types = FALSE)
GRAMBANKdf_full <- read_csv(here("data", "grammar", "GRAMBANKdf_full.csv"),
                             show_col_types = FALSE)

GRAMMAR_pred <- GRAMMAR_cossim |>
  left_join(
    GRAMBANKdf_full |> select(language, glottocode = Language_ID) |>
      distinct(language, .keep_all = TRUE),
    by = "language"
  ) |>
  mutate(span_delta = cossim_span - cossim_unr) |>
  select(glottocode, language, span_delta)

grammar_merged <- join_and_filter(GRAMMAR_pred, "span_delta", "Grammar")

write.csv(grammar_merged, here("data", "correlate", "CORRELATE_grammar_admixture.csv"),
          row.names = FALSE)


# ── 5. Phoneme: Spanish-similarity delta ─────────────────────────────────────
# Same span_delta construction as grammar. PHONEME_cossim_marked carries no
# glottocode either; bridge via RUHLENdf_PH's iso6393, resolved to glottocode
# with lingtypology::gltc.iso().
PHONEME_cossim <- read_csv(here("data", "cosine_distribution", "PHONEME_cossim_marked.csv"),
                            show_col_types = FALSE)
RUHLENdf_PH    <- read_csv(here("data", "phoneme", "RUHLENdf_PH.csv"),
                            show_col_types = FALSE)

PHONEME_pred <- PHONEME_cossim |>
  left_join(
    RUHLENdf_PH |> select(language, iso6393) |> distinct(language, .keep_all = TRUE),
    by = "language"
  ) |>
  mutate(glottocode = gltc.iso(iso6393),
         span_delta = cossim_span - cossim_unr) |>
  select(glottocode, language, span_delta)

phoneme_merged <- join_and_filter(PHONEME_pred, "span_delta", "Phoneme")

write.csv(phoneme_merged, here("data", "correlate", "CORRELATE_phoneme_admixture.csv"),
          row.names = FALSE)


# ── 6. Models: admixture_exists ~ <domain predictor> ─────────────────────────
m_admx_cognate <- glm(admixture_exists ~ number_of_loans,
                       data = cognate_merged, family = binomial)
m_admx_grammar <- glm(admixture_exists ~ span_delta,
                       data = grammar_merged, family = binomial)
m_admx_phoneme <- glm(admixture_exists ~ span_delta,
                       data = phoneme_merged, family = binomial)

print(summary(m_admx_cognate))
print(summary(m_admx_grammar))
print(summary(m_admx_phoneme))

saveRDS(m_admx_cognate, here("data", "correlate", "CORRELATE_cognate_admixture_glm.rds"))
saveRDS(m_admx_grammar, here("data", "correlate", "CORRELATE_grammar_admixture_glm.rds"))
saveRDS(m_admx_phoneme, here("data", "correlate", "CORRELATE_phoneme_admixture_glm.rds"))


# ── 7. Probability score + figure per domain ─────────────────────────────────
# p(x) = plogis(b0 + b1*x) via predict(..., type = "response"), which applies
# the inverse link for us. With a continuous predictor the fitted curve is drawn
# over the observed range of x, so — unlike the earlier binary coding, where the
# curve interpolated between two points — every part of it is data-supported.
# Only y is binary now, so jitter is applied on y alone; x positions stay
# truthful, apart from a small nudge for cognate's tied integer counts.
#
# Subtitle anchors the effect at one SD of the predictor: the relative change
# (p_hi - p_lo) / p_lo between x = mean(x) and x = mean(x) + sd(x), plus the
# Wald p-value for the coefficient (H0: b_pred = 0), formatted with the same
# "< 0.001" / "= x.xxx" idiom used for the Wilcoxon boxplot's sig_labels in
# [2]_*_cosine_distribution_analysis.R.
plot_admixture_glm <- function(df, model, predictor, predictor_label, domain_label,
                                jitter_w = 0) {
  x   <- df[[predictor]]
  mu  <- mean(x)
  sdx <- sd(x)

  p_lo <- predict(model, newdata = tibble(!!predictor := mu),       type = "response")
  p_hi <- predict(model, newdata = tibble(!!predictor := mu + sdx), type = "response")

  rel_pct   <- (p_hi - p_lo) / p_lo * 100
  direction <- if (rel_pct >= 0) "more" else "less"

  p_value <- summary(model)$coefficients[predictor, "Pr(>|z|)"]
  p_label <- ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))

  pad      <- 0.04 * diff(range(x))
  curve_df <- tibble(x = seq(min(x) - pad, max(x) + pad, length.out = 200)) |>
    mutate(p = predict(model, newdata = tibble(!!predictor := x), type = "response"))

  ggplot(df, aes(x = .data[[predictor]], y = admixture_exists)) +
    geom_line(data = curve_df, aes(x = x, y = p), inherit.aes = FALSE,
              color = "#2ca6a4", linewidth = 1) +
    geom_jitter(width = jitter_w, height = 0.04, alpha = 0.6) +
    scale_y_continuous(breaks = c(0, 1), labels = c("No admixture", "Admixture")) +
    labs(title = paste0(domain_label, " influence vs. Spanish admixture"),
         subtitle = sprintf("Admixture is %.0f%% %s likely per SD of %s (%s)  |  N = %d",
                             abs(rel_pct), direction, predictor_label, p_label, nrow(df)),
         x = predictor_label, y = "Spanish admixture") +
    theme_minimal()
}

# cognate's predictor is an integer count with many ties, so it gets a small
# horizontal nudge; the cosine deltas are effectively continuous and get none.
p_cognate <- plot_admixture_glm(cognate_merged, m_admx_cognate, "number_of_loans",
                                 "Spanish loanwords", "Cognate", jitter_w = 0.12)
p_grammar <- plot_admixture_glm(grammar_merged, m_admx_grammar, "span_delta",
                                 "grammar similarity to Spanish", "Grammar")
p_phoneme <- plot_admixture_glm(phoneme_merged, m_admx_phoneme, "span_delta",
                                 "phoneme similarity to Spanish", "Phoneme")

p_cognate
p_grammar
p_phoneme

ggsave(here("figures", "correlate", "cognate_admixture_glm.png"), p_cognate,
       width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "correlate", "grammar_admixture_glm.png"), p_grammar,
       width = 7, height = 4.5, units = "in", dpi = 300)
ggsave(here("figures", "correlate", "phoneme_admixture_glm.png"), p_phoneme,
       width = 7, height = 4.5, units = "in", dpi = 300)


# ── 8. Combined probability table across domains ─────────────────────────────
# With continuous predictors on different scales (loan counts vs cosine deltas),
# the raw beta is not comparable across domains — or_per_sd is, so both are
# reported alongside the fitted probability at the mean and one SD above it.
predict_prob <- function(model, predictor, df, domain_label) {
  x   <- df[[predictor]]
  mu  <- mean(x)
  sdx <- sd(x)
  cf  <- summary(model)$coefficients[predictor, ]

  tibble(
    domain            = domain_label,
    predictor         = predictor,
    beta              = cf[["Estimate"]],
    se                = cf[["Std. Error"]],
    p_value           = cf[["Pr(>|z|)"]],
    or_per_sd         = exp(cf[["Estimate"]] * sdx),
    p_at_mean         = predict(model, newdata = tibble(!!predictor := mu),       type = "response"),
    p_at_mean_plus_sd = predict(model, newdata = tibble(!!predictor := mu + sdx), type = "response"),
    N                 = nrow(df),
    n_positive        = sum(df$admixture_exists)
  )
}

CORRELATE_prob_by_domain <- bind_rows(
  predict_prob(m_admx_cognate, "number_of_loans", cognate_merged, "Cognate"),
  predict_prob(m_admx_grammar, "span_delta",      grammar_merged, "Grammar"),
  predict_prob(m_admx_phoneme, "span_delta",      phoneme_merged, "Phoneme")
)

print(CORRELATE_prob_by_domain)

write.csv(CORRELATE_prob_by_domain,
          here("data", "correlate", "CORRELATE_admixture_prob_by_domain.csv"),
          row.names = FALSE)
