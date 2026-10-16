# =============================================================================
# [4] All domains — scrub-step sensitivity
#
# Refits the [4] scrub with each geography term alone (signed geographic axis;
# distance from Manila) against both together, and compares the residual pure
# phylogeny and its Spanish-signal regression across the three versions.
#
# Inputs:  R/shared/[4]_ALL_PVR.R's per-domain script lists (sourced via run_domain())
# Outputs: data/pvr/pvr_scrub_sensitivity.csv
# =============================================================================

library(tidyverse)
library(here)
source(here("R", "shared", "consolidate_domains.R"))

SCRIPTS <- list(
  phoneme = list(part_a(here("R", "phoneme_analysis", "[0]_CREANZA_RUHLENdatabase.R")),
                 here("R", "phoneme_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "phoneme_analysis", "[4]_PHONEME_PVR.R")),
  grammar = list(part_a(here("R", "grammar_analysis", "[0]_GRAMBANKdatabase.R")),
                 here("R", "grammar_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "grammar_analysis", "[4]_GRAMMAR_PVR.R")),
  cognate = list(here("R", "cognate_analysis", "[0]_Phylogenetic_Tree.R"),
                 here("R", "cognate_analysis", "[4]_COGNATE_PVR.R")),
  genetic = list(here("R", "genetic_analysis", "[4]_GENETIC_PVR.R"))
)

VARIANTS <- list(both   = axis1 ~ s + h,
                 signed = axis1 ~ s,
                 dist   = axis1 ~ h)

term_or_na <- function(cf, term, col) if (term %in% rownames(cf)) cf[term, col] else NA_real_

# scrub_variants(): one row per scrub variant for one [4] analysis frame.
scrub_variants <- function(frame, domain, zeros = NA_character_) {
  d <- data.frame(axis1 = frame$phylo_axis1, s = frame$geo_signed_km / 100,
                  h = frame$x_km / 100, y = frame$y)
  r_both <- resid(lm(VARIANTS$both, data = d))

  imap_dfr(VARIANTS, function(f, variant) {
    m  <- lm(f, data = d)
    cf <- summary(m)$coefficients
    d$pure <- resid(m)
    m_y    <- summary(lm(y ~ pure, data = d))
    tibble(
      domain, zeros, variant, n = nrow(d),
      cor_s_h         = cor(d$s, d$h),
      vif             = 1 / (1 - cor(d$s, d$h)^2),
      scrub_r2        = summary(m)$r.squared,
      slope_s         = term_or_na(cf, "s", "Estimate"),
      se_s            = term_or_na(cf, "s", "Std. Error"),
      slope_h         = term_or_na(cf, "h", "Estimate"),
      se_h            = term_or_na(cf, "h", "Std. Error"),
      cor_pure_s      = cor(d$pure, d$s),
      cor_pure_h      = cor(d$pure, d$h),
      cor_pure_both   = cor(d$pure, r_both),
      y_pure_slope    = m_y$coefficients["pure", "Estimate"],
      y_pure_r2       = m_y$r.squared,
      y_pure_p        = m_y$coefficients["pure", "Pr(>|t|)"],
      # > 0 means y ~ h and y ~ pure overlap, so their R^2s double-count.
      additivity_gap  = summary(lm(y ~ h, data = d))$r.squared + m_y$r.squared -
                        summary(lm(y ~ h + pure, data = d))$r.squared
    )
  })
}

sensitivity <- imap_dfr(SCRIPTS, function(scripts, domain) {
  message("\n[", domain, "]")
  env <- run_domain(scripts)
  if (exists("kept", envir = env, inherits = FALSE)) {
    bind_rows(scrub_variants(env$removed$frame, domain, "removed"),
              scrub_variants(env$kept$frame,    domain, "kept"))
  } else {
    scrub_variants(env$fit$frame, domain)
  }
})

dir.create(here("data", "pvr"), recursive = TRUE, showWarnings = FALSE)
write.csv(sensitivity, here("data", "pvr", "pvr_scrub_sensitivity.csv"), row.names = FALSE)
print(sensitivity, n = Inf, width = Inf)
