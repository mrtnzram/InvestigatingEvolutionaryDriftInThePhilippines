# =============================================================================
# [1] Phoneme Evolutionary Maximum Likelihood Models
#
# Runs a two-stage Markov model selection over every analysable Austronesian
# phoneme on the pulses tree — ER vs. ARD, then 1-regime vs. 2-regime with the
# Philippine clade painted as its own regime — selecting by LRT at alpha = 0.05
# and falling back to the simpler model when the complex fit returns a
# non-finite standard error. Phonemes present in <=2 tips or in all of them are
# unidentifiable and go to the dropped table.
#
# Rates are instantaneous CTMC rates in [0, Inf), not probabilities; saturated
# fits have an arbitrary magnitude but a meaningful q01:q10 ratio.
#
# Inputs:  data/RUHLENdf_PH_evolutionary.csv, data/phoneme_evolutionary_tip_mapping.csv,
#          data/phoneme_evolutionary_analysis_tree.rds, data/phoneme_key_pnas.csv
#
# Outputs: data/PHONEME_evolutionary_model_selection.csv
#          data/PHONEME_evolutionary_dropped_phonemes.csv
# =============================================================================

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(phytools)
library(numDeriv)

ALPHA <- 0.05

# ── 1. Load data ─────────────────────────────────────────────

RUHLENdf <- read_csv(
  here("data", "phoneme", "evolutionary", "initial_datasets", "RUHLENdf_PH_evolutionary.csv"),
  show_col_types = FALSE
)

tip_mapping <- read_csv(
  here("data", "phoneme", "evolutionary", "initial_datasets", "phoneme_evolutionary_tip_mapping.csv"),
  show_col_types = FALSE
)

# Tree pruning, conflict-tip resolution, and the Philippine regime are decided
# once in [0]_Phylogenetic_Tree.R and loaded here, not rebuilt.
tree_bundle     <- readRDS(here("data", "phoneme", "evolutionary", "initial_datasets", "phoneme_evolutionary_analysis_tree.rds"))
analysis_tree   <- tree_bundle$analysis_tree
regime_tree     <- tree_bundle$regime_tree
philippine_tips <- tree_bundle$philippine_tips

plotSimmap(
  regime_tree,
  type  = "fan",
  fsize = 0.5,
  lwd = 2
)

phoneme_cols <- grep("^phoneme_", names(RUHLENdf), value = TRUE)

# phoneme_n -> IPA keymap, built and exported by [0]_CREANZA_RUHLENdatabase.R
pnas_key <- read_csv(
  here("data", "phoneme", "evolutionary", "initial_datasets", "phoneme_key_pnas.csv"),
  show_col_types = FALSE
) |>
  select(phoneme = phoneme_id, ipa, class)

# ── 2. Tip-level table ───────────────────────────────────────
# Restricted to the loaded analysis tree's own tip set, so conflict-tip
# resolution never needs to be re-applied here.

tips <- tip_mapping |>
  filter(tip %in% analysis_tree$tip.label) |>
  left_join(
    RUHLENdf,
    by = c("iso6393", "ruhlen_language" = "language")
  )

stopifnot(
  "tip_mapping should join 1:1 onto the analysis tree" =
    nrow(tips) == Ntip(analysis_tree),
  "every tip must resolve a Language_type" =
    !any(is.na(tips$Language_type)),
  "phoneme matrix must be complete for all tips" =
    !any(is.na(tips[phoneme_cols])),
  "analysis tree must be binary" =
    is.binary(analysis_tree),
  "[0]'s Philippine tips must match tip-level Language_type" =
    setequal(philippine_tips, tips$tip[tips$Language_type == "Philippine Language"])
)

message(
  "Loaded analysis tree: ", Ntip(analysis_tree), " tips (",
  length(philippine_tips), " Philippine / ",
  Ntip(analysis_tree) - length(philippine_tips), " background)."
)

# ── 3. Prevalence table (all 728 phonemes) ───────────────────
# Denominators are the tree-matched languages, so prevalences describe the
# languages the models actually see; unrelated controls use the full Ruhlen set.

prevalence_of <- function(data, prefix) {
  n <- nrow(data)
  data |>
    summarise(across(all_of(phoneme_cols), \(x) sum(x))) |>
    pivot_longer(everything(), names_to = "phoneme", values_to = "n") |>
    transmute(
      phoneme,
      "n_{prefix}"          := n,
      "{prefix}_prevalence" := n / .env$n
    )
}

unrelated_langs <- RUHLENdf |> filter(Language_type == "Unrelated Language")
interest_langs  <- RUHLENdf |> filter(Language_type == "Language of Interest")

has_phoneme <- function(lang) {
  row <- interest_langs |> filter(language == lang)
  stopifnot("language of interest not found" = nrow(row) == 1L)
  as.integer(row[phoneme_cols])
}

phoneme_summary <- prevalence_of(tips, "austronesian") |>
  left_join(
    prevalence_of(filter(tips, Language_type == "Philippine Language"), "ph"),
    by = "phoneme"
  ) |>
  left_join(
    prevalence_of(filter(tips, Language_type == "Austronesian"), "austronesian_nonph"),
    by = "phoneme"
  ) |>
  left_join(prevalence_of(unrelated_langs, "unrelated"), by = "phoneme") |>
  mutate(
    span_has = has_phoneme("Spanish"),
    eng_has  = has_phoneme("English"),
    jap_has  = has_phoneme("Japanese")
  ) |>
  left_join(pnas_key, by = "phoneme")

n_tree <- nrow(tips)

analysable <- phoneme_summary |>
  filter(n_austronesian > 2, n_austronesian < n_tree) |>
  pull(phoneme)

message(
  length(analysable), " phonemes analysable (present in >2 and <",
  n_tree, " tree tips)."
)

# ── 4. Model-fitting helpers ─────────────────────────────────

# Downgrades errors to NULL and collects warnings so one bad phoneme cannot abort
# the loop.
try_fit <- function(expr) {
  notes <- character(0)
  value <- withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) {
        notes <<- c(notes, paste0("ERROR: ", conditionMessage(e)))
        NULL
      }
    ),
    warning = function(w) {
      notes <<- c(notes, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, notes = notes)
}

# Look rates up by transition via index.matrix, never by position: ape orders a
# 2-state ARD as rates[1] = loss, rates[2] = gain.
rate_01 <- function(fit, offset = 0L) unname(fit$rates[offset + fit$index.matrix[1, 2]])
rate_10 <- function(fit, offset = 0L) unname(fit$rates[offset + fit$index.matrix[2, 1]])

# fitmultiMk stacks its per-regime rate blocks in fit$regimes order.
regime_offset <- function(fit, regime, k) {
  (match(regime, fit$regimes) - 1L) * k
}

# 2-regime standard errors come from a numeric Hessian of fitmultiMk's $lik,
# because fitMk's own $lik errors in phytools 2.5.2.
se_from_lik <- function(fit) {
  n_par <- length(fit$rates)
  h <- tryCatch(
    numDeriv::hessian(fit$lik, fit$rates),
    error = function(e) NULL
  )
  if (is.null(h)) return(rep(NA_real_, n_par))
  suppressWarnings(
    tryCatch(
      sqrt(diag(solve(-h))),
      error = function(e) rep(NA_real_, n_par)
    )
  )
}

state_vector <- function(phoneme_col) {
  v <- tips[[phoneme_col]][match(analysis_tree$tip.label, tips$tip)]
  v <- factor(v, levels = c(0, 1))
  names(v) <- analysis_tree$tip.label
  v
}

# ── 5. Loop ──────────────────────────────────────────────────

results <- vector("list", length(analysable))
t_start <- Sys.time()

for (i in seq_along(analysable)) {
  phoneme_col <- analysable[i]
  v <- state_vector(phoneme_col)
  notes <- character(0)

  row <- tibble(
    phoneme         = phoneme_col,
    winning_model   = NA_character_,
    er_rate         = NA_real_,
    ard_rate_01     = NA_real_,
    ard_rate_10     = NA_real_,
    er2_rate_bg     = NA_real_,
    er2_rate_ph     = NA_real_,
    ard2_rate_bg_01 = NA_real_,
    ard2_rate_bg_10 = NA_real_,
    ard2_rate_ph_01 = NA_real_,
    ard2_rate_ph_10 = NA_real_,
    p_er_vs_ard     = NA_real_,
    p_1reg_vs_2reg  = NA_real_,
    loglik_er       = NA_real_,
    loglik_ard      = NA_real_,
    loglik_1reg     = NA_real_,
    loglik_2reg     = NA_real_,
    ard_se_na       = NA,
    reg2_se_na      = NA
  )

  # --- Stage 1: ER vs. ARD ---
  f_er  <- try_fit(ace(v, analysis_tree, type = "discrete", model = "ER"))
  f_ard <- try_fit(ace(v, analysis_tree, type = "discrete", model = "ARD"))
  notes <- c(notes, f_er$notes, f_ard$notes)

  if (is.null(f_er$value) || is.null(f_ard$value)) {
    row$fit_note <- paste(unique(notes), collapse = " | ")
    results[[i]] <- row
    next
  }

  row$er_rate     <- unname(f_er$value$rates[1])
  row$ard_rate_01 <- rate_01(f_ard$value)
  row$ard_rate_10 <- rate_10(f_ard$value)
  row$loglik_er   <- f_er$value$loglik
  row$loglik_ard  <- f_ard$value$loglik

  cmp <- try_fit(anova(f_er$value, f_ard$value))
  notes <- c(notes, cmp$notes)
  row$p_er_vs_ard <- if (is.null(cmp$value)) NA_real_ else cmp$value[2, "Pr(>|Chi|)"]

  row$ard_se_na <- !all(is.finite(f_ard$value$se))
  winner_1 <- if (isTRUE(row$p_er_vs_ard < ALPHA) && !row$ard_se_na) "ARD" else "ER"
  if (row$ard_se_na && isTRUE(row$p_er_vs_ard < ALPHA)) {
    notes <- c(notes, "ARD favoured by LRT but SE non-finite; forced ER")
  }

  k <- if (winner_1 == "ARD") 2L else 1L

  # --- Stage 2: 1-regime vs. 2-regime, under the stage-1 winner ---
  # fitMk and fitmultiMk share a root prior; ace does not, so the LRT must
  # compare these two and never reuse the ace log-likelihood above.
  f_1reg <- try_fit(fitMk(analysis_tree, v, model = winner_1))
  f_2reg <- try_fit(fitmultiMk(regime_tree, v, model = winner_1))
  notes  <- c(notes, f_1reg$notes, f_2reg$notes)

  if (is.null(f_1reg$value) || is.null(f_2reg$value)) {
    row$winning_model <- paste0(winner_1, "1")
    row$fit_note <- paste(unique(notes), collapse = " | ")
    results[[i]] <- row
    next
  }

  row$loglik_1reg <- as.numeric(logLik(f_1reg$value))
  row$loglik_2reg <- f_2reg$value$logLik

  lr <- 2 * (row$loglik_2reg - row$loglik_1reg)
  row$p_1reg_vs_2reg <- pchisq(lr, df = k, lower.tail = FALSE)

  se_2 <- se_from_lik(f_2reg$value)
  row$reg2_se_na <- !all(is.finite(se_2))

  winner_2 <- if (isTRUE(row$p_1reg_vs_2reg < ALPHA) && !row$reg2_se_na) "2" else "1"
  if (row$reg2_se_na && isTRUE(row$p_1reg_vs_2reg < ALPHA)) {
    notes <- c(notes, "2-regime favoured by LRT but SE non-finite; forced 1-regime")
  }

  off_bg <- regime_offset(f_2reg$value, "Background", k)
  off_ph <- regime_offset(f_2reg$value, "Philippine", k)

  if (winner_1 == "ER") {
    row$er2_rate_bg <- rate_01(f_2reg$value, off_bg)
    row$er2_rate_ph <- rate_01(f_2reg$value, off_ph)
  } else {
    row$ard2_rate_bg_01 <- rate_01(f_2reg$value, off_bg)
    row$ard2_rate_bg_10 <- rate_10(f_2reg$value, off_bg)
    row$ard2_rate_ph_01 <- rate_01(f_2reg$value, off_ph)
    row$ard2_rate_ph_10 <- rate_10(f_2reg$value, off_ph)
  }

  row$winning_model <- paste0(winner_1, winner_2)
  row$fit_note <- paste(unique(notes), collapse = " | ")
  results[[i]] <- row

  if (i %% 10 == 0 || i == length(analysable)) {
    message(
      "  fitted ", i, "/", length(analysable), " (",
      round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min)"
    )
  }
}

# ── 6. Assemble and write outputs ────────────────────────────

model_selection <- bind_rows(results) |>
  left_join(phoneme_summary, by = "phoneme") |>
  select(
    phoneme, ipa, class, winning_model,
    er_rate, ard_rate_01, ard_rate_10,
    er2_rate_bg, er2_rate_ph,
    ard2_rate_bg_01, ard2_rate_bg_10, ard2_rate_ph_01, ard2_rate_ph_10,
    span_has, eng_has, jap_has,
    ph_prevalence, unrelated_prevalence, austronesian_prevalence,
    n_austronesian, n_ph, austronesian_nonph_prevalence,
    p_er_vs_ard, p_1reg_vs_2reg,
    loglik_er, loglik_ard, loglik_1reg, loglik_2reg,
    ard_se_na, reg2_se_na, fit_note
  ) |>
  arrange(phoneme)

# Excluded phonemes that are still attested on the tree, i.e. real but
# unanalysable rather than simply absent.
dropped_phonemes <- phoneme_summary |>
  filter(n_austronesian > 0, !phoneme %in% analysable) |>
  mutate(
    drop_reason = if_else(
      n_austronesian == n_tree, "100% prevalence", "<=2 languages"
    )
  ) |>
  select(
    phoneme, ipa, class,
    span_has, eng_has, jap_has,
    ph_prevalence, unrelated_prevalence, austronesian_prevalence,
    n_austronesian, drop_reason
  ) |>
  arrange(phoneme)

write_csv(model_selection, here("data", "phoneme", "evolutionary", "ML_model", "PHONEME_evolutionary_model_selection.csv"))
write_csv(dropped_phonemes, here("data", "phoneme", "evolutionary", "ML_model", "PHONEME_evolutionary_dropped_phonemes.csv"))

message(
  "\nDone in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min.",
  "\n  model selection : ", nrow(model_selection), " phonemes -> ",
  "data/PHONEME_evolutionary_model_selection.csv",
  "\n  dropped         : ", nrow(dropped_phonemes), " phonemes -> ",
  "data/PHONEME_evolutionary_dropped_phonemes.csv\n"
)

print(table(model_selection$winning_model, useNA = "ifany"))
