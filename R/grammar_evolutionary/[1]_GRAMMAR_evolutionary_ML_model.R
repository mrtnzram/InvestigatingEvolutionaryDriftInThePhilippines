# =============================================================================
# [1] Grammar Evolutionary Maximum Likelihood Models
#
# Runs a two-stage Markov model selection over every analysable GRAMBANK feature
# on the pulses tree — ER vs. SYM vs. ARD, then 1-regime vs. 2-regime with the
# Philippine clade painted as its own regime — selecting by LRT at alpha = 0.05.
#
# Characters are MULTISTATE: k is read off each feature's own observed state set
# rather than assumed to be 2, and the free-rate count that sets every LRT df
# comes from the fitted index.matrix (ER = 1, SYM = k(k-1)/2, ARD = k(k-1)).
# GB065 and GB130 carry states {1,2,3} with no 0 state, so state labels are
# never treated as 0-indexed.
#
# Missing cells become ambiguous tip priors — an all-ones row in the prior
# matrix, which fitMk and fitmultiMk both consume multiplicatively — so a
# language is never dropped or imputed on account of one unanswered feature.
#
# Rates are instantaneous CTMC rates in [0, Inf), not probabilities. Grammar
# rates saturate far harder than phoneme rates do (max ~1e5 here against ~20
# there), and a saturated 2-regime fit is unmappable in [2] however well it
# fits, so an explicit segment-budget gate forces those back to 1 regime.
#
# Inputs:  data/grammar/grammar_evolutionary_tip_matrix.csv,
#          data/grammar/grammar_evolutionary_analysis_tree.rds
#
# Outputs: data/evolutionary/ml_model/GRAMMAR_evolutionary_model_selection.csv
#          data/evolutionary/ml_model/GRAMMAR_evolutionary_model_rates.csv
#          data/evolutionary/ml_model/GRAMMAR_evolutionary_dropped_features.csv
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(ape)
library(here)
library(phytools)
library(numDeriv)

source(here("R", "shared", "grambank_features.R"))

ALPHA <- 0.05

# A state seen in fewer than this many tips cannot identify its own rates.
MIN_STATE_TIPS <- 3L

# Expected mapped segments per simulated history, above which make.simmap in [2]
# becomes intractable. A fit at 2600 transitions/unit time already costs ~2.3 s
# and 54,000 segments per draw, so anything past this budget is forced to the
# 1-regime model rather than being handed to [2] to hang on.
MAX_SEGMENTS_PER_SIM <- 20000

# ── 1. Load data ─────────────────────────────────────────────

tip_matrix <- read_csv(
  here("data", "grammar", "grammar_evolutionary_tip_matrix.csv"),
  show_col_types = FALSE
)

tree_bundle     <- readRDS(here("data", "grammar", "grammar_evolutionary_analysis_tree.rds"))
analysis_tree   <- tree_bundle$analysis_tree
regime_tree     <- tree_bundle$regime_tree
philippine_tips <- tree_bundle$philippine_tips

stopifnot(
  "tip matrix should be 1:1 with the analysis tree" =
    nrow(tip_matrix) == Ntip(analysis_tree) &&
    setequal(tip_matrix$tip, analysis_tree$tip.label),
  "analysis tree must be binary" = is.binary(analysis_tree),
  "regime tree must share the analysis tree's edge matrix" =
    identical(regime_tree$edge, analysis_tree$edge)
)

# Total time available to each regime, used by the saturation gate below.
regime_branch_length <- colSums(regime_tree$mapped.edge)

message(
  "Loaded analysis tree: ", Ntip(analysis_tree), " tips (",
  length(philippine_tips), " Philippine / ",
  Ntip(analysis_tree) - length(philippine_tips), " background).\n",
  "Regime branch length: Background ", round(regime_branch_length[["Background"]], 1),
  ", Philippine ", round(regime_branch_length[["Philippine"]], 1), "."
)

# ── 2. Per-feature state inventory ───────────────────────────
# States are read from the data, not assumed, and trimmed to those attested
# often enough to identify a rate. A trimmed state's tips become ambiguous
# rather than dropped, so the topology stays identical across every feature and
# [2] can build its tree split once.

feature_values <- function(feature) {
  v <- tip_matrix[[feature]][match(analysis_tree$tip.label, tip_matrix$tip)]
  names(v) <- analysis_tree$tip.label
  v
}

state_inventory <- lapply(GRAMBANK_FEATURE_IDS, function(feature) {
  v  <- feature_values(feature)
  tb <- table(as.character(v[!is.na(v)]))
  keep <- names(tb)[tb >= MIN_STATE_TIPS]
  tibble(
    feature       = feature,
    n_observed    = sum(!is.na(v)),
    n_missing     = sum(is.na(v)),
    n_ph_observed = sum(!is.na(v[philippine_tips])),
    states_raw    = paste(names(tb), tb, sep = ":", collapse = " "),
    n_states_raw  = length(tb),
    states        = paste(keep, collapse = "/"),
    k             = length(keep)
  )
}) |> bind_rows()

analysable <- state_inventory |> filter(k >= 2) |> pull(feature)

message(
  "\n", length(analysable), " of ", length(GRAMBANK_FEATURE_IDS),
  " features analysable (>= 2 states, each attested in >= ", MIN_STATE_TIPS, " tips).\n",
  "  multistate (k > 2): ",
  paste(state_inventory$feature[state_inventory$k > 2], collapse = ", ")
)

# ── 3. Model-fitting helpers ─────────────────────────────────

# Downgrades errors to NULL and collects warnings so one bad feature cannot
# abort the loop.
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

# Tip priors over `states`: a point mass where the state is observed and kept, an
# all-ones row where it is missing or was trimmed. fitMk indexes rows by
# tip.label and reads states off colnames, so both must be set.
state_priors <- function(v, states) {
  X <- matrix(0, length(v), length(states), dimnames = list(names(v), states))
  obs <- !is.na(v) & as.character(v) %in% states
  X[!obs, ] <- 1
  X[cbind(which(obs), match(as.character(v[obs]), states))] <- 1
  X
}

# Free rates in a fitted model, from its own index.matrix: ER = 1,
# SYM = k(k-1)/2, ARD = k(k-1), without hard-coding any of them.
npar_of <- function(fit) max(fit$index.matrix, na.rm = TRUE)

# fitmultiMk stacks its per-regime rate blocks in fit$regimes order.
regime_offset <- function(fit, regime) (match(regime, fit$regimes) - 1L) * npar_of(fit)

# Q from a fitted model, looked up by transition via index.matrix rather than by
# position. Diagonals come from the row sums so no downstream consumer can be
# handed an invalid generator.
q_from_fit <- function(fit, offset, states) {
  k <- length(states)
  Q <- matrix(0, k, k, dimnames = list(states, states))
  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i != j) Q[i, j] <- unname(fit$rates[offset + fit$index.matrix[i, j]])
    }
  }
  diag(Q) <- -rowSums(Q)
  Q
}

# 2-regime standard errors come from a numeric Hessian of fitmultiMk's $lik,
# because fitMk's own $lik errors in phytools 2.5.2.
se_from_lik <- function(fit) {
  n_par <- length(fit$rates)
  h <- tryCatch(numDeriv::hessian(fit$lik, fit$rates), error = function(e) NULL)
  if (is.null(h)) return(rep(NA_real_, n_par))
  suppressWarnings(
    tryCatch(sqrt(diag(solve(-h))), error = function(e) rep(NA_real_, n_par))
  )
}

max_offdiag <- function(Q) max(Q[row(Q) != col(Q)])

# Long-format rate rows for one fitted Q; the wide er2_/ard2_ schema the phoneme
# script uses cannot hold a k x k generator.
rate_rows <- function(feature, model, regime, Q, se = NA_real_) {
  states <- rownames(Q)
  expand_grid(from_state = states, to_state = states) |>
    filter(from_state != to_state) |>
    mutate(
      feature       = feature,
      k             = length(states),
      winning_model = model,
      regime        = regime,
      rate          = Q[cbind(match(from_state, states), match(to_state, states))],
      se            = se
    ) |>
    dplyr::select(feature, k, winning_model, regime, from_state, to_state, rate, se)
}

# ── 4. Loop ──────────────────────────────────────────────────

results   <- vector("list", length(analysable))
rate_list <- vector("list", length(analysable))
t_start   <- Sys.time()

for (i in seq_along(analysable)) {
  feature <- analysable[i]
  inv     <- state_inventory |> filter(feature == !!feature)
  states  <- strsplit(inv$states, "/")[[1]]
  v       <- feature_values(feature)
  X       <- state_priors(v, states)
  notes   <- character(0)

  row <- tibble(
    feature         = feature,
    k               = length(states),
    states          = inv$states,
    stage1_winner   = NA_character_,
    winning_model   = NA_character_,
    npar_per_regime = NA_integer_,
    loglik_er       = NA_real_,
    loglik_sym      = NA_real_,
    loglik_ard      = NA_real_,
    loglik_1reg     = NA_real_,
    loglik_2reg     = NA_real_,
    p_er_vs_sym     = NA_real_,
    p_sym_vs_ard    = NA_real_,
    p_1reg_vs_2reg  = NA_real_,
    reg2_se_na      = NA,
    max_rate_bg     = NA_real_,
    max_rate_ph     = NA_real_,
    exp_transitions = NA_real_,
    saturated       = NA,
    fit_note        = NA_character_
  )

  # --- Stage 1: ER -> SYM -> ARD, promoting only on a significant LRT ---
  # SYM is identical to ER when k = 2, so it is only fitted for multistate
  # features. Each step is nested in the one before, so the chain is valid.
  f_er <- try_fit(fitMk(analysis_tree, X, model = "ER"))
  notes <- c(notes, f_er$notes)

  if (is.null(f_er$value)) {
    row$fit_note <- paste(unique(notes), collapse = " | ")
    results[[i]] <- row
    next
  }

  row$loglik_er <- as.numeric(logLik(f_er$value))
  best     <- f_er$value
  winner_1 <- "ER"

  if (length(states) > 2) {
    f_sym <- try_fit(fitMk(analysis_tree, X, model = "SYM"))
    notes <- c(notes, f_sym$notes)
    if (!is.null(f_sym$value)) {
      row$loglik_sym <- as.numeric(logLik(f_sym$value))
      row$p_er_vs_sym <- pchisq(
        2 * (row$loglik_sym - row$loglik_er),
        df = npar_of(f_sym$value) - npar_of(best), lower.tail = FALSE
      )
      if (isTRUE(row$p_er_vs_sym < ALPHA)) {
        best     <- f_sym$value
        winner_1 <- "SYM"
      }
    }
  }

  f_ard <- try_fit(fitMk(analysis_tree, X, model = "ARD"))
  notes <- c(notes, f_ard$notes)
  if (!is.null(f_ard$value)) {
    row$loglik_ard <- as.numeric(logLik(f_ard$value))
    row$p_sym_vs_ard <- pchisq(
      2 * (row$loglik_ard - as.numeric(logLik(best))),
      df = npar_of(f_ard$value) - npar_of(best), lower.tail = FALSE
    )
    if (isTRUE(row$p_sym_vs_ard < ALPHA)) {
      best     <- f_ard$value
      winner_1 <- "ARD"
    }
  }

  row$stage1_winner <- winner_1

  # --- Stage 2: 1-regime vs. 2-regime, under the stage-1 winner ---
  # Both stages use fitMk, so every log-likelihood in this script shares one
  # root prior and the LRT never crosses conventions.
  f_1reg <- try_fit(fitMk(analysis_tree, X, model = winner_1))
  f_2reg <- try_fit(fitmultiMk(regime_tree, X, model = winner_1))
  notes  <- c(notes, f_1reg$notes, f_2reg$notes)

  if (is.null(f_1reg$value) || is.null(f_2reg$value)) {
    row$winning_model <- paste0(winner_1, "1")
    row$fit_note <- paste(unique(notes), collapse = " | ")
    results[[i]] <- row
    next
  }

  npar <- npar_of(f_1reg$value)
  row$npar_per_regime <- npar
  row$loglik_1reg <- as.numeric(logLik(f_1reg$value))
  row$loglik_2reg <- f_2reg$value$logLik

  row$p_1reg_vs_2reg <- pchisq(
    2 * (row$loglik_2reg - row$loglik_1reg), df = npar, lower.tail = FALSE
  )

  se_2 <- se_from_lik(f_2reg$value)
  row$reg2_se_na <- !all(is.finite(se_2))

  Q_bg <- q_from_fit(f_2reg$value, regime_offset(f_2reg$value, "Background"), states)
  Q_ph <- q_from_fit(f_2reg$value, regime_offset(f_2reg$value, "Philippine"), states)

  row$max_rate_bg <- max_offdiag(Q_bg)
  row$max_rate_ph <- max_offdiag(Q_ph)
  row$exp_transitions <-
    row$max_rate_bg * regime_branch_length[["Background"]] +
    row$max_rate_ph * regime_branch_length[["Philippine"]]
  row$saturated <- row$exp_transitions > MAX_SEGMENTS_PER_SIM

  # The finite-SE gate alone does not catch saturation: GB071 and GB121 both
  # return finite standard errors while implying millions of segments per draw.
  winner_2 <- if (isTRUE(row$p_1reg_vs_2reg < ALPHA) &&
                  !row$reg2_se_na && !row$saturated) "2" else "1"

  if (isTRUE(row$p_1reg_vs_2reg < ALPHA)) {
    if (row$reg2_se_na) {
      notes <- c(notes, "2-regime favoured by LRT but SE non-finite; forced 1-regime")
    }
    if (row$saturated) {
      notes <- c(notes, paste0(
        "2-regime favoured by LRT but saturated (",
        signif(row$exp_transitions, 3), " segments/sim); forced 1-regime"
      ))
    }
  }

  row$winning_model <- paste0(winner_1, winner_2)

  rate_list[[i]] <- if (winner_2 == "2") {
    bind_rows(
      rate_rows(feature, row$winning_model, "Background", Q_bg),
      rate_rows(feature, row$winning_model, "Philippine", Q_ph)
    )
  } else {
    rate_rows(feature, row$winning_model, "Single", q_from_fit(f_1reg$value, 0L, states))
  }

  row$fit_note <- paste(unique(notes), collapse = " | ")
  results[[i]] <- row

  if (i %% 10 == 0 || i == length(analysable)) {
    message(
      "  fitted ", i, "/", length(analysable), " (",
      round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min)"
    )
  }
}

# ── 5. Assemble and write outputs ────────────────────────────

model_selection <- bind_rows(results) |>
  left_join(
    state_inventory |> dplyr::select(feature, n_observed, n_missing, n_ph_observed),
    by = "feature"
  ) |>
  dplyr::select(
    feature, k, states, n_observed, n_missing, n_ph_observed,
    stage1_winner, winning_model, npar_per_regime,
    loglik_er, loglik_sym, loglik_ard, loglik_1reg, loglik_2reg,
    p_er_vs_sym, p_sym_vs_ard, p_1reg_vs_2reg,
    reg2_se_na, max_rate_bg, max_rate_ph, exp_transitions, saturated, fit_note
  ) |>
  arrange(feature)

model_rates <- bind_rows(rate_list) |> arrange(feature, regime, from_state, to_state)

dropped_features <- state_inventory |>
  filter(k < 2) |>
  mutate(
    drop_reason = if_else(
      n_states_raw <= 1,
      "monomorphic",
      paste0("no second state attested in >= ", MIN_STATE_TIPS, " tips")
    )
  ) |>
  dplyr::select(
    feature, n_observed, n_missing, n_states_raw, states_raw, k, drop_reason
  ) |>
  arrange(feature)

# fitMk returns -1e+50 as a failure sentinel rather than erroring, so a finite
# log-likelihood everywhere is what proves the sibling-conflict drop in [0] did
# its job.
n_sentinel <- sum(model_selection$loglik_1reg < -1e49, na.rm = TRUE)
stopifnot(
  "fitMk returned its -1e+50 failure sentinel; check [0]'s sibling-conflict drop" =
    n_sentinel == 0L
)

dir.create(here("data", "evolutionary", "ml_model"), recursive = TRUE, showWarnings = FALSE)
write_csv(model_selection, here("data", "evolutionary", "ml_model", "GRAMMAR_evolutionary_model_selection.csv"))
write_csv(model_rates, here("data", "evolutionary", "ml_model", "GRAMMAR_evolutionary_model_rates.csv"))
write_csv(dropped_features, here("data", "evolutionary", "ml_model", "GRAMMAR_evolutionary_dropped_features.csv"))

message(
  "\nDone in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min.",
  "\n  model selection : ", nrow(model_selection), " features -> ",
  "data/evolutionary/ml_model/GRAMMAR_evolutionary_model_selection.csv",
  "\n  rates (long)    : ", nrow(model_rates), " rows -> ",
  "data/evolutionary/ml_model/GRAMMAR_evolutionary_model_rates.csv",
  "\n  dropped         : ", nrow(dropped_features), " features -> ",
  "data/evolutionary/ml_model/GRAMMAR_evolutionary_dropped_features.csv\n"
)

print(table(model_selection$winning_model, useNA = "ifany"))

message(
  "\nSaturated fits forced back to 1 regime: ",
  paste(model_selection$feature[which(model_selection$saturated)], collapse = ", ")
)
message(
  "2-regime winners: ",
  paste(model_selection$feature[grepl("2$", model_selection$winning_model)], collapse = ", ")
)
