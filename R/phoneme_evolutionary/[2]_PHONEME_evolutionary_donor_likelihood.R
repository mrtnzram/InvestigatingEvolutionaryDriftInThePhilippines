# =============================================================================
# [2] Phoneme Evolutionary — Donor rate-multiplier models
#
# [1] asks WHICH phonemes change at a different rate inside the Philippine clade.
# It does not test a mechanism: its 2-regime model just grants the clade a free
# rate, with no reference to what Spanish, English or Japanese actually contain.
#
# [2] supplies the mechanism. For each donor a single global Beta scales the
# phoneme-GAIN rate inside the Philippine clade, but only for phonemes that donor
# has:
#
#     q01_philippine = q01_base * exp(Beta * PH_clade * D)
#
# PH_clade in {0,1} is edge membership in the Philippine clade; D in {0,1} is
# donor presence (D = prevalence in [0,1] for the unrelated control set). The
# loss rate q10 is untouched, so Beta > 0 reads as "contact makes Philippine
# languages more likely to ACQUIRE phonemes this donor has".
#
# Five models are compared by a leave-one-out predictive score per Philippine
# tip: the [1] baseline plus one model per donor.
#
# Base rates are the 1-regime MLEs of [1]'s winning rate structure (ER or ARD),
# refit here and then held fixed, so Beta is the only free parameter and the LRT
# against baseline has df = 1. [1] is consulted for ER-vs-ARD only — its
# er2_*/ard2_* Philippine-regime rates are deliberately unused.
#
# CAVEAT: D is confounded with phoneme commonness (correlation with
# austronesian_prevalence: unrelated +0.93, Spanish +0.74, Japanese +0.66,
# English +0.39), so the unrelated set acts as a null control rather than a donor,
# and the four Betas are not independent.
#
# Rates are not capped; saturated fits have an arbitrary magnitude but a
# meaningful q01:q10 ratio, and beta_excl_saturated reports their leverage.
#
# Inputs:  data/RUHLENdf_PH_evolutionary.csv, data/phoneme_evolutionary_tip_mapping.csv,
#          data/phoneme_evolutionary_analysis_tree.rds, data/phoneme_key_pnas.csv,
#          data/PHONEME_evolutionary_model_selection.csv
#
# Outputs: data/PHONEME_evolutionary_donor_scores.csv     (Philippine tips)
#          data/PHONEME_evolutionary_donor_betas.csv      (5 models)
#          data/PHONEME_evolutionary_donor_base_rates.csv (84 phonemes)
#
# Runtime: ~2 minutes.
# =============================================================================

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(phytools)
library(numDeriv)

LOGP_FLOOR   <- log(1e-12)  # guards the log score against a single ~0 probability
SATURATION_N <- 10          # expected transitions over the tree before "saturated"
BETA_RANGE   <- c(-5, 5)

# ── 1. Load data ─────────────────────────────────────────────

RUHLENdf <- read_csv(here("data", "RUHLENdf_PH_evolutionary.csv"), show_col_types = FALSE)

tip_mapping <- read_csv(here("data", "phoneme_evolutionary_tip_mapping.csv"), show_col_types = FALSE)
pnas_key    <- read_csv(here("data", "phoneme_key_pnas.csv"), show_col_types = FALSE) |>
  select(phoneme = phoneme_id, ipa, class)

model_selection <- read_csv(
  here("data", "PHONEME_evolutionary_model_selection.csv"),
  show_col_types = FALSE
)

# Tree pruning, conflict-tip resolution, and the Philippine regime are decided
# once in [0]_Phylogenetic_Tree.R and loaded here, not rebuilt.
tree_bundle     <- readRDS(here("data", "phoneme_evolutionary_analysis_tree.rds"))
analysis_tree   <- tree_bundle$analysis_tree
regime_tree     <- tree_bundle$regime_tree
philippine_tips <- tree_bundle$philippine_tips

tips <- tip_mapping |>
  filter(tip %in% analysis_tree$tip.label) |>
  left_join(RUHLENdf, by = c("iso6393", "ruhlen_language" = "language"))

phoneme_cols <- grep("^phoneme_", names(RUHLENdf), value = TRUE)

stopifnot(
  "tip_mapping should join 1:1 onto the analysis tree" = nrow(tips) == Ntip(analysis_tree),
  "every tip must resolve a Language_type"             = !any(is.na(tips$Language_type)),
  "phoneme matrix must be complete"                    = !any(is.na(tips[phoneme_cols])),
  "[0]'s Philippine tips must match tip-level Language_type" =
    setequal(philippine_tips, tips$tip[tips$Language_type == "Philippine Language"])
)

# ── 2. Likelihood-engine setup ────────────────────────────────
# reorder() preserves the simmap, so read the per-edge regime straight off $maps.

tree_po    <- reorder(regime_tree, "postorder")
is_ph_edge <- vapply(tree_po$maps, \(m) names(m)[1] == "Philippine", logical(1))

stopifnot(
  "reorder() must preserve the simmap"       = !is.null(tree_po$maps),
  "every edge must carry one regime segment" = all(lengths(tree_po$maps) == 1L)
)

n_tip   <- Ntip(tree_po)
n_node  <- tree_po$Nnode
root_id <- n_tip + 1L
parents <- unique(tree_po$edge[, 1])
mean_depth <- mean(node.depth.edgelength(analysis_tree)[seq_len(Ntip(analysis_tree))])

message(
  "Loaded analysis tree: ", n_tip, " tips (", length(philippine_tips), " Philippine), ",
  sum(is_ph_edge), " Philippine edges of ", nrow(tree_po$edge),
  "; mean root-to-tip depth ", round(mean_depth, 2)
)

# ── 3. Core likelihood ───────────────────────────────────────
# Postorder pruning for a 2-state trait, using the closed-form 2x2 transition
# probabilities (no matrix exponential). NA tips integrate out, which is what makes
# the leave-one-out score possible.

mk_loglik <- function(state_vec, q01_bg, q10_bg, q01_ph, q10_ph, pi = c(0.5, 0.5)) {
  L  <- matrix(0, n_tip + n_node, 2L)
  ok <- !is.na(state_vec)
  L[cbind(which(ok), state_vec[ok] + 1L)] <- 1
  L[which(!ok), ] <- 1

  log_scale <- 0

  for (nd in parents) {
    ii  <- which(tree_po$edge[, 1] == nd)
    acc <- c(1, 1)
    for (i in ii) {
      child <- tree_po$edge[i, 2]
      len   <- tree_po$edge.length[i]
      q01   <- if (is_ph_edge[i]) q01_ph else q01_bg
      q10   <- if (is_ph_edge[i]) q10_ph else q10_bg
      s <- q01 + q10
      e <- exp(-s * len)
      P <- matrix(
        c((q10 + q01 * e) / s, q01 * (1 - e) / s,
          q10 * (1 - e) / s,   (q01 + q10 * e) / s),
        nrow = 2L, byrow = TRUE
      )
      acc <- acc * (P %*% L[child, ])[, 1]
    }
    sc <- sum(acc)
    if (sc > 0) {
      L[nd, ]   <- acc / sc
      log_scale <- log_scale + log(sc)
    } else {
      L[nd, ] <- acc
    }
  }

  log_scale + log(sum(pi * L[root_id, ]))
}

state_of <- function(phoneme_col) {
  tips[[phoneme_col]][match(tree_po$tip.label, tips$tip)]
}

# ---- Validation gate: must match phytools before anything else runs ---------
# Also pins down the rate ordering: fitmultiMk uses [q10_bg, q01_bg, q10_ph, q01_ph].
validate_likelihood <- function() {
  probes     <- c("phoneme_141", "phoneme_012", "phoneme_294")
  rate_sets  <- list(c(0.2, 0.1, 0.5, 0.3), c(0.05, 0.05, 0.05, 0.05), c(1.5, 0.4, 0.9, 2.2))
  worst <- 0
  for (p in probes) {
    v  <- state_of(p)
    fv <- factor(v, levels = c(0, 1)); names(fv) <- tree_po$tip.label
    ref <- fitmultiMk(regime_tree, fv, model = "ARD")
    for (r in c(rate_sets, list(ref$rates))) {
      got <- mk_loglik(v, q01_bg = r[2], q10_bg = r[1], q01_ph = r[4], q10_ph = r[3])
      worst <- max(worst, abs(got - ref$lik(r)))
    }
  }
  message("Likelihood validation vs phytools::fitmultiMk: max abs diff = ", signif(worst, 3))
  stopifnot("custom likelihood does not match fitmultiMk" = worst < 1e-8)
}
validate_likelihood()

# ── 4. Phoneme set, states, and donor D vectors ──────────────

phonemes       <- model_selection$phoneme
rate_structure <- if_else(startsWith(model_selection$winning_model, "ARD"), "ARD", "ER")
states         <- lapply(phonemes, state_of)
n_phon         <- length(phonemes)

D <- list(
  span      = model_selection$span_has,
  eng       = model_selection$eng_has,
  jap       = model_selection$jap_has,
  unrelated = model_selection$unrelated_prevalence
)

message(n_phon, " phonemes; D=1 counts -- ",
        paste(names(D), vapply(D, \(d) sum(d > 0), integer(1)), sep = ": ", collapse = ", "))

# ── 5. Base rates, per phoneme ───────────────────────────────
# The 1-regime model is the constrained case q_ph == q_bg, fitted through the same
# engine as Beta so the baseline is exactly the Beta = 0 case of every donor model.
# Refit rather than reuse [1]'s ace() rates, which use a different root prior.

# Optimise on the log scale so positivity is automatic and the small and saturated
# optima are equally well conditioned.
LOG_RATE_RANGE <- c(log(1e-8), log(1e6))

fit_base_rates <- function(i) {
  v <- states[[i]]
  if (rate_structure[i] == "ER") {
    o <- optimize(\(lq) mk_loglik(v, exp(lq), exp(lq), exp(lq), exp(lq)),
                  interval = LOG_RATE_RANGE, maximum = TRUE)
    q <- exp(o$maximum)
    list(q01 = q, q10 = q, loglik = o$objective, start_used = "bounded_1d")
  } else {
    # Two starts: ape's ip = 0.1 default, and [1]'s ace estimates. The 2-D
    # surface is where [1] showed pathologies, so one start risks a local optimum.
    starts <- list(
      ip  = c(0.1, 0.1),
      ace = c(model_selection$ard_rate_01[i], model_selection$ard_rate_10[i])
    )
    best <- NULL
    for (nm in names(starts)) {
      st <- log(pmin(pmax(starts[[nm]], 1e-8), 1e6))
      o <- tryCatch(
        nlminb(st, \(lp) -mk_loglik(v, exp(lp[1]), exp(lp[2]), exp(lp[1]), exp(lp[2])),
               lower = rep(LOG_RATE_RANGE[1], 2), upper = rep(LOG_RATE_RANGE[2], 2)),
        error = function(e) NULL
      )
      if (!is.null(o) && (is.null(best) || -o$objective > best$loglik)) {
        best <- list(q01 = exp(o$par[1]), q10 = exp(o$par[2]),
                     loglik = -o$objective, start_used = nm)
      }
    }
    best
  }
}

base <- lapply(seq_len(n_phon), fit_base_rates)
q01_base <- vapply(base, \(b) b$q01, numeric(1))
q10_base <- vapply(base, \(b) b$q10, numeric(1))

base_rates <- tibble(
  phoneme              = phonemes,
  rate_structure       = rate_structure,
  q01_base             = q01_base,
  q10_base             = q10_base,
  stationary_p_present = q01_base / (q01_base + q10_base),
  loglik_1reg          = vapply(base, \(b) b$loglik, numeric(1)),
  start_used           = vapply(base, \(b) b$start_used, character(1)),
  base_rate_saturated  = pmax(q01_base, q10_base) * mean_depth > SATURATION_N
) |>
  left_join(pnas_key, by = "phoneme") |>
  left_join(model_selection |> select(phoneme, austronesian_prevalence), by = "phoneme") |>
  select(phoneme, ipa, class, rate_structure, q01_base, q10_base,
         stationary_p_present, austronesian_prevalence, loglik_1reg,
         base_rate_saturated, start_used)

message("Base rates fitted. Saturated: ", sum(base_rates$base_rate_saturated),
        " of ", n_phon, "; ARD starts won by ace: ",
        sum(base_rates$start_used == "ace"))

# ── 6. Beta, per donor ───────────────────────────────────────
# Joint objective: summed log-likelihood over phonemes, base rates fixed, one
# shared Beta scaling the Philippine gain rate in proportion to D.

joint_loglik <- function(beta, d, idx = seq_len(n_phon)) {
  sum(vapply(idx, function(i) {
    mk_loglik(
      states[[i]],
      q01_bg = q01_base[i], q10_bg = q10_base[i],
      q01_ph = q01_base[i] * exp(beta * d[i]), q10_ph = q10_base[i]
    )
  }, numeric(1)))
}

loglik_baseline <- joint_loglik(0, D$span)  # Beta = 0 -> any d gives the baseline

fit_beta <- function(name) {
  d   <- D[[name]]
  o   <- optimize(\(b) joint_loglik(b, d), BETA_RANGE, maximum = TRUE)
  se  <- tryCatch({
    h <- numDeriv::hessian(\(b) joint_loglik(b, d), o$maximum)
    suppressWarnings(sqrt(-1 / h[1, 1]))
  }, error = function(e) NA_real_)

  keep <- which(!base_rates$base_rate_saturated)
  o_ex <- optimize(\(b) joint_loglik(b, d, keep), BETA_RANGE, maximum = TRUE)

  lr <- 2 * (o$objective - loglik_baseline)
  tibble(
    model               = name,
    beta                = o$maximum,
    beta_se             = se,
    beta_excl_saturated = o_ex$maximum,
    n_phonemes_D1       = sum(d > 0),
    loglik              = o$objective,
    loglik_baseline     = loglik_baseline,
    LR                  = lr,
    df                  = 1L,
    p_value             = pchisq(pmax(lr, 0), df = 1, lower.tail = FALSE),
    beta_at_bound       = isTRUE(abs(o$maximum - BETA_RANGE[1]) < 1e-3 |
                                 abs(o$maximum - BETA_RANGE[2]) < 1e-3),
    se_na               = !is.finite(se)
  )
}

betas <- bind_rows(
  tibble(model = "austronesian", beta = 0, beta_se = NA_real_, beta_excl_saturated = 0,
         n_phonemes_D1 = NA_integer_, loglik = loglik_baseline,
         loglik_baseline = loglik_baseline, LR = 0, df = 0L, p_value = NA_real_,
         beta_at_bound = FALSE, se_na = NA),
  bind_rows(lapply(names(D), fit_beta))
)

message("Betas fitted:\n", paste0(
  "  ", betas$model, ": beta = ", round(betas$beta, 4),
  " (excl. saturated ", round(betas$beta_excl_saturated, 4), ")",
  ", p = ", signif(betas$p_value, 3), collapse = "\n"
))

# ── 7. Leave-one-out predictive scores ───────────────────────
# P(observed | all other tips) = L(tip=obs) / sum_s L(tip=s), computed on the log
# scale so it stays stable when one state dominates.

ph_tip_idx <- match(philippine_tips, tree_po$tip.label)

# Per-model Philippine gain rates; NULL beta means the baseline (no multiplier).
model_rates <- c(
  list(austronesian = q01_base),
  lapply(names(D), \(nm) q01_base * exp(betas$beta[betas$model == nm] * D[[nm]]))
)
names(model_rates) <- c("austronesian", names(D))

# logp[[model]][phoneme, tip]
logp <- lapply(model_rates, \(.) matrix(NA_real_, n_phon, length(ph_tip_idx)))

loo_logp <- function(v, q01_bg, q10_bg, q01_ph) {
  vapply(ph_tip_idx, function(ti) {
    l1 <- mk_loglik(replace(v, ti, 1L), q01_bg, q10_bg, q01_ph, q10_bg)
    l0 <- mk_loglik(replace(v, ti, 0L), q01_bg, q10_bg, q01_ph, q10_bg)
    lp <- if (v[ti] == 1L) -log1p(exp(l0 - l1)) else -log1p(exp(l1 - l0))
    if (!is.finite(lp)) LOGP_FLOOR else max(lp, LOGP_FLOOR)
  }, numeric(1))
}

for (i in seq_len(n_phon)) {
  v  <- states[[i]]
  qs <- vapply(model_rates, \(q) q[i], numeric(1))
  # Models sharing a Philippine gain rate share predictions, so evaluate only the
  # distinct rates.
  uq    <- unique(qs)
  cache <- lapply(uq, \(q01_ph) loo_logp(v, q01_base[i], q10_base[i], q01_ph))
  for (m in names(model_rates)) logp[[m]][i, ] <- cache[[match(qs[[m]], uq)]]

  if (i %% 10 == 0 || i == n_phon) message("  scored ", i, "/", n_phon, " phonemes")
}

floored_by_tip <- vapply(
  seq_along(ph_tip_idx),
  \(j) sum(vapply(logp, \(M) sum(M[, j] <= LOGP_FLOOR), integer(1))),
  integer(1)
)
floored <- sum(floored_by_tip)

# Macro-average on the log scale, then exponentiate to (0, 1] — a geometric mean,
# deliberately not mean(exp(logp)), which would not be a proper score.
balanced_score <- function(lp, present, idx = seq_along(present)) {
  lp <- lp[idx]; present <- present[idx]
  halves <- c(
    if (any(present))  mean(lp[present])  else NA_real_,
    if (any(!present)) mean(lp[!present]) else NA_real_
  )
  c(score   = exp(mean(halves, na.rm = TRUE)),
    present = exp(halves[1]),
    absent  = exp(halves[2]))
}

score_rows <- lapply(seq_along(ph_tip_idx), function(j) {
  ti      <- ph_tip_idx[j]
  present <- vapply(states, \(v) v[ti] == 1L, logical(1))

  out <- tibble(
    tip            = tree_po$tip.label[ti],
    n_present      = sum(present),
    n_absent       = sum(!present),
    n_logp_floored = floored_by_tip[j]
  )

  for (m in names(model_rates)) {
    b <- balanced_score(logp[[m]][, j], present)
    out[[paste0(m, "_score")]]         <- b[["score"]]
    out[[paste0(m, "_score_present")]] <- b[["present"]]
    out[[paste0(m, "_score_absent")]]  <- b[["absent"]]
  }
  # D=1 companions: each donor and the baseline on the SAME phoneme subset,
  # so the pair is directly comparable.
  for (nm in names(D)) {
    idx <- which(D[[nm]] > 0)
    out[[paste0(nm, "_score_D1")]] <- balanced_score(logp[[nm]][, j], present, idx)[["score"]]
    out[[paste0("austronesian_score_", nm, "_D1")]] <-
      balanced_score(logp[["austronesian"]][, j], present, idx)[["score"]]
  }
  out
})

donor_scores <- bind_rows(score_rows) |>
  left_join(
    tips |>
      add_count(ruhlen_language, name = "n_tips_for_language") |>
      select(tip, ruhlen_language, iso6393, glottocode, inventory_size,
             source, latitude, longitude, Language_type, n_tips_for_language),
    by = "tip"
  ) |>
  select(tip, ruhlen_language, iso6393, glottocode, inventory_size, source,
         latitude, longitude, Language_type, n_tips_for_language,
         austronesian_score, span_score, eng_score, jap_score, unrelated_score,
         ends_with("_score_present"), ends_with("_score_absent"),
         ends_with("_D1"), n_present, n_absent, n_logp_floored) |>
  arrange(tip)

# ── 8. Write outputs ─────────────────────────────────────────

write_csv(donor_scores, here("data", "PHONEME_evolutionary_donor_scores.csv"))
write_csv(betas,        here("data", "PHONEME_evolutionary_donor_betas.csv"))
write_csv(base_rates,   here("data", "PHONEME_evolutionary_donor_base_rates.csv"))

message(
  "\nDone.",
  "\n  donor scores : ", nrow(donor_scores), " Philippine tips -> data/PHONEME_evolutionary_donor_scores.csv",
  "\n  betas        : ", nrow(betas), " models -> data/PHONEME_evolutionary_donor_betas.csv",
  "\n  base rates   : ", nrow(base_rates), " phonemes -> data/PHONEME_evolutionary_donor_base_rates.csv",
  "\n  logp floored : ", floored, " of ", n_phon * length(ph_tip_idx) * length(model_rates)
)

print(
  donor_scores |>
    summarise(across(c(austronesian_score, span_score, eng_score, jap_score, unrelated_score),
                     \(x) round(mean(x), 4)))
)
