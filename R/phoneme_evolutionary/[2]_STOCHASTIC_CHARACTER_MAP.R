# =============================================================================
# [2] Phoneme Stochastic Character Maps
#
# Runs 1000 stochastic character maps per phoneme for the phonemes that won a
# two-regime model in [1] (ER2 or ARD2), and plots the resulting ancestral
# reconstruction on the pulses tree with the Philippine clade painted red.
#
# Inputs:  data/PHONEME_evolutionary_model_selection.csv,
#          data/RUHLENdf_PH_evolutionary.csv, data/phoneme_evolutionary_tip_mapping.csv,
#          data/phoneme_evolutionary_analysis_tree.rds
#
# Outputs: data/PHONEME_evolutionary_simmap_node_posteriors.csv
#          data/PHONEME_evolutionary_simmap_transitions.csv
#          data/PHONEME_evolutionary_simmap_maps.rds
#          figures/phoneme/evolutionary/simmap_<phoneme>_<model>.png
# =============================================================================

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(phytools)

NSIM           <- 1000
SEED           <- 20260727
TARGET_MODELS  <- c("ER2", "ARD2")
STATES         <- c("0", "1")
ROOT_PRIOR     <- "equal"   # matches the root prior fitMk/fitmultiMk used in [1]
PH_PLACEHOLDER <- "PH_CLADE"

# Figure geometry. The tree is 172 tips, so it has to be tall and narrow; at
# actual size this is roughly half a 1512 x 982 pt laptop screen.
FIG_W   <- 7.5
FIG_H   <- 10
FIG_DPI <- 300
FSIZE   <- 0.38   # tip label size, relative to the plotting default
PIE_CEX <- 0.17   # radius of both the node pies and the observed-state tip circles
TIP_OFFSET <- 1.2 # gap between the tip circles and the labels, in branch units

REGIME_COLORS <- c(Background = "black", Philippine = "red")
PIE_COLORS    <- c("white", "black")   # absent, present

# Toggle for the expensive single-Q equivalence check in section 6.
RUN_SPLIT_VALIDATION <- TRUE

set.seed(SEED)

# ── 1. Load data ─────────────────────────────────────────────

RUHLENdf <- read_csv(
  here("data", "RUHLENdf_PH_evolutionary.csv"),
  show_col_types = FALSE
)

tip_mapping <- read_csv(
  here("data", "phoneme_evolutionary_tip_mapping.csv"),
  show_col_types = FALSE
)

model_selection <- read_csv(
  here("data", "PHONEME_evolutionary_model_selection.csv"),
  show_col_types = FALSE
)

tree_bundle     <- readRDS(here("data", "phoneme_evolutionary_analysis_tree.rds"))
analysis_tree   <- tree_bundle$analysis_tree
regime_tree     <- tree_bundle$regime_tree
philippine_tips <- tree_bundle$philippine_tips

phoneme_cols <- grep("^phoneme_", names(RUHLENdf), value = TRUE)

# ── 2. Tip-level table ───────────────────────────────────────
# Same join and assertions as [1], so both scripts see identical tip states.

tips <- tip_mapping |>
  filter(tip %in% analysis_tree$tip.label) |>
  left_join(
    RUHLENdf,
    by = c("iso6393", "ruhlen_language" = "language")
  )

stopifnot(
  "tip_mapping should join 1:1 onto the analysis tree" =
    nrow(tips) == Ntip(analysis_tree),
  "phoneme matrix must be complete for all tips" =
    !any(is.na(tips[phoneme_cols])),
  "analysis tree must be binary" =
    is.binary(analysis_tree),
  "[0]'s Philippine tips must match tip-level Language_type" =
    setequal(philippine_tips, tips$tip[tips$Language_type == "Philippine Language"]),
  # Node pies and tip squares are drawn on regime_tree's coordinates but indexed
  # by analysis_tree's node numbers, so the two must be the same tree.
  "regime_tree must share analysis_tree's edge matrix" =
    identical(regime_tree$edge, analysis_tree$edge)
)

targets <- model_selection |>
  filter(winning_model %in% TARGET_MODELS) |>
  arrange(winning_model, phoneme)

stopifnot("no ER2/ARD2 phonemes found" = nrow(targets) > 0)

message(
  "Loaded analysis tree: ", Ntip(analysis_tree), " tips (",
  length(philippine_tips), " Philippine / ",
  Ntip(analysis_tree) - length(philippine_tips), " background).\n",
  nrow(targets), " two-regime phonemes to map (",
  paste(names(table(targets$winning_model)), table(targets$winning_model),
        sep = " x ", collapse = ", "), ")."
)

# ── 3. Rate matrices and the pruning algorithm ───────────────

# Q[i, j] is the instantaneous rate i -> j with states ordered ("0", "1"), i.e.
# Q[1, 2] is gain and Q[2, 1] is loss. Diagonals are set from the row sums so
# make.simmap never sees an invalid generator.
q_matrix <- function(rate_01, rate_10) {
  Q <- matrix(
    c(-rate_01, rate_01, rate_10, -rate_10),
    nrow = 2, byrow = TRUE, dimnames = list(STATES, STATES)
  )
  Q
}

# Regime-specific Q pair for one row of the model-selection table.
regime_q <- function(row) {
  if (row$winning_model == "ER2") {
    list(
      bg = q_matrix(row$er2_rate_bg, row$er2_rate_bg),
      ph = q_matrix(row$er2_rate_ph, row$er2_rate_ph)
    )
  } else {
    list(
      bg = q_matrix(row$ard2_rate_bg_01, row$ard2_rate_bg_10),
      ph = q_matrix(row$ard2_rate_ph_01, row$ard2_rate_ph_10)
    )
  }
}

# Closed-form transition probabilities for a 2-state generator, used instead of
# expm() to stay exact and dependency-free; validated against expm::expm in section 6.
p_matrix <- function(Q, t) {
  a <- Q[1, 2]
  b <- Q[2, 1]
  s <- a + b
  if (t == 0 || s == 0) return(diag(2))
  e <- exp(-s * t)
  matrix(
    c(b + a * e, a - a * e,
      b - b * e, a + b * e),
    nrow = 2, byrow = TRUE
  ) / s
}

# Felsenstein pruning for a binary character; tip_prior is one row per tip
# (rownames = tip labels), columns ("0", "1"). Node vectors are rescaled to sum 1
# with the discarded factors accumulated in logscale, so deep trees cannot
# underflow; both are returned, since section 6 needs the unnormalised likelihood.
partial_lik <- function(tree, tip_prior, Q) {
  n_tip <- Ntip(tree)
  L     <- matrix(1, n_tip + Nnode(tree), 2)
  L[seq_len(n_tip), ] <- tip_prior[tree$tip.label, , drop = FALSE]

  tr       <- reorder(tree, "postorder")
  logscale <- 0

  for (i in seq_len(nrow(tr$edge))) {
    parent <- tr$edge[i, 1]
    child  <- tr$edge[i, 2]
    L[parent, ] <- L[parent, ] * as.vector(p_matrix(Q, tr$edge.length[i]) %*% L[child, ])
    s <- sum(L[parent, ])
    if (s > 0) {
      L[parent, ] <- L[parent, ] / s
      logscale    <- logscale + log(s)
    }
  }

  list(lik = L[n_tip + 1L, ], logscale = logscale)
}

# Point-mass tip priors from an observed 0/1 vector.
prior_matrix <- function(states, labels) {
  m <- matrix(0, length(labels), 2, dimnames = list(labels, STATES))
  m[cbind(seq_along(labels), as.integer(as.character(states[labels])) + 1L)] <- 1
  m
}

state_vector <- function(phoneme_col) {
  v <- tips[[phoneme_col]][match(analysis_tree$tip.label, tips$tip)]
  v <- factor(v, levels = c(0, 1))
  names(v) <- analysis_tree$tip.label
  v
}

# ── 4. Split the tree at the Philippine MRCA ─────────────────
# Topology-only, so this is built once and reused for every phoneme.

philippine_node <- getMRCA(analysis_tree, philippine_tips)
ph_clade        <- extract.clade(analysis_tree, philippine_node)
stem_len        <- analysis_tree$edge.length[match(philippine_node, analysis_tree$edge[, 2])]

# drop.tip collapses the singleton path through the MRCA, so the surviving
# representative's terminal edge comes out as stem + path-to-that-tip. Overwrite
# it with the stem alone: on the background tree that edge *is* the clade's stem.
bg_tree <- drop.tip(analysis_tree, setdiff(philippine_tips, philippine_tips[1]))
bg_tree$tip.label[bg_tree$tip.label == philippine_tips[1]] <- PH_PLACEHOLDER
ph_edge_bg <- match(match(PH_PLACEHOLDER, bg_tree$tip.label), bg_tree$edge[, 2])
bg_tree$edge.length[ph_edge_bg] <- stem_len

n_bg_tip   <- Ntip(analysis_tree) - length(philippine_tips) + 1L
n_full_edge <- nrow(analysis_tree$edge)

stopifnot(
  "background tree should keep every non-Philippine tip plus one placeholder" =
    Ntip(bg_tree) == n_bg_tip,
  "clade should hold exactly the Philippine tips" =
    setequal(ph_clade$tip.label, philippine_tips),
  "split must conserve edges" =
    nrow(bg_tree$edge) + nrow(ph_clade$edge) == n_full_edge
)

# Descendant tip set of every node, used to line the two halves' edges back up
# with the full tree's.
descendant_tips <- function(tree) {
  n_tip <- Ntip(tree)
  out   <- vector("list", n_tip + Nnode(tree))
  out[seq_len(n_tip)] <- as.list(tree$tip.label)
  tr <- reorder(tree, "postorder")
  for (i in seq_len(nrow(tr$edge))) {
    out[[tr$edge[i, 1]]] <- c(out[[tr$edge[i, 1]]], out[[tr$edge[i, 2]]])
  }
  out
}

clade_key <- function(x) paste(sort(x), collapse = "\r")

full_desc <- descendant_tips(analysis_tree)
full_keys <- vapply(analysis_tree$edge[, 2], \(v) clade_key(full_desc[[v]]), character(1))

bg_desc <- descendant_tips(bg_tree)
bg_keys <- vapply(bg_tree$edge[, 2], function(v) {
  labs <- bg_desc[[v]]
  # The placeholder stands for the whole clade, so expand it before matching.
  if (PH_PLACEHOLDER %in% labs) labs <- c(setdiff(labs, PH_PLACEHOLDER), philippine_tips)
  clade_key(labs)
}, character(1))

ph_desc <- descendant_tips(ph_clade)
ph_keys <- vapply(ph_clade$edge[, 2], \(v) clade_key(ph_desc[[v]]), character(1))

bg_to_full <- match(bg_keys, full_keys)
ph_to_full <- match(ph_keys, full_keys)

stopifnot(
  "every background edge must map onto a full-tree edge" = !anyNA(bg_to_full),
  "every clade edge must map onto a full-tree edge"      = !anyNA(ph_to_full),
  "edge correspondence must be a bijection" =
    setequal(c(bg_to_full, ph_to_full), seq_len(n_full_edge)) &&
    !anyDuplicated(c(bg_to_full, ph_to_full)),
  "mapped edges must keep their lengths" =
    max(abs(c(
      bg_tree$edge.length - analysis_tree$edge.length[bg_to_full],
      ph_clade$edge.length - analysis_tree$edge.length[ph_to_full]
    ))) < 1e-8
)

ph_edge_full <- ph_to_full                       # edges inside the Philippine regime
bg_edge_full <- bg_to_full                       # everything else, incl. the stem
regime_of_edge <- rep("Background", n_full_edge)
regime_of_edge[ph_edge_full] <- "Philippine"

regime_branch_length <- c(
  Background = sum(analysis_tree$edge.length[bg_edge_full]),
  Philippine = sum(analysis_tree$edge.length[ph_edge_full])
)

internal_nodes <- seq(Ntip(analysis_tree) + 1L, Ntip(analysis_tree) + Nnode(analysis_tree))
node_first_edge <- match(internal_nodes, analysis_tree$edge[, 1])
node_n_desc <- vapply(internal_nodes, \(v) length(full_desc[[v]]), integer(1))
node_in_clade <- vapply(
  internal_nodes,
  \(v) all(full_desc[[v]] %in% philippine_tips),
  logical(1)
)

# ── 5. Regime-split stochastic mapping ───────────────────────
#
# make.simmap takes a single Q, so the two-regime model is mapped in two halves
# and spliced — exact because the likelihood factorises at the clade MRCA. This
# relies on make.simmap treating matrix x as a multiplicative prior, which
# checks 6b/6c test.

# make.simmap returns a bare simmap when nsim == 1 and a multiSimmap otherwise.
as_map_list <- function(x, n) if (n == 1L) list(x) else lapply(seq_len(n), \(i) x[[i]])

# State at the tipward end of an edge's map (phytools stores maps rootward first).
end_state <- function(m) names(m)[length(m)]

regime_simmap <- function(v, Q_bg, Q_ph, nsim = NSIM) {
  x_ph <- prior_matrix(v, ph_clade$tip.label)

  # Step 1: collapse the clade to a 2-vector of partial likelihoods at its root.
  pl_ph <- partial_lik(ph_clade, x_ph, Q_ph)

  # Step 2: background map, with that vector as the placeholder tip's prior.
  x_bg <- prior_matrix(v, setdiff(bg_tree$tip.label, PH_PLACEHOLDER))
  x_bg <- rbind(x_bg, matrix(pl_ph$lik, 1, 2, dimnames = list(PH_PLACEHOLDER, STATES)))
  x_bg <- x_bg[bg_tree$tip.label, , drop = FALSE]

  bg_maps <- as_map_list(
    make.simmap(bg_tree, x_bg, Q = Q_bg, nsim = nsim, pi = ROOT_PRIOR, message = FALSE),
    nsim
  )

  # Step 3: the state each draw sampled at the Philippine MRCA.
  root_states <- vapply(bg_maps, \(m) end_state(m$maps[[ph_edge_bg]]), character(1))

  # Clade histories are i.i.d. given that state, and there are only two possible
  # states, so batch into two calls and deal the draws back out in order rather
  # than making 1000 separate conditional calls.
  ph_maps <- vector("list", nsim)
  for (s in STATES) {
    idx <- which(root_states == s)
    if (!length(idx)) next
    pi_s <- setNames(as.numeric(STATES == s), STATES)
    ph_maps[idx] <- as_map_list(
      make.simmap(ph_clade, x_ph, Q = Q_ph, nsim = length(idx), pi = pi_s, message = FALSE),
      length(idx)
    )
  }

  # Step 4: splice the two halves onto the full tree.
  stitched <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    maps <- vector("list", n_full_edge)
    maps[bg_edge_full] <- bg_maps[[i]]$maps
    maps[ph_edge_full] <- ph_maps[[i]]$maps

    full <- analysis_tree
    full$maps <- maps
    full$mapped.edge <- t(vapply(maps, function(m) {
      out <- setNames(numeric(2), STATES)
      agg <- tapply(m, names(m), sum)
      out[names(agg)] <- agg
      out
    }, numeric(2)))
    rownames(full$mapped.edge) <- paste(analysis_tree$edge[, 1], ",",
                                        analysis_tree$edge[, 2], sep = "")
    class(full) <- c("simmap", "phylo")
    stitched[[i]] <- full
  }

  stopifnot(
    "stitched edge times must equal the full tree's branch lengths" =
      max(abs(vapply(stitched[[1]]$maps, sum, numeric(1)) - analysis_tree$edge.length)) < 1e-8
  )

  list(maps = stitched, root_states = root_states, pl_ph = pl_ph, x_bg = x_bg)
}

# Posterior P(state = 1) at every internal node, read off the start of each
# node's first outgoing edge.
node_posterior <- function(maps) {
  states <- vapply(
    maps,
    \(m) vapply(node_first_edge, \(e) names(m$maps[[e]])[1], character(1)),
    character(length(internal_nodes))
  )
  rowMeans(states == "1")
}

# Gains (0 -> 1) and losses (1 -> 0) per draw, split by regime.
transition_counts <- function(maps) {
  out <- matrix(
    0L, length(maps), 4,
    dimnames = list(NULL, c("Background_01", "Background_10",
                            "Philippine_01", "Philippine_10"))
  )
  for (i in seq_along(maps)) {
    m   <- maps[[i]]$maps
    hit <- which(lengths(m) > 1L)
    for (e in hit) {
      nm  <- names(m[[e]])
      key <- paste0(regime_of_edge[e], "_", nm[-length(nm)], nm[-1])
      tb  <- table(key)
      out[i, names(tb)] <- out[i, names(tb)] + as.integer(tb)
    }
  }
  out
}

# ── 6. Verification ──────────────────────────────────────────

# 6a. Closed-form P(t) against a general matrix exponential.
local({
  set.seed(1)
  worst <- 0
  for (rep in 1:50) {
    Q <- q_matrix(runif(1, 0, 20), runif(1, 0, 20))
    for (t in c(0, 1e-6, 0.1, 1, 4.5)) {
      worst <- max(worst, max(abs(p_matrix(Q, t) - as.matrix(expm::expm(Q * t)))))
    }
  }
  stopifnot("closed-form P(t) must match expm::expm" = worst < 1e-10)
  message("check 6a: closed-form P(t) matches expm to ", signif(worst, 3), ".")
})

# 6b. make.simmap must treat a matrix x as a multiplicative prior (not a fixed
# state) and honour a point-mass pi. If either fails, the split is invalid.
local({
  v      <- state_vector(targets$phoneme[1])
  Q      <- q_matrix(0.3, 0.3)
  x_soft <- prior_matrix(v, bg_tree$tip.label[bg_tree$tip.label != PH_PLACEHOLDER])
  x_soft <- rbind(x_soft, matrix(c(0.5, 0.5), 1, 2,
                                 dimnames = list(PH_PLACEHOLDER, STATES)))
  x_soft <- x_soft[bg_tree$tip.label, , drop = FALSE]

  probe <- as_map_list(
    make.simmap(bg_tree, x_soft, Q = Q, nsim = 50, pi = ROOT_PRIOR, message = FALSE), 50
  )
  varies <- length(unique(vapply(probe, \(m) end_state(m$maps[[ph_edge_bg]]), character(1)))) > 1

  pinned <- as_map_list(
    make.simmap(ph_clade, prior_matrix(v, ph_clade$tip.label), Q = Q, nsim = 20,
                pi = setNames(c(1, 0), STATES), message = FALSE), 20
  )
  root_edge <- which(pinned[[1]]$edge[, 1] == Ntip(ph_clade) + 1L)[1]
  pinned_ok <- all(vapply(pinned, \(m) names(m$maps[[root_edge]])[1], character(1)) == "0")

  stopifnot(
    "make.simmap must sample, not fix, a tip with a 50/50 prior" = varies,
    "make.simmap must honour a point-mass pi at the root"        = pinned_ok
  )
  message("check 6b: tip priors are multiplicative and point-mass pi is honoured.")
})

# 6c. With Q_bg == Q_ph the split must be invisible, so split-and-stitch node
# posteriors have to match a plain single-Q map — the one check that exercises
# tip priors, root sampling, clade draws, stem continuity and stitching together.
if (RUN_SPLIT_VALIDATION) local({
  probe_row <- targets |> filter(phoneme == "phoneme_123")
  if (!nrow(probe_row)) probe_row <- targets[1, ]
  v  <- state_vector(probe_row$phoneme)
  Q1 <- q_matrix(probe_row$er_rate, probe_row$er_rate)

  split_post <- node_posterior(regime_simmap(v, Q1, Q1, nsim = NSIM)$maps)
  plain <- as_map_list(
    make.simmap(analysis_tree, v, Q = Q1, nsim = NSIM, pi = ROOT_PRIOR, message = FALSE),
    NSIM
  )
  plain_post <- node_posterior(plain)

  # Both estimates are binomial proportions over NSIM draws, so the tolerance is
  # set on the pooled two-proportion z scale (+1 smoothed for nodes pinned at 0
  # or 1) rather than in fixed posterior units. Mean |delta| catches a systematic
  # shift; max |z| catches a single badly-placed node.
  d      <- abs(split_post - plain_post)
  p_pool <- (split_post * NSIM + plain_post * NSIM + 1) / (2 * NSIM + 2)
  z      <- d / sqrt(2 * p_pool * (1 - p_pool) / NSIM)

  message(
    "check 6c: split vs single-Q node posteriors on ", probe_row$phoneme,
    " at nsim = ", NSIM, " -- mean |delta| = ", round(mean(d), 4),
    ", max |delta| = ", round(max(d), 4), ", max |z| = ", round(max(z), 2),
    " (", sum(z > 3), "/", length(z), " nodes over |z| = 3)."
  )
  stopifnot(
    "regime split must be invisible when Q_bg == Q_ph (systematic shift)" =
      mean(d) < 0.02,
    "regime split must be invisible when Q_bg == Q_ph (outlier node)" =
      max(z) < 5
  )
})

# 6d. Exact posterior marginal at the Philippine MRCA under the real
# Q_bg != Q_ph, versus what the background maps sampled. Computed from the
# pruning algorithm, so it carries no Monte Carlo error.
exact_root_marginal <- function(v, Q_bg, Q_ph) {
  pl_ph <- partial_lik(ph_clade, prior_matrix(v, ph_clade$tip.label), Q_ph)
  x_obs <- prior_matrix(v, setdiff(bg_tree$tip.label, PH_PLACEHOLDER))

  joint <- vapply(seq_along(STATES), function(s) {
    pinned <- rbind(
      x_obs,
      matrix(as.numeric(seq_along(STATES) == s), 1, 2,
             dimnames = list(PH_PLACEHOLDER, STATES))
    )[bg_tree$tip.label, , drop = FALSE]
    pl_bg <- partial_lik(bg_tree, pinned, Q_bg)
    # log P(clade data | MRCA = s) + log P(background data, MRCA = s)
    pl_ph$logscale + log(pl_ph$lik[s]) +
      pl_bg$logscale + log(sum(pl_bg$lik / length(STATES)))
  }, numeric(1))

  p <- exp(joint - max(joint))
  setNames(p / sum(p), STATES)
}

# ── 7. Map every two-regime phoneme ──────────────────────────

results <- vector("list", nrow(targets))
names(results) <- targets$phoneme
t_start <- Sys.time()

for (i in seq_len(nrow(targets))) {
  row <- targets[i, ]
  v   <- state_vector(row$phoneme)
  Q   <- regime_q(row)

  fit <- regime_simmap(v, Q$bg, Q$ph, nsim = NSIM)

  # check 6d, per phoneme: sampled root-state frequency vs its exact marginal.
  exact  <- exact_root_marginal(v, Q$bg, Q$ph)
  observed <- mean(fit$root_states == "1")
  # Clamped because a marginal of exactly 0 or 1 is a valid outcome for a phoneme
  # whose MRCA state is certain, and binom.test rejects p on the boundary.
  ci <- binom.test(
    sum(fit$root_states == "1"), NSIM,
    min(max(exact[["1"]], 1e-9), 1 - 1e-9)
  )$p.value

  results[[i]] <- list(
    phoneme    = row$phoneme,
    posterior  = node_posterior(fit$maps),
    trans      = transition_counts(fit$maps),
    root_exact = exact[["1"]],
    root_obs   = observed,
    root_p     = ci,
    maps       = fit$maps
  )

  message(
    "  [", i, "/", nrow(targets), "] ", row$phoneme, " (", row$ipa, ", ",
    row$winning_model, ") -- P(MRCA present): exact ", round(exact[["1"]], 3),
    " vs sampled ", round(observed, 3), " (p = ", signif(ci, 2), "); ",
    round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min"
  )
}

bad_marginal <- vapply(results, \(r) r$root_p, numeric(1)) < 0.001
if (any(bad_marginal)) {
  stop(
    "check 6d failed -- sampled MRCA states depart from their exact marginal for: ",
    paste(names(results)[bad_marginal], collapse = ", ")
  )
}
message("check 6d: sampled MRCA states match their exact marginals for all ",
        nrow(targets), " phonemes.")

# ── 8. Summary tables ────────────────────────────────────────

node_posteriors <- bind_rows(lapply(results, function(r) {
  tibble(
    phoneme            = r$phoneme,
    node               = internal_nodes,
    p_present          = r$posterior,
    n_desc_tips        = node_n_desc,
    in_philippine_clade = node_in_clade
  )
})) |>
  left_join(select(targets, phoneme, ipa, class, winning_model), by = "phoneme") |>
  select(phoneme, ipa, class, winning_model, node, p_present,
         n_desc_tips, in_philippine_clade)

transitions <- bind_rows(lapply(results, function(r) {
  as_tibble(r$trans) |>
    pivot_longer(everything(), names_to = "key", values_to = "n") |>
    separate(key, into = c("regime", "direction"), sep = "_") |>
    group_by(regime, direction) |>
    summarise(
      mean   = mean(n),
      median = median(n),
      q025   = quantile(n, 0.025),
      q975   = quantile(n, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      phoneme             = r$phoneme,
      direction           = if_else(direction == "01", "gain (0->1)", "loss (1->0)"),
      total_branch_length = regime_branch_length[regime],
      per_unit_time       = mean / total_branch_length
    )
})) |>
  left_join(select(targets, phoneme, ipa, class, winning_model), by = "phoneme") |>
  select(phoneme, ipa, class, winning_model, regime, direction,
         mean, median, q025, q975, total_branch_length, per_unit_time) |>
  arrange(phoneme, regime, direction)

write_csv(node_posteriors, here("data", "PHONEME_evolutionary_simmap_node_posteriors.csv"))
write_csv(transitions, here("data", "PHONEME_evolutionary_simmap_transitions.csv"))
saveRDS(
  lapply(results, \(r) r$maps),
  here("data", "PHONEME_evolutionary_simmap_maps.rds")
)

# ── 9. Figures ───────────────────────────────────────────────

fig_dir <- here("figures", "phoneme", "evolutionary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(targets))) {
  row <- targets[i, ]
  r   <- results[[i]]
  obs <- as.integer(as.character(state_vector(row$phoneme)))

  pies <- cbind(1 - r$posterior, r$posterior)

  png(
    file.path(fig_dir, paste0("simmap_", row$phoneme, "_", row$winning_model, ".png")),
    width = FIG_W, height = FIG_H, units = "in", res = FIG_DPI
  )

  plotSimmap(
    regime_tree,
    colors = REGIME_COLORS,
    fsize  = FSIZE,
    ftype  = "reg",
    lwd    = 1.0,
    offset = TIP_OFFSET,
    mar    = c(0.4, 0.2, 2.2, 0.2)
  )
  nodelabels(pie = pies, piecol = PIE_COLORS, cex = PIE_CEX)
  # Tip states are drawn as pies, not pch symbols: tiplabels/nodelabels size pies
  # by the same radius formula, so one cex makes tip and node markers match by
  # construction, whereas pch cex is a font-size scale.
  tiplabels(pie = cbind(1 - obs, obs), piecol = PIE_COLORS, cex = PIE_CEX)

  title(main = sprintf("/%s/ — %s", row$ipa, row$class), cex.main = 1.1)

  dev.off()
}

# ── 10. Report ───────────────────────────────────────────────
# Transitions per unit branch length inside vs outside the Philippine clade.
# ph_over_bg will not equal the fitted rate ratios: a mapped count is rate x time
# occupying the originating state, and occupancy differs between regimes.

rate_check <- transitions |>
  select(phoneme, ipa, winning_model, regime, direction, per_unit_time) |>
  pivot_wider(names_from = regime, values_from = per_unit_time) |>
  mutate(ph_over_bg = Philippine / Background)

message(
  "\nDone in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min.",
  "\n  node posteriors : ", nrow(node_posteriors), " rows -> ",
  "data/PHONEME_evolutionary_simmap_node_posteriors.csv",
  "\n  transitions     : ", nrow(transitions), " rows -> ",
  "data/PHONEME_evolutionary_simmap_transitions.csv",
  "\n  cached maps     : data/PHONEME_evolutionary_simmap_maps.rds",
  "\n  figures         : ", nrow(targets), " PNGs -> figures/phoneme/evolutionary/\n"
)

print(as.data.frame(rate_check))
