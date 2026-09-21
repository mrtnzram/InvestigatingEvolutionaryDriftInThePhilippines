# =============================================================================
# [2] Grammar Stochastic Character Maps
#
# Runs 1000 stochastic character maps per target GRAMBANK feature and plots the
# resulting ancestral reconstruction on the pulses tree with the Philippine
# clade painted red. Targets are the two multistate features plus every feature
# that won a two-regime model in [1].
#
# Characters are MULTISTATE throughout: k comes from each feature's own state
# set, transition probabilities are a general matrix exponential rather than the
# closed-form two-state formula, and node posteriors, transition counts and pies
# are all k-wide. Missing tip values stay ambiguous (an all-ones prior row) and
# are drawn with no tip marker at all.
#
# Inputs:  data/evolutionary/ml_model/GRAMMAR_evolutionary_model_selection.csv,
#          data/evolutionary/ml_model/GRAMMAR_evolutionary_model_rates.csv,
#          data/grammar/grammar_evolutionary_tip_matrix.csv,
#          data/grammar/grammar_evolutionary_analysis_tree.rds
#
# Outputs: data/evolutionary/stochastic_character_map/GRAMMAR_evolutionary_simmap_node_posteriors.csv
#          data/evolutionary/stochastic_character_map/GRAMMAR_evolutionary_simmap_transitions.csv
#          data/evolutionary/stochastic_character_map/GRAMMAR_evolutionary_simmap_maps.rds
#          figures/evolutionary/simmap_<feature>_<model>.png
# =============================================================================

library(ape)
library(dplyr)
library(tidyr)
library(readr)
library(here)
library(phytools)
library(expm)

source(here("R", "shared", "grambank_features.R"))

NSIM           <- 1000
SEED           <- 20260727
ALWAYS_MAP     <- c("GB065", "GB130")   # the multistate features, mapped whatever [1] chose
ROOT_PRIOR     <- "equal"   # matches the root prior fitMk/fitmultiMk used in [1]
PH_PLACEHOLDER <- "PH_CLADE"

# Figure geometry. The tree is 151 tips; laid out as a fan the tips spread
# around a circle instead of down a column, so the canvas is wide/square
# rather than tall and narrow.
FIG_W   <- 9
FIG_H   <- 8.5
FIG_DPI <- 300
FSIZE   <- 0.38   # tip label size, relative to the plotting default
PIE_CEX <- 0.17   # radius of both the node pies and the observed-state tip pies
TIP_OFFSET <- 1.2 # gap between the tip pies and the labels, in branch units

REGIME_COLORS <- c(Background = "black", Philippine = "red")

# Sliced to k and used for BOTH node and tip pies, so the two match by
# construction and the figure needs no legend. Colours key on a state's position
# in the feature's own sorted state set, never on the label: GB065 and GB130 run
# {1,2,3} with no 0 state.
STATE_COLORS <- c("white", "#1B9E77", "#7570B3")

# The stitched maps are ~10 MB per feature; set FALSE to skip persisting them.
SAVE_MAPS <- TRUE

# Toggle for the expensive single-Q equivalence check in section 6.
RUN_SPLIT_VALIDATION <- TRUE

set.seed(SEED)

# ── 1. Load data ─────────────────────────────────────────────

tip_matrix <- read_csv(
  here("data", "grammar", "grammar_evolutionary_tip_matrix.csv"),
  show_col_types = FALSE
)

model_selection <- read_csv(
  here("data", "evolutionary", "ml_model", "GRAMMAR_evolutionary_model_selection.csv"),
  show_col_types = FALSE
)

# from_state / to_state are state LABELS, not numbers: "0"/"1" would otherwise
# parse as numeric and stop matching the state vectors they index.
model_rates <- read_csv(
  here("data", "evolutionary", "ml_model", "GRAMMAR_evolutionary_model_rates.csv"),
  col_types = cols(from_state = col_character(), to_state = col_character(), .default = col_guess())
)

tree_bundle     <- readRDS(here("data", "grammar", "grammar_evolutionary_analysis_tree.rds"))
analysis_tree   <- tree_bundle$analysis_tree
regime_tree     <- tree_bundle$regime_tree
philippine_tips <- tree_bundle$philippine_tips
regime_tips     <- tree_bundle$regime_tips

stopifnot(
  "tip matrix should be 1:1 with the analysis tree" =
    nrow(tip_matrix) == Ntip(analysis_tree) &&
    setequal(tip_matrix$tip, analysis_tree$tip.label),
  "analysis tree must be binary" = is.binary(analysis_tree),
  # Node pies are drawn on regime_tree's coordinates but indexed by
  # analysis_tree's node numbers, so the two must be the same tree.
  "regime_tree must share analysis_tree's edge matrix" =
    identical(regime_tree$edge, analysis_tree$edge),
  "the regime clade must contain every Philippine tip" =
    all(philippine_tips %in% regime_tips)
)

targets <- model_selection |>
  filter(feature %in% ALWAYS_MAP | grepl("2$", winning_model)) |>
  arrange(desc(k), feature)

stopifnot("no target features found" = nrow(targets) > 0)

message(
  "Loaded analysis tree: ", Ntip(analysis_tree), " tips (",
  length(philippine_tips), " Philippine / ",
  Ntip(analysis_tree) - length(philippine_tips), " background); regime clade ",
  length(regime_tips), " tips.\n",
  nrow(targets), " features to map: ",
  paste0(targets$feature, " (", targets$winning_model, ", k=", targets$k, ")", collapse = ", ")
)

# ── 2. Feature helpers ───────────────────────────────────────

feature_values <- function(feature) {
  v <- tip_matrix[[feature]][match(analysis_tree$tip.label, tip_matrix$tip)]
  names(v) <- analysis_tree$tip.label
  v
}

feature_states <- function(feature) {
  strsplit(model_selection$states[model_selection$feature == feature], "/")[[1]]
}

# Point-mass tip priors where the state is observed and kept, an all-ones row
# where it is missing or was trimmed in [1]. Both fitMk and make.simmap consume
# the row multiplicatively, so an all-ones row is an uninformative prior.
prior_matrix <- function(v, labels, states) {
  vv <- v[labels]
  X  <- matrix(0, length(labels), length(states), dimnames = list(labels, states))
  obs <- !is.na(vv) & as.character(vv) %in% states
  X[!obs, ] <- 1
  X[cbind(which(obs), match(as.character(vv[obs]), states))] <- 1
  X
}

# ── 3. Rate matrices and the pruning algorithm ───────────────

# Q from [1]'s long rate table. A 1-regime winner emits a single "Single" block,
# so both regimes get the same generator and the split below collapses to the
# plain single-Q case.
regime_q <- function(feature, states) {
  rr <- model_rates |> filter(feature == !!feature)
  build <- function(sub) {
    Q <- matrix(0, length(states), length(states), dimnames = list(states, states))
    Q[cbind(match(sub$from_state, states), match(sub$to_state, states))] <- sub$rate
    diag(Q) <- -rowSums(Q)
    Q
  }
  if (all(rr$regime == "Single")) {
    Q <- build(rr)
    list(bg = Q, ph = Q)
  } else {
    list(
      bg = build(rr |> filter(regime == "Background")),
      ph = build(rr |> filter(regime == "Philippine"))
    )
  }
}

# General matrix exponential rather than the two-state closed form. partial_lik
# calls this a few thousand times per feature on a k x k matrix with k <= 3, not
# once per simulation, so the cost is immaterial next to the fitting in [1] and
# it stays correct for any k.
p_matrix <- function(Q, t) {
  if (t == 0) return(diag(nrow(Q)))
  as.matrix(expm::expm(Q * t))
}

# Felsenstein pruning; tip_prior is one row per tip (rownames = tip labels),
# one column per state. Node vectors are rescaled to sum 1 with the discarded
# factors accumulated in logscale, so deep trees cannot underflow; both are
# returned, since section 6 needs the unnormalised likelihood.
partial_lik <- function(tree, tip_prior, Q) {
  k     <- ncol(Q)
  n_tip <- Ntip(tree)
  L     <- matrix(1, n_tip + Nnode(tree), k)
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

# Log-likelihood of a whole tree under one Q with an equal root prior.
tree_loglik <- function(tree, tip_prior, Q) {
  pl <- partial_lik(tree, tip_prior, Q)
  pl$logscale + log(sum(pl$lik / ncol(Q)))
}

# ── 4. Split the tree at the Philippine MRCA ─────────────────
# Topology-only, so this is built once and reused for every feature.
#
# The split keys on `regime_tips`, the clade paintSubTree actually painted, NOT
# on `philippine_tips`. The two differ by the tips the MRCA swept in, and
# splitting on the smaller set would leave a swept-in tip on the background tree
# and inside the clade at once, double-counting its data.

philippine_node <- getMRCA(analysis_tree, philippine_tips)
ph_clade        <- extract.clade(analysis_tree, philippine_node)
stem_len        <- analysis_tree$edge.length[match(philippine_node, analysis_tree$edge[, 2])]

# drop.tip collapses the singleton path through the MRCA, so the surviving
# representative's terminal edge comes out as stem + path-to-that-tip. Overwrite
# it with the stem alone: on the background tree that edge *is* the clade's stem.
bg_tree <- drop.tip(analysis_tree, setdiff(regime_tips, regime_tips[1]))
bg_tree$tip.label[bg_tree$tip.label == regime_tips[1]] <- PH_PLACEHOLDER
ph_edge_bg <- match(match(PH_PLACEHOLDER, bg_tree$tip.label), bg_tree$edge[, 2])
bg_tree$edge.length[ph_edge_bg] <- stem_len

n_bg_tip    <- Ntip(analysis_tree) - length(regime_tips) + 1L
n_full_edge <- nrow(analysis_tree$edge)

stopifnot(
  "background tree should keep every non-clade tip plus one placeholder" =
    Ntip(bg_tree) == n_bg_tip,
  "clade should hold exactly the regime tips" =
    setequal(ph_clade$tip.label, regime_tips),
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
  if (PH_PLACEHOLDER %in% labs) labs <- c(setdiff(labs, PH_PLACEHOLDER), regime_tips)
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

internal_nodes  <- seq(Ntip(analysis_tree) + 1L, Ntip(analysis_tree) + Nnode(analysis_tree))
node_first_edge <- match(internal_nodes, analysis_tree$edge[, 1])
node_n_desc     <- vapply(internal_nodes, \(v) length(full_desc[[v]]), integer(1))
node_in_clade   <- vapply(
  internal_nodes,
  \(v) all(full_desc[[v]] %in% regime_tips),
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

regime_simmap <- function(v, states, Q_bg, Q_ph, nsim = NSIM) {
  k    <- length(states)
  x_ph <- prior_matrix(v, ph_clade$tip.label, states)

  # Step 1: collapse the clade to a k-vector of partial likelihoods at its root.
  pl_ph <- partial_lik(ph_clade, x_ph, Q_ph)

  # Step 2: background map, with that vector as the placeholder tip's prior.
  x_bg <- prior_matrix(v, setdiff(bg_tree$tip.label, PH_PLACEHOLDER), states)
  x_bg <- rbind(x_bg, matrix(pl_ph$lik, 1, k, dimnames = list(PH_PLACEHOLDER, states)))
  x_bg <- x_bg[bg_tree$tip.label, , drop = FALSE]

  bg_maps <- as_map_list(
    make.simmap(bg_tree, x_bg, Q = Q_bg, nsim = nsim, pi = ROOT_PRIOR, message = FALSE),
    nsim
  )

  # Step 3: the state each draw sampled at the Philippine MRCA.
  root_states <- vapply(bg_maps, \(m) end_state(m$maps[[ph_edge_bg]]), character(1))

  # Clade histories are i.i.d. given that state, so batch into one call per state
  # and deal the draws back out in order rather than making nsim conditional calls.
  ph_maps <- vector("list", nsim)
  for (s in states) {
    idx <- which(root_states == s)
    if (!length(idx)) next
    pi_s <- setNames(as.numeric(states == s), states)
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
      out <- setNames(numeric(k), states)
      agg <- tapply(m, names(m), sum)
      out[names(agg)] <- agg
      out
    }, numeric(k)))
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

# Posterior over states at every internal node, read off the start of each
# node's first outgoing edge. Returns a node x k matrix.
node_posterior <- function(maps, states) {
  st <- vapply(
    maps,
    \(m) vapply(node_first_edge, \(e) names(m$maps[[e]])[1], character(1)),
    character(length(internal_nodes))
  )
  out <- t(apply(st, 1, \(r) tabulate(match(r, states), length(states)) / length(maps)))
  dimnames(out) <- list(NULL, states)
  out
}

# Transitions per draw, split by regime and by ordered state pair. The key is
# "|"-separated because state labels are not guaranteed to be single characters
# the way the phoneme script's fixed "01"/"10" directions were.
transition_counts <- function(maps, states) {
  keys <- expand_grid(regime = names(REGIME_COLORS), from_state = states, to_state = states) |>
    filter(from_state != to_state) |>
    mutate(key = paste(regime, from_state, to_state, sep = "|")) |>
    pull(key)

  out <- matrix(0L, length(maps), length(keys), dimnames = list(NULL, keys))
  for (i in seq_along(maps)) {
    m   <- maps[[i]]$maps
    hit <- which(lengths(m) > 1L)
    for (e in hit) {
      nm  <- names(m[[e]])
      key <- paste(regime_of_edge[e], nm[-length(nm)], nm[-1], sep = "|")
      tb  <- table(key)
      out[i, names(tb)] <- out[i, names(tb)] + as.integer(tb)
    }
  }
  out
}

# ── 6. Verification ──────────────────────────────────────────

# 6a. The split likelihood machinery against phytools' own. fitMk with a fixed Q
# evaluates the likelihood without optimising, so this is cheap. One assertion
# covers p_matrix, partial_lik, prior_matrix, regime_q, the equal root prior,
# and — via the second comparison — the clade split and stem assignment.
local({
  worst_phy <- 0
  worst_split <- 0
  for (i in seq_len(nrow(targets))) {
    feature <- targets$feature[i]
    states  <- feature_states(feature)
    v       <- feature_values(feature)
    X       <- prior_matrix(v, analysis_tree$tip.label, states)
    Q       <- regime_q(feature, states)$bg

    mine <- tree_loglik(analysis_tree, X, Q)
    phy  <- as.numeric(logLik(fitMk(analysis_tree, X, fixedQ = Q)))
    worst_phy <- max(worst_phy, abs(mine - phy))

    # With one Q on both halves the split must be invisible at the likelihood level.
    pl_ph <- partial_lik(ph_clade, X[ph_clade$tip.label, , drop = FALSE], Q)
    x_bg  <- rbind(
      X[setdiff(bg_tree$tip.label, PH_PLACEHOLDER), , drop = FALSE],
      matrix(pl_ph$lik, 1, length(states), dimnames = list(PH_PLACEHOLDER, states))
    )[bg_tree$tip.label, , drop = FALSE]
    pl_bg <- partial_lik(bg_tree, x_bg, Q)
    split <- pl_ph$logscale + pl_bg$logscale + log(sum(pl_bg$lik / length(states)))
    worst_split <- max(worst_split, abs(split - mine))
  }
  stopifnot(
    "own likelihood must match fitMk at a fixed Q" = worst_phy < 1e-6,
    "split likelihood must equal the whole-tree likelihood under one Q" = worst_split < 1e-6
  )
  message(
    "check 6a: likelihood matches fitMk to ", signif(worst_phy, 3),
    "; split matches whole-tree to ", signif(worst_split, 3), "."
  )
})

# 6e. The edge-level regime must agree with the painted regime tree. A
# descendant-set predicate would also claim the clade's stem edge, which
# paintSubTree(stem = FALSE) assigns to Background; neither 6c nor 6d would
# catch that, because both are self-consistent under either convention.
local({
  painted <- unname(vapply(regime_tree$maps, \(z) names(z)[1], character(1)))
  stopifnot(
    "regime_of_edge must match the painted regime tree" =
      identical(regime_of_edge, painted)
  )
  message(
    "check 6e: edge regimes match the painted tree (",
    sum(regime_of_edge == "Philippine"), " Philippine / ",
    sum(regime_of_edge == "Background"), " Background edges)."
  )
})

# 6b. make.simmap must treat a matrix x as a multiplicative prior (not a fixed
# state) and honour a point-mass pi. If either fails, the split is invalid.
local({
  feature <- targets$feature[1]
  states  <- feature_states(feature)
  k       <- length(states)
  v       <- feature_values(feature)
  Q       <- regime_q(feature, states)$bg

  x_soft <- prior_matrix(v, setdiff(bg_tree$tip.label, PH_PLACEHOLDER), states)
  x_soft <- rbind(x_soft, matrix(1 / k, 1, k, dimnames = list(PH_PLACEHOLDER, states)))
  x_soft <- x_soft[bg_tree$tip.label, , drop = FALSE]

  probe <- as_map_list(
    make.simmap(bg_tree, x_soft, Q = Q, nsim = 50, pi = ROOT_PRIOR, message = FALSE), 50
  )
  varies <- length(unique(vapply(probe, \(m) end_state(m$maps[[ph_edge_bg]]), character(1)))) > 1

  x_ph <- prior_matrix(v, ph_clade$tip.label, states)
  pinned_ok <- all(vapply(states, function(s) {
    pinned <- as_map_list(
      make.simmap(ph_clade, x_ph, Q = Q, nsim = 10,
                  pi = setNames(as.numeric(states == s), states), message = FALSE), 10
    )
    root_edge <- which(pinned[[1]]$edge[, 1] == Ntip(ph_clade) + 1L)[1]
    all(vapply(pinned, \(m) names(m$maps[[root_edge]])[1], character(1)) == s)
  }, logical(1)))

  stopifnot(
    "make.simmap must sample, not fix, a tip with a uniform prior" = varies,
    "make.simmap must honour a point-mass pi at the root for every state" = pinned_ok
  )
  message("check 6b: tip priors are multiplicative and point-mass pi is honoured for all k states.")
})

# 6c. With Q_bg == Q_ph the split must be invisible, so split-and-stitch node
# posteriors have to match a plain single-Q map. Multinomial cell proportions
# are marginally binomial, so the phoneme script's pooled two-proportion z
# applies unchanged — just evaluated over every node x state cell.
if (RUN_SPLIT_VALIDATION) local({
  probe   <- targets |> filter(k == max(k)) |> slice(1)
  feature <- probe$feature
  states  <- feature_states(feature)
  v       <- feature_values(feature)
  Q1      <- regime_q(feature, states)$bg

  split_post <- node_posterior(regime_simmap(v, states, Q1, Q1, nsim = NSIM)$maps, states)
  plain <- as_map_list(
    make.simmap(analysis_tree, prior_matrix(v, analysis_tree$tip.label, states),
                Q = Q1, nsim = NSIM, pi = ROOT_PRIOR, message = FALSE),
    NSIM
  )
  plain_post <- node_posterior(plain, states)

  d      <- abs(split_post - plain_post)
  p_pool <- (split_post * NSIM + plain_post * NSIM + 1) / (2 * NSIM + 2)
  z      <- d / sqrt(2 * p_pool * (1 - p_pool) / NSIM)

  message(
    "check 6c: split vs single-Q node posteriors on ", feature,
    " (k = ", length(states), ") at nsim = ", NSIM,
    " -- mean |delta| = ", round(mean(d), 4),
    ", max |delta| = ", round(max(d), 4), ", max |z| = ", round(max(z), 2),
    " (", sum(z > 3), "/", length(z), " cells over |z| = 3)."
  )
  stopifnot(
    "regime split must be invisible when Q_bg == Q_ph (systematic shift)" = mean(d) < 0.02,
    "regime split must be invisible when Q_bg == Q_ph (outlier cell)"     = max(z) < 5
  )
})

# 6d. Exact posterior marginal at the Philippine MRCA under the real
# Q_bg != Q_ph, versus what the background maps sampled. Computed from the
# pruning algorithm, so it carries no Monte Carlo error.
exact_root_marginal <- function(v, states, Q_bg, Q_ph) {
  k     <- length(states)
  x_ph  <- prior_matrix(v, ph_clade$tip.label, states)
  pl_ph <- partial_lik(ph_clade, x_ph, Q_ph)
  x_obs <- prior_matrix(v, setdiff(bg_tree$tip.label, PH_PLACEHOLDER), states)

  joint <- vapply(seq_along(states), function(s) {
    pinned <- rbind(
      x_obs,
      matrix(as.numeric(seq_along(states) == s), 1, k,
             dimnames = list(PH_PLACEHOLDER, states))
    )[bg_tree$tip.label, , drop = FALSE]
    pl_bg <- partial_lik(bg_tree, pinned, Q_bg)
    # log P(clade data | MRCA = s) + log P(background data, MRCA = s)
    pl_ph$logscale + log(pl_ph$lik[s]) + pl_bg$logscale + log(sum(pl_bg$lik / k))
  }, numeric(1))

  p <- exp(joint - max(joint))
  setNames(p / sum(p), states)
}

# ── 7. Map every target feature ──────────────────────────────

results <- vector("list", nrow(targets))
names(results) <- targets$feature
t_start <- Sys.time()

for (i in seq_len(nrow(targets))) {
  row     <- targets[i, ]
  feature <- row$feature
  states  <- feature_states(feature)
  v       <- feature_values(feature)
  Q       <- regime_q(feature, states)

  fit <- regime_simmap(v, states, Q$bg, Q$ph, nsim = NSIM)

  # check 6d, per feature: sampled root-state frequencies vs their exact
  # marginal. A k-cell table can hold several near-zero cells, so the chi-square
  # falls back to a simulated p-value rather than the phoneme script's clamp.
  exact    <- exact_root_marginal(v, states, Q$bg, Q$ph)
  observed <- tabulate(match(fit$root_states, states), length(states))
  p_exact  <- pmax(exact, 1e-9)
  p_exact  <- p_exact / sum(p_exact)
  gof <- suppressWarnings(chisq.test(
    observed, p = p_exact,
    simulate.p.value = any(p_exact * NSIM < 5), B = 2000
  ))

  results[[i]] <- list(
    feature    = feature,
    states     = states,
    posterior  = node_posterior(fit$maps, states),
    trans      = transition_counts(fit$maps, states),
    root_exact = exact,
    root_obs   = observed / NSIM,
    root_p     = gof$p.value,
    maps       = fit$maps
  )

  message(
    "  [", i, "/", nrow(targets), "] ", feature, " (", row$winning_model,
    ", k = ", length(states), ") -- MRCA exact [",
    paste(round(exact, 3), collapse = ", "), "] vs sampled [",
    paste(round(observed / NSIM, 3), collapse = ", "), "] (p = ",
    signif(gof$p.value, 2), "); ",
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
        nrow(targets), " features.")

# ── 8. Summary tables ────────────────────────────────────────

node_posteriors <- bind_rows(lapply(results, function(r) {
  as_tibble(r$posterior) |>
    mutate(node = internal_nodes, n_desc_tips = node_n_desc, in_philippine_clade = node_in_clade) |>
    pivot_longer(all_of(r$states), names_to = "state", values_to = "p_state") |>
    mutate(feature = r$feature, k = length(r$states))
})) |>
  left_join(dplyr::select(targets, feature, winning_model), by = "feature") |>
  dplyr::select(feature, k, winning_model, node, state, p_state,
                n_desc_tips, in_philippine_clade) |>
  arrange(feature, node, state)

transitions <- bind_rows(lapply(results, function(r) {
  as_tibble(r$trans) |>
    pivot_longer(everything(), names_to = "key", values_to = "n") |>
    separate(key, into = c("regime", "from_state", "to_state"), sep = "\\|") |>
    group_by(regime, from_state, to_state) |>
    summarise(
      mean   = mean(n),
      median = median(n),
      q025   = quantile(n, 0.025),
      q975   = quantile(n, 0.975),
      .groups = "drop"
    ) |>
    mutate(
      feature             = r$feature,
      k                   = length(r$states),
      total_branch_length = regime_branch_length[regime],
      per_unit_time       = mean / total_branch_length
    )
})) |>
  left_join(dplyr::select(targets, feature, winning_model), by = "feature") |>
  dplyr::select(feature, k, winning_model, regime, from_state, to_state,
                mean, median, q025, q975, total_branch_length, per_unit_time) |>
  arrange(feature, regime, from_state, to_state)

out_dir <- here("data", "evolutionary", "stochastic_character_map")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(node_posteriors, file.path(out_dir, "GRAMMAR_evolutionary_simmap_node_posteriors.csv"))
write_csv(transitions, file.path(out_dir, "GRAMMAR_evolutionary_simmap_transitions.csv"))
if (SAVE_MAPS) {
  saveRDS(
    lapply(results, \(r) r$maps),
    file.path(out_dir, "GRAMMAR_evolutionary_simmap_maps.rds")
  )
}

# ── 9. Figures ───────────────────────────────────────────────

fig_dir <- here("figures", "evolutionary")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

for (i in seq_len(nrow(targets))) {
  row     <- targets[i, ]
  r       <- results[[i]]
  states  <- r$states
  k       <- length(states)
  piecol  <- STATE_COLORS[seq_len(k)]

  v       <- feature_values(row$feature)
  obs_idx <- which(!is.na(v) & as.character(v) %in% states)
  obs_pie <- matrix(0, length(obs_idx), k)
  obs_pie[cbind(seq_along(obs_idx), match(as.character(v[obs_idx]), states))] <- 1

  png(
    file.path(fig_dir, paste0("simmap_", row$feature, "_", row$winning_model, ".png")),
    width = FIG_W, height = FIG_H, units = "in", res = FIG_DPI
  )

  plotSimmap(
    regime_tree,
    type   = "fan",
    colors = REGIME_COLORS,
    fsize  = FSIZE,
    ftype  = "reg",
    lwd    = 1.0,
    offset = TIP_OFFSET,
    mar    = c(2.2, 0.2, 2.2, 0.2)
  )
  nodelabels(pie = r$posterior, piecol = piecol, cex = PIE_CEX)
  # Tips whose value is missing get no marker: an all-ones prior drawn as an
  # even k-way pie would be indistinguishable from a genuinely uncertain node.
  tiplabels(pie = obs_pie, tip = obs_idx, piecol = piecol, cex = PIE_CEX)

  title(main = sprintf("%s — %s", row$feature, row$winning_model), cex.main = 1.1)
  # mtext, not title(sub =): title places a subtitle at line 4, outside a margin
  # this shallow, and clips it away silently.
  mtext(
    sprintf("states %s   |   %d tips missing",
            paste(states, collapse = "/"), sum(is.na(v))),
    side = 1, line = 0.7, cex = 0.8
  )

  dev.off()
}

# ── 10. Report ───────────────────────────────────────────────
# Transitions per unit branch length inside vs outside the Philippine clade.
# ph_over_bg will not equal the fitted rate ratios: a mapped count is rate x time
# occupying the originating state, and occupancy differs between regimes.

rate_check <- transitions |>
  dplyr::select(feature, winning_model, regime, from_state, to_state, per_unit_time) |>
  pivot_wider(names_from = regime, values_from = per_unit_time) |>
  mutate(ph_over_bg = Philippine / Background)

message(
  "\nDone in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min.",
  "\n  node posteriors : ", nrow(node_posteriors), " rows -> ",
  "data/evolutionary/stochastic_character_map/GRAMMAR_evolutionary_simmap_node_posteriors.csv",
  "\n  transitions     : ", nrow(transitions), " rows -> ",
  "data/evolutionary/stochastic_character_map/GRAMMAR_evolutionary_simmap_transitions.csv",
  "\n  cached maps     : ", if (SAVE_MAPS) "GRAMMAR_evolutionary_simmap_maps.rds" else "(skipped)",
  "\n  figures         : ", nrow(targets), " PNGs -> figures/evolutionary/\n"
)

print(as.data.frame(rate_check))
