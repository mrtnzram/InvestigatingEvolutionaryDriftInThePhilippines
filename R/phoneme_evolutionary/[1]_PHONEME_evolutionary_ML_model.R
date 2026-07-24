# =============================================================================
# [0] Phoneme Evolutionary Maximum Likelihood Models
#
# first pass test on runnning an ape() state-based markov model on one phoneme 
# 1. compare ER and ADR model
# 2. compare 1 regime vs 2 regime model
#
# Inputs:data/RUHLENdf_PH_evolutionary.csv, data/phoneme_evolutionary_tree.nwk
#         data/phoneme_evolutionary_tip_mapping.csv
#
# Outputs:
# =============================================================================

library(ape)
library(dplyr)
library(readr)
library(here)
library(phytools)

# ── 1. Load data ─────────────────────────────────────────────

RUHLENdf <- read_csv(
  here("data", "RUHLENdf_PH_evolutionary.csv"),
  show_col_types = FALSE
)

tree <- read.tree(
  here("data", "phoneme_evolutionary_tree.nwk")
)

tip_mapping <- read_csv(
  here("data", "phoneme_evolutionary_tip_mapping.csv"),
  show_col_types = FALSE
)

phoneme_col <- "phoneme_086"

# Attach the selected Ruhlen phoneme state to each tree tip
phoneme_data <- tip_mapping %>%
  left_join(
    RUHLENdf %>%
      select(
        iso6393,
        language,
        all_of(phoneme_col)
      ),
    by = c("iso6393", "ruhlen_language" = "language")
  ) %>%
  transmute(
    tip,
    iso6393,
    ruhlen_language,
    phoneme_state = as.integer(.data[[phoneme_col]])
  ) %>%
  filter(
    !is.na(phoneme_state),
    phoneme_state %in% c(0L, 1L)
  )


# ── Prune and align the tree and trait vector ─────────────

phoneme_tree_test <- keep.tip(
  tree,
  phoneme_data$tip
)

phoneme_vec <- phoneme_data$phoneme_state
names(phoneme_vec) <- phoneme_data$tip

# Reorder the trait vector to exactly match tree-tip order
phoneme_vec <- phoneme_vec[
  phoneme_tree_test$tip.label
]

# Treat 0 and 1 as discrete states
phoneme_vec <- factor(
  phoneme_vec,
  levels = c(0, 1)
)

names(phoneme_vec) <- phoneme_tree_test$tip.label


#rooting the tree
phoneme_tree_test$root.edge <- 0
is.rooted(phoneme_tree_test)


# fix polytomies : consult with dr. creanza
# Number of descendants from each internal node
node_degree <- table(phoneme_tree_test$edge[, 1])

# Nodes with >2 descendants = polytomies
node_degree[node_degree > 2]

phoneme_tree_test <- multi2di(
  phoneme_tree_test,
  random = FALSE
)

is.binary(phoneme_tree_test)
is.rooted(phoneme_tree_test)

# Fit model

fit_ER <- ace(
  x = phoneme_vec,
  phy = phoneme_tree_test,
  type = "discrete",
  model = "ER"
)

fit_ARD <- ace(
  x = phoneme_vec,
  phy = phoneme_tree_test,
  type = "discrete",
  model = "ARD"
)

anova(fit_ER, fit_ARD)

#ARD is not significantly better choose ER

fit_ER$rates
fit_ARD$rates


# --- 1 vs 2 regime test ------

philippine_tips <- tip_mapping %>%
  filter(
    ruhlen_language %in% Ph_Languages
  ) %>%
  pull(tip) %>%
  intersect(
    phoneme_tree_test$tip.label
  )

length(philippine_tips)

philippine_tip_numbers <- match(
  philippine_tips,
  phoneme_tree_test$tip.label
)

# grab clade
philippine_node <- getMRCA(
  phoneme_tree_test,
  ph_tip_names
)

#paint clade
regime_tree <- paintSubTree(
  tree = phoneme_tree_test,
  node = philippine_node,
  state = "Philippine",
  anc.state = "Background",
  stem = FALSE
)

plotSimmap(
  regime_tree,
  colors = c(
    Background = "grey70",
    Philippine = "red"
  ),
  ftype = "i",
  fsize=0.4
)

#fit models

fit_1regime <- fitMk(
  tree = phoneme_tree_test,
  x = phoneme_vec,
  model = "ER"
)
fit_2regime <- fitmultiMk(
  tree = regime_tree,
  x = phoneme_vec,
  model = "ER"
)

#which one is better?
data.frame(
  Model = c(
    "1-regime ER",
    "2-regime ER"
  ),
  LogLik = c(
    logLik_1regime,
    logLik_2regime
  ),
  Parameters = c(
    1,
    2
  )
)


logLik_1regime <- as.numeric(
  logLik(fit_1regime)
)

logLik_2regime <- fit_2regime$logLik

LR <- 2 * (
  logLik_2regime -
    logLik_1regime
)

p_value <- pchisq(
  LR,
  df = 1,
  lower.tail = FALSE
)

cat(
  "\nLikelihood-ratio test",
  "\nChi-square =", LR,
  "\ndf = 1",
  "\np =", p_value,
  "\n\nBackground rate =", q_background,
  "\nPhilippine rate =", q_philippine,
  "\n"
)


