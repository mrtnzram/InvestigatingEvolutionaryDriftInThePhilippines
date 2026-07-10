library(caper)
library(tidyverse)
library(here)

PHONEME_PGLS <- PHONEME_cossim |> 
  select(language,cossim_span,geodist_H1_span) |> 
  left_join(
    select(tree_df_matched, original, ph),
    by = c("language" = "ph")
  )

PHONEME_PGLS |> 
  select(original, cossim_span, distance = geodist_H1_span)

PHONEME_PGLS_df <- PHONEME_PGLS |> 
  as.data.frame()

# 'df' must have a column matching tip labels in the tree
comp_data <- comparative.data(
  phy = tree_pruned,
  data = PHONEME_PGLS_df,
  names.col = "original",
  vcv = TRUE
)

# Estimate lambda by ML
pgls_fit <- pgls(
  cossim_span ~ geodist_H1_span,
  data = comp_data,
  lambda = "ML" # estimates lambda from data
)
summary(pgls_fit)
# Key outputs: slope on distance_z, estimated lambda, AIC

pgls_ols <- pgls(cossim_span ~ geodist_H1_span, data = comp_data, lambda = 1e-6)
pgls_ml  <- pgls_fit  # lambda = "ML", already fit

AIC(pgls_ols, pgls_ml)
