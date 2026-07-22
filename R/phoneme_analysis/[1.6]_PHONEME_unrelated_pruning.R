# =============================================================================
# [1.6] Phoneme Analysis — Unrelated-language pruning candidates
# For each Philippine language, finds the right-skew cutoff (median + 2*sd) of
# its own distribution of cosine similarity to the 137 unrelated controls, then
# flags which unrelated languages exceed that cutoff. Tallying the flags across
# all 58 Philippine languages ranks the unrelated controls by how often they
# show up as an anomalous outlier — candidates for further pruning of the
# unrelated set (likely contact-contamination leaks that survived the
# contact_contaminated filter in [0]_CREANZA_RUHLENdatabase.R).
#
# This script only identifies and ranks candidates; it does not modify
# RUHLENdf_PH.csv, PHONEME_cosine_matrix.csv, or any downstream pipeline file.
#
# Inputs:  data/RUHLENdf_PH.csv, data/PHONEME_cosine_matrix.csv
# Output:  data/PHONEME_unrelated_skew_candidates.csv (137 rows: language, n_exceed)
# =============================================================================

library(tidyverse)
library(ggplot2)
library(here)

# --- Loading data -------------------------------------------------------------
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"))
ph_lang  <- RUHLENdf |> filter(Language_type == "Philippine Language") |> pull(language)
unr_lang <- RUHLENdf |> filter(Language_type == "Unrelated Language")  |> pull(language)

# Per-PH-language null: rows = Philippine languages, cols = unrelated controls
# (same slice as [2]_PHONEME_cosine_distribution_analysis.R's load_null_matrix).
load_null_matrix <- function(ph_lang, unr_lang) {
  read.csv(here("data", "PHONEME_cosine_matrix.csv"),
           row.names = 1, check.names = FALSE) |>
    as.matrix() |>
    (\(m) m[ph_lang, unr_lang, drop = FALSE])()
}
null_mat <- load_null_matrix(ph_lang, unr_lang)

# --- Per-PH-language cutoff (median + 2*sd) -----------------------------------
# Each row is one PH language's distribution of similarity to the 137 unrelated
# controls. The cutoff is row-specific: a PH language's own median/sd, not a
# pooled/global statistic.
cutoff_tbl <- tibble(
  language     = rownames(null_mat),
  row_median   = apply(null_mat, 1, median),
  row_sd       = apply(null_mat, 1, sd),
) |>
  mutate(cutoff = row_median + 2 * row_sd)

# Same-shape logical matrix: exceed[i, j] = null_mat[i, j] > cutoff_i (row-wise).
exceed <- null_mat > matrix(cutoff_tbl$cutoff, nrow = nrow(null_mat), ncol = ncol(null_mat))

# --- Visualize the cutoff: representative PH languages ------------------------
# Sample 6 PH languages spread across low/mid/high skew (tertiles of cutoff,
# similar sampling idea to [2]'s Shapiro-Wilk rep_langs), 2 per tertile.
rep_langs <- cutoff_tbl |>
  mutate(tertile = ntile(cutoff, 3) |>
           factor(labels = c("low cutoff", "mid cutoff", "high cutoff"))) |>
  arrange(tertile, cutoff) |>
  group_by(tertile) |>
  slice(round(quantile(seq_len(n()), probs = c(0.33, 0.67)))) |>
  ungroup()

null_sel_long <- as_tibble(null_mat[rep_langs$language, , drop = FALSE], rownames = "language") |>
  pivot_longer(-language, names_to = "unr_language", values_to = "similarity") |>
  left_join(rep_langs |> select(language, tertile, cutoff), by = "language") |>
  mutate(language = factor(language, levels = rep_langs$language))

cutoff_example_plot <- ggplot(null_sel_long, aes(x = similarity)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(data = rep_langs |> mutate(language = factor(language, levels = rep_langs$language)),
             aes(xintercept = cutoff),
             linetype = "dashed", color = "#e07a5f", linewidth = 1) +
  facet_wrap(~ tertile + language, scales = "free", ncol = 2,
             labeller = label_wrap_gen(multi_line = FALSE)) +
  labs(title = "Unrelated-control similarity per Philippine language",
       subtitle = "Dashed line = that language's own cutoff (median + 2×sd)",
       x = "Cosine similarity to unrelated controls", y = "Count") +
  theme_bw()

cutoff_example_plot

ggsave(
  filename = here("figures", "phoneme", "distributions", "phoneme_unrelated_skew_cutoff_examples.png"),
  plot = cutoff_example_plot,
  width = 8, height = 8, units = "in", dpi = 300
)

# --- Ranked candidate tibble ---------------------------------------------------
# Tally exceed[] by column (by unrelated language): how many of the 58 PH-
# language distributions each unrelated language exceeded the cutoff in.
skew_candidates <- tibble(
  language  = colnames(null_mat),
  n_exceed  = colSums(exceed)
) |>
  arrange(desc(n_exceed))

print(skew_candidates, n = Inf)

stopifnot(nrow(skew_candidates) == length(unr_lang))

write.csv(skew_candidates, file = here("data", "PHONEME_unrelated_skew_candidates.csv"), row.names = FALSE)
