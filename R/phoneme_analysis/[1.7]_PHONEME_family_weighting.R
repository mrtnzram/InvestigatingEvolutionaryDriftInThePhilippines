# =============================================================================
# [1.7] Phoneme Analysis — Genealogical family weighting of the unrelated null
# The 137 unrelated controls are unevenly spread across 26 language families
# (Niger-Congo 19, Afro-Asiatic 17, Sino-Tibetan 15, ... 8 singletons), so a few
# clades dominate the "unrelated baseline" that Philippine-language similarity is
# compared against. This down-weights that imbalance: each control gets a
# frequency weight w = 1/family_size, so every family contributes equal mass to
# the unrelated null's mean and density — without altering any cosine value.
#
# This is a standalone diagnostic of how family weighting reshapes the unrelated
# null distribution + mean (and [1.6]'s skew cutoff). It does NOT feed the
# downstream pipeline: [2]-[5] use the raw [1] baselines.
#
# Inputs:  data/RUHLENdf_PH.csv, data/PHONEME_cosine_matrix.csv,
#          data/PHONEME_unrelated_families.csv, data/PHONEME_cossim.csv
# Outputs: data/PHONEME_family_weight_comparison.csv (per unrelated language:
#            n_exceed_raw, n_exceed_family, family, size)
#          figures/phoneme/corrections/phoneme_unrelated_family_distribution.png
#          figures/phoneme/corrections/phoneme_ridge_raw_vs_corrected.png
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggridges)
library(here)

fig_dir <- here("figures", "phoneme", "corrections")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Loading data --------------------------------------------------------------
RUHLENdf <- read_csv(here("data", "RUHLENdf_PH.csv"))
ph_lang  <- RUHLENdf |> filter(Language_type == "Philippine Language") |> pull(language)
unr_lang <- RUHLENdf |> filter(Language_type == "Unrelated Language")  |> pull(language)

phoneme_cols <- RUHLENdf |>
  dplyr::select(-language, -source, -iso6393, -Language_type, -latitude, -longitude) |>
  colnames()
inventory_size <- RUHLENdf |>
  transmute(language, size = rowSums(dplyr::across(dplyr::all_of(phoneme_cols)))) |>
  deframe()

cosine_matrix <- read.csv(here("data", "PHONEME_cosine_matrix.csv"),
                          row.names = 1, check.names = FALSE) |>
  as.matrix()
null_mat <- cosine_matrix[ph_lang, unr_lang, drop = FALSE]  # 58 x 137

# Family lookup + weights: w = 1/family_size, aligned to the 137 unrelated cols.
fam_lookup <- read_csv(here("data", "PHONEME_unrelated_families.csv")) |>
  add_count(language_family, name = "fam_n") |>
  mutate(w = 1 / fam_n)
fam_w <- setNames(fam_lookup$w, fam_lookup$language)[unr_lang]
stopifnot(!anyNA(fam_w))

# Weighted helpers (frequency weights; values are never altered).
weighted_median <- function(x, w) { o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) / sum(w) >= 0.5)[1]] }
weighted_sd <- function(x, w) { w <- w / sum(w); mu <- sum(w * x)
  sqrt(sum(w * (x - mu)^2) / (1 - sum(w^2))) }
wmean <- function(x, w) sum((w / sum(w)) * x)

# =============================================================================
# 1. Genealogical distribution of the unrelated set
# =============================================================================
fam_counts <- fam_lookup |> count(language_family, sort = TRUE)
cat(sprintf("Unrelated controls: %d languages across %d families.\n",
            nrow(fam_lookup), nrow(fam_counts)))
cat(sprintf("Top-3 families = %d / %d (%.0f%%); %d singleton families.\n",
            sum(sort(fam_counts$n, decreasing = TRUE)[1:3]), nrow(fam_lookup),
            100 * sum(sort(fam_counts$n, decreasing = TRUE)[1:3]) / nrow(fam_lookup),
            sum(fam_counts$n == 1)))

fam_bar <- ggplot(fam_counts |> mutate(language_family = fct_reorder(language_family, n)),
                  aes(x = n, y = language_family)) +
  geom_col(fill = "#2ca6a4", alpha = 0.85) +
  geom_text(aes(label = n), hjust = -0.3, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "Genealogical distribution of the 137 unrelated controls",
       subtitle = "Ruhlen language_family; the null baseline is dominated by a few large families",
       x = "Number of languages", y = NULL) +
  theme_minimal()
ggsave(file.path(fig_dir, "phoneme_unrelated_family_distribution.png"),
       fam_bar, width = 8, height = 6, units = "in", dpi = 300)

# =============================================================================
# 2. Family weighting adjustment to the unrelated null (mean + cutoff)
# =============================================================================
# n_exceed cutoff (median + 2*sd per PH-language row), unweighted vs family-
# weighted over the 137 unrelated columns. Count per unrelated language.
tally_skew <- function(m, w = NULL) {
  cutoff <- if (is.null(w)) apply(m, 1, \(r) median(r) + 2 * sd(r))
            else            apply(m, 1, \(r) weighted_median(r, w) + 2 * weighted_sd(r, w))
  exceed <- m > matrix(cutoff, nrow = nrow(m), ncol = ncol(m))
  tibble(language = colnames(m), n_exceed = colSums(exceed))
}

fw_comparison <- tally_skew(null_mat) |> rename(n_exceed_raw = n_exceed) |>
  left_join(tally_skew(null_mat, fam_w) |> rename(n_exceed_family = n_exceed), by = "language") |>
  mutate(family = fam_lookup$language_family[match(language, fam_lookup$language)],
         size   = inventory_size[language],
         delta  = n_exceed_family - n_exceed_raw) |>
  arrange(desc(n_exceed_raw))

print(fw_comparison, n = Inf)
write.csv(fw_comparison, file = here("data", "PHONEME_family_weight_comparison.csv"), row.names = FALSE)

# Unrelated null mean per PH language: raw rowMeans vs family-weighted mean.
cat(sprintf("\nUnrelated null mean (over 58 PH languages): raw = %.4f, family-weighted = %.4f\n",
            mean(rowMeans(null_mat)),
            mean(apply(null_mat, 1, \(r) wmean(r, fam_w)))))

# =============================================================================
# 3. Ridge comparison: raw vs family-weighted baseline distributions
# =============================================================================
# Only the Unrelated baseline changes under family weighting (interest languages
# are single languages, not weighted). cossim_unr(raw) = plain row mean;
# cossim_unr(family) = family-weighted row mean. Interest baselines identical.
PHONEME_cossim <- read_csv(here("data", "PHONEME_cossim.csv")) |>
  dplyr::select(-any_of("...1"))

bl <- PHONEME_cossim |>
  transmute(language,
            Unrelated_raw       = cossim_unr,
            Unrelated_corrected = apply(null_mat[language, , drop = FALSE], 1, \(r) wmean(r, fam_w)),
            Spanish_raw = cossim_span,  Spanish_corrected  = cossim_span,
            Japanese_raw = cossim_jap,  Japanese_corrected = cossim_jap,
            English_raw = cossim_eng,   English_corrected  = cossim_eng)

means_tbl <- bl |> summarise(across(-language, mean)) |>
  pivot_longer(everything(), names_to = c("baseline", "version"), names_sep = "_",
               values_to = "mean_cossim") |>
  pivot_wider(names_from = version, values_from = mean_cossim) |>
  mutate(shift = corrected - raw)
cat("\nBaseline means (over 58 PH languages), raw vs family-weighted:\n")
print(means_tbl)

rl <- bl |>
  pivot_longer(-language, names_to = c("baseline", "version"), names_sep = "_",
               values_to = "similarity") |>
  mutate(baseline = factor(baseline, levels = c("Unrelated", "English", "Spanish", "Japanese")),
         version  = factor(version,  levels = c("raw", "corrected")))
rm_means <- rl |> summarise(m = mean(similarity), .by = c(baseline, version))

ridge_compare <- ggplot(rl, aes(x = similarity, y = baseline, fill = version)) +
  geom_density_ridges(alpha = 0.45, scale = 1.1, color = "grey30", rel_min_height = 0.01) +
  geom_segment(data = rm_means,
               aes(x = m, xend = m, y = as.numeric(baseline), yend = as.numeric(baseline) + 0.9,
                   color = version),
               linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE) +
  scale_fill_manual(values = c(raw = "grey70", corrected = "#2ca6a4"),
                    labels = c(raw = "raw", corrected = "family-weighted")) +
  scale_color_manual(values = c(raw = "grey40", corrected = "#1b6f6d"),
                     labels = c(raw = "raw", corrected = "family-weighted")) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05)) +
  scale_y_discrete(expand = c(0.01, 0)) +
  labs(title = "Baseline cosine distributions: raw vs. family-weighted unrelated null",
       subtitle = "Family weighting (w = 1/family size) affects the Unrelated baseline only",
       x = "Cosine similarity", y = NULL, fill = NULL, color = NULL) +
  theme_minimal()
ggsave(file.path(fig_dir, "phoneme_ridge_raw_vs_corrected.png"),
       ridge_compare, width = 8, height = 5, units = "in", dpi = 300)
