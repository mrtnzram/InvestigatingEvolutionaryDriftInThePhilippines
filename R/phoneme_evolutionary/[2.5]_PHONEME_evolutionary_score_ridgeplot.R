# =============================================================================
# [2.5] Phoneme Evolutionary — Ridgeplot of the donor model scores
#
# Plots the leave-one-out predictive scores from [2] across the 46 Philippine
# tips, one ridge per model.
#
# Two panels, because one alone misleads:
#   A. Raw scores. Between-language SD is ~0.108, about 17x the largest
#      between-model difference (~0.006), so the five ridges sit almost exactly
#      on top of each other. That is the honest picture of the effect size.
#   B. Paired differences vs the baseline (same language, same phoneme set), which
#      is where the model differences are actually visible. Pairing removes the
#      between-language variance that swamps panel A.
#
# The unrelated set is a NULL CONTROL, not a donor: its D is a phoneme-commonness
# proxy (r = +0.93 with austronesian_prevalence). It is drawn with the donors so
# the real donors can be read against it — if it keeps pace, the signal is
# commonness rather than contact.
#
# Inputs:  data/PHONEME_evolutionary_donor_scores.csv
# Outputs: figures/phoneme/distributions/donor_score_ridges.png
# =============================================================================

library(dplyr)
library(tidyr)
library(readr)
library(here)
library(ggplot2)
library(ggridges)
library(patchwork)

scores <- read_csv(
  here("data", "PHONEME_evolutionary_donor_scores.csv"),
  show_col_types = FALSE
)

# Ridge order is bottom-to-top; baseline first so the donors stack above it.
model_levels <- c("austronesian", "span", "eng", "jap", "unrelated")
model_labels <- c(
  austronesian = "Austronesian\n(baseline)",
  span         = "Spanish",
  eng          = "English",
  jap          = "Japanese",
  unrelated    = "Unrelated set\n(null control)"
)

# Categorical hues from the validated reference palette (slots 1-4), with a
# neutral for the baseline since it is a reference rather than a peer. Adjacent
# pairs in this stacking order clear both hard gates: worst normal-vision
# dE 17.8 (>=15), worst CVD dE 16.7 (>=8). Two slots fall under 3:1 contrast, so
# the relief rule applies — every ridge is direct-labelled on the y axis, and
# identity is never carried by colour alone.
model_colors <- c(
  austronesian = "#6b6a65",
  span         = "#2a78d6",
  eng          = "#008300",
  jap          = "#e87ba4",
  unrelated    = "#eda100"
)

ink_primary   <- "#0b0b0b"
ink_secondary <- "#52514e"
surface       <- "#fcfcfb"

long_scores <- scores |>
  select(tip, all_of(paste0(model_levels, "_score"))) |>
  pivot_longer(-tip, names_to = "model", values_to = "score") |>
  mutate(
    model = factor(sub("_score$", "", model), levels = model_levels)
  )

long_deltas <- long_scores |>
  left_join(
    scores |> select(tip, base = austronesian_score),
    by = "tip"
  ) |>
  filter(model != "austronesian") |>
  mutate(
    delta = score - base,
    model = factor(as.character(model), levels = setdiff(model_levels, "austronesian"))
  )

base_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.background   = element_rect(fill = surface, colour = NA),
    panel.background  = element_rect(fill = surface, colour = NA),
    panel.grid.minor  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "#e6e5e0", linewidth = 0.3),
    axis.text         = element_text(colour = ink_secondary),
    axis.text.y       = element_text(colour = ink_primary, hjust = 1),
    axis.title        = element_text(colour = ink_secondary),
    plot.title        = element_text(colour = ink_primary, face = "bold", size = 12),
    plot.subtitle     = element_text(colour = ink_secondary, size = 9.5, lineheight = 1.15),
    legend.position   = "none"   # every ridge is direct-labelled on the y axis
  )

sd_between <- sd(scores$austronesian_score)
max_delta  <- max(abs(long_deltas$delta))

p_raw <- ggplot(long_scores, aes(x = score, y = model, fill = model)) +
  geom_density_ridges(
    scale = 1.15, rel_min_height = 0.01,
    colour = surface, linewidth = 0.6, alpha = 0.9,
    jittered_points = TRUE, point_size = 0.7, point_alpha = 0.45,
    point_colour = ink_primary,
    position = position_points_jitter(width = 0, height = 0)
  ) +
  scale_fill_manual(values = model_colors) +
  scale_y_discrete(labels = model_labels, expand = expansion(add = c(0.2, 1.3))) +
  labs(
    title    = "A. Predictive score by model",
    subtitle = sprintf(
      "Balanced leave-one-out score per Philippine language (n = 46); 1 = perfect, 0.5 = chance.\nRidges nearly coincide: between-language SD (%.3f) is ~%.0fx the largest model difference (%.3f).",
      sd_between, sd_between / 0.0064, 0.0064
    ),
    x = "Score  (geometric mean probability assigned to the truth)", y = NULL
  ) +
  base_theme

p_delta <- ggplot(long_deltas, aes(x = delta, y = model, fill = model)) +
  geom_vline(xintercept = 0, colour = ink_secondary, linewidth = 0.4, linetype = "22") +
  geom_density_ridges(
    scale = 0.95, rel_min_height = 0.01,
    colour = surface, linewidth = 0.6, alpha = 0.9,
    jittered_points = TRUE, point_size = 0.7, point_alpha = 0.45,
    point_colour = ink_primary,
    position = position_points_jitter(width = 0, height = 0)
  ) +
  scale_fill_manual(values = model_colors) +
  scale_y_discrete(labels = model_labels, expand = expansion(add = c(0.2, 1.1))) +
  # coord_cartesian, not scale_x_continuous(limits=): the latter drops rows
  # outside the range BEFORE the density is computed, which would reshape the
  # ridges rather than just framing them.
  coord_cartesian(xlim = c(-max_delta * 0.55, max_delta * 1.35)) +
  labs(
    title    = "B. Difference from baseline, paired within language",
    subtitle = "Same language and same 84 phonemes on both sides, so between-language variance cancels.\nDashed line = no improvement over the Austronesian baseline.",
    x = "Score difference vs baseline", y = NULL
  ) +
  base_theme

ridge_plot <- p_raw / p_delta +
  plot_annotation(
    theme = theme(plot.background = element_rect(fill = surface, colour = NA))
  )

out_dir <- here("figures", "phoneme", "distributions")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(out_dir, "donor_score_ridges.png"),
  ridge_plot,
  width = 8, height = 8.5, dpi = 300, bg = surface
)

message(
  "Wrote figures/phoneme/distributions/donor_score_ridges.png\n",
  paste0(
    "  ", model_levels, ": mean ",
    sprintf("%.4f", tapply(long_scores$score, long_scores$model, mean)),
    collapse = "\n"
  )
)
