# =============================================================================
# Phoneme Analysis — pure-phylogeny methods figure
#
# Four-panel methods figure tracing the [4] construction on the phoneme data:
# (A) phylogenetic Axis 1, (B) the Manila-anchored signed geographic axis against
# distance from Manila, (C) the geography fit of Axis 1, (D) the residual pure phylogeny.
#
# Input:   `df`, `pvr_geo_lm` from [4]_PHONEME_PVR.R; data/phoneme/PHONEME_subgroup_lookup.csv
# Outputs: figures/regression/phoneme_pure_phylogeny_steps.png
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggrepel)
library(here)

stopifnot(
  "Run [4]_PHONEME_PVR.R first: `pvr_geo_lm` is not defined." = exists("pvr_geo_lm"),
  "Run [4]_PHONEME_PVR.R first: `df` lacks the pure-phylogeny columns." =
    exists("df") && all(c("phylo_axis1", "geo_signed_km", "pure_phylogeny") %in% names(df)),
  "`df` and `pvr_geo_lm` disagree on row count." = nrow(df) == length(fitted(pvr_geo_lm))
)

FOCUS <- c(Tagalog = "#D62728", Kapampangan = "#FF7F0E", Ilokano = "#1F77B4",
           Cebuano = "#2CA02C", Sangil = "#9467BD")

sub <- read.csv(here("data", "phoneme", "PHONEME_subgroup_lookup.csv"))
d <- df |>
  left_join(sub |> dplyr::select(language, colour), by = "language") |>
  mutate(s = geo_signed_km / 100, h = x_km / 100, pred = fitted(pvr_geo_lm))
set.seed(2)
d$lane <- runif(nrow(d), -0.25, 0.25)                # panel A jitter, fixed so focus markers sit on their point
foc <- d |> filter(language %in% names(FOCUS)) |> distinct(language, .keep_all = TRUE) |>
  mutate(fc = FOCUS[language])

stopifnot("Focus language missing from the phoneme frame." = setequal(foc$language, names(FOCUS)))

th <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.tag = element_text(face = "bold", size = 22))
x_signed <- "Geographic Axis 1 (100 km, Manila = 0)"
manila_v <- geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35")
stat_box <- function(label, top = FALSE) {
  annotate("label", x = Inf, y = if (top) Inf else -Inf, hjust = 1.05,
           vjust = if (top) 1.3 else -0.3, size = 3.3, label = label)
}

# Focus languages keep one colour and marker in every panel; A labels them, B carries the legend.
focus <- function(xv, yv, ...) {
  list(geom_point(data = foc, aes(.data[[xv]], .data[[yv]]), colour = foc$fc, size = 3.2),
       geom_label_repel(data = foc, aes(.data[[xv]], .data[[yv]], label = language),
                        colour = foc$fc, size = 3.2, label.size = 0.2, min.segment.length = 0,
                        box.padding = 0.5, seed = 1, ...))
}
focus_pts <- function(xv, yv, legend = FALSE) {
  list(geom_point(data = foc, aes(.data[[xv]], .data[[yv]], colour = language), size = 3.2),
       scale_colour_manual(values = FOCUS, breaks = names(FOCUS), name = NULL,
                           guide = if (legend) "legend" else "none"))
}


# ── A. Phylogenetic Axis 1 (colour = subgroup) ────────────────────────────
p_a <- ggplot(d, aes(phylo_axis1, lane)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_point(aes(colour = colour), size = 2.2, alpha = 0.8) +
  scale_colour_identity() +
  annotate("text", x = min(d$phylo_axis1), y = -0.42, hjust = 0, size = 3.2,
           label = "Northern Luzon subgroups\n(Meso-Cordilleran, Cagayan Valley)") +
  annotate("text", x = max(d$phylo_axis1), y = -0.42, hjust = 1, size = 3.2,
           label = "Central & southern subgroups\n(Central Philippine, Palawanic, Manobo, ...)") +
  focus("phylo_axis1", "lane", nudge_y = 0.3, direction = "x") +
  scale_y_continuous(limits = c(-0.55, 0.65)) +
  th +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(x = "Phylogenetic Axis 1", y = NULL)


# ── B. Signed geographic axis vs. distance from Manila ────────────────────
u <- d |> distinct(language, .keep_all = TRUE)
p_b <- ggplot(u, aes(s, h)) +
  geom_abline(slope = c(-1, 1), intercept = 0, linetype = "dotted", colour = "grey50") +
  manila_v +
  geom_point(colour = "grey55", size = 2, alpha = 0.7) +
  annotate("point", x = 0, y = 0, shape = 8, size = 4, stroke = 1.2) +
  focus_pts("s", "h", legend = TRUE) +
  stat_box(sprintf("r(distance, |axis 1|) = %.2f", cor(u$h, abs(u$s)))) +
  th +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1), legend.key.size = unit(0.9, "lines"),
        legend.background = element_rect(colour = "black", linewidth = 0.3)) +
  labs(x = x_signed, y = "Distance from Manila (100 km)")


# ── C. Geography fit of Axis 1 (Axis 1 ~ signed axis + distance from Manila) ──
p_c <- ggplot(d, aes(s, phylo_axis1)) +
  manila_v +
  geom_segment(aes(xend = s, yend = pred), colour = "grey75", linewidth = 0.4) +
  geom_segment(data = foc, aes(xend = s, yend = pred), colour = foc$fc, linewidth = 1.1) +
  geom_line(data = arrange(d, s), aes(y = pred), colour = "black", linewidth = 1.1) +
  geom_point(aes(y = pred), colour = "black", size = 1.2) +
  geom_point(colour = "grey55", size = 2, alpha = 0.7) +
  focus_pts("s", "phylo_axis1") +
  stat_box(sprintf("r^2 = %.2f", summary(pvr_geo_lm)$r.squared)) +
  th +
  labs(x = x_signed, y = "Phylogenetic Axis 1 (oriented)")


# ── D. Pure phylogeny (residual of C) ─────────────────────────────────────
p_d <- ggplot(d, aes(s, pure_phylogeny)) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 1.1) +
  manila_v +
  geom_segment(aes(xend = s, yend = 0), colour = "grey75", linewidth = 0.4) +
  geom_segment(data = foc, aes(xend = s, yend = 0), colour = foc$fc, linewidth = 1.1) +
  geom_point(colour = "grey55", size = 2, alpha = 0.7) +
  focus_pts("s", "pure_phylogeny") +
  stat_box(sprintf("r^2 = %.2f", summary(lm(pure_phylogeny ~ s + h, data = d))$r.squared), top = TRUE) +
  th +
  labs(x = x_signed, y = "Pure phylogeny (residual)")


steps <- (p_a | p_b) / (p_c | p_d) + plot_annotation(tag_levels = "A")
print(steps)

ggsave(here("figures", "regression", "phoneme_pure_phylogeny_steps.png"),
       steps, width = 12, height = 9.5, units = "in", dpi = 300)
