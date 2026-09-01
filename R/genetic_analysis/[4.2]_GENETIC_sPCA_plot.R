# =============================================================================
# [4.2] Genetic Analysis — sPCA plot (LOCAL companion to the remote fit)
#
# [4.2]_GENETIC_sPCA_remote.R runs headless on the server and writes its
# results as CSVs (see README_GENETIC_sPCA_remote.md's "Files to bring back");
# this reads them back and builds the same point-symbol map the other three
# domains' [4.2] scripts produce inline, since a map needs no server access.
#
# Input:   data/genetic/PVR/GENETIC_sPCA_scores.csv, GENETIC_sPCA_results.csv
#          (both scp'd back from the remote run)
# Outputs: figures/genetic/regression/genetic_sPCA_surface.png,
#          data/genetic/PVR/base_plot_genetic_sPCA.rds
# =============================================================================

library(tidyverse)
library(here)
library(maps)

scores_df <- read.csv(here("data", "genetic", "PVR", "GENETIC_sPCA_scores.csv"))
results   <- read.csv(here("data", "genetic", "PVR", "GENETIC_sPCA_results.csv"))


# ── Plot: sPCA point-symbol map (size = |sPC1|, colour = sign) ──────────────
dir.create(here("figures", "genetic", "regression"), recursive = TRUE, showWarnings = FALSE)

world_map  <- map_data("world")
map_subset <- world_map |> filter(region %in% c("Philippines", "Malaysia"))

mag_threshold <- median(abs(scores_df$sPC1))
scores_df <- scores_df |>
  mutate(
    sign      = if_else(sPC1 >= 0, "Positive", "Negative"),
    magnitude = if_else(abs(sPC1) >= mag_threshold, "Large", "Small"),
    key       = factor(paste(magnitude, sign),
                       levels = c("Large Positive", "Small Positive",
                                  "Small Negative", "Large Negative"))
  )

key_sizes  <- c("Large Positive" = 10, "Small Positive" = 4,
                "Small Negative" = 4,  "Large Negative" = 10)
key_fills  <- c("Large Positive" = "black", "Small Positive" = "black",
                "Small Negative" = "white", "Large Negative" = "white")
key_labels <- c("Large Positive" = "Large (+)", "Small Positive" = "Small (+)",
                "Small Negative" = "Small (-)", "Large Negative" = "Large (-)")

p_surface <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray97", color = "black") +
  geom_point(data = scores_df,
             aes(x = longitude, y = latitude, size = key, fill = key),
             shape = 22, colour = "black", stroke = 0.6) +
  scale_size_manual(values = key_sizes, labels = key_labels, name = "sPC1", drop = FALSE) +
  scale_fill_manual(values = key_fills, labels = key_labels, name = "sPC1", drop = FALSE) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Genetic sPCA",
       subtitle = sprintf("p = %.3f, r^2 = %.3f", results$perm_p, results$variance_explained))
print(p_surface)

ggsave(here("figures", "genetic", "regression", "genetic_sPCA_surface.png"),
       p_surface, width = 7.5, height = 6, units = "in", dpi = 300)
saveRDS(p_surface, file = here("data", "genetic", "PVR", "base_plot_genetic_sPCA.rds"))
