# =============================================================================
# [3] Phoneme Analysis — Overlay waypoint routes on the EEMS base map
# Combines the EEMS migration surface (base_plot_PA) with the historical
# migration-route arrows and marks the colonial capital (Manila).
# =============================================================================

library(ggplot2)
library(here)

# ---- Load the EEMS base map and the waypoint-route arrow layers ----
base_plot <- readRDS(here("data", "base_plot_PA.rds"))
arrows    <- readRDS(here("data", "phoneme_waypoint_plot.rds"))

# ---- Overlay the arrow layers onto the cropped base map ----
final_plot_PA <- base_plot + coord_sf(xlim = c(116, 127), ylim = c(4, 21))
for (layer in arrows$layers) {
  final_plot_PA <- final_plot_PA + layer
}

# ---- Mark Manila (colonial capital) and render the combined figure ----
ref_coords1 <- c(121, 14.6)
capital_df <- data.frame(x = ref_coords1[1], y = ref_coords1[2], label = "Capital")

final_plot_PA_ <- final_plot_PA +
  geom_point(data = capital_df, aes(x = x, y = y, shape = label), color = "red", size = 4) +
  scale_shape_manual(values = c("Capital" = 18)) +
  guides(shape = guide_legend(title = NULL)) +
  theme(legend.position = "right")

print(final_plot_PA_)
