library(ggplot2)
library(dplyr)
library(here)
library(scales)

phon_surface <- read.csv(here("data", "phoneme_surface_raster.csv"))
nodepos      <- read.csv(here("data", "nodepos_phoneme.csv"))
PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))

global_lim <- c(0, max(PHONEME_cossim$cossim_span, GRAMMAR_cossim$cossim_span, na.rm = TRUE))

vmax <- max(abs(phon_surface$log_w_ratio), na.rm = TRUE)
lims <- c(-vmax, vmax)

world_map  <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_tile(data = phon_surface, aes(x = lon, y = lat, fill = log_w_ratio), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = PHONEME_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 6, alpha = 0.7) +
  geom_point(data = PHONEME_cossim, aes(x = longitude, y = latitude),
             size = 6, shape = 21, color = "black") +
  scale_fill_gradientn(colors = c("orange", "white", "cyan"),
                       limits = range(phon_surface$log_w_ratio, na.rm = TRUE),
                       na.value = "transparent") +
  scale_color_gradient(low = "white", high = "navy", limits = global_lim) +
  guides(
    fill = guide_colorbar(title = expression(log[10](w/bar(w))), title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title = "Cosine Similarity", title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_x_continuous(breaks = seq(115, 130, by = 2)) +
  scale_y_continuous(breaks = seq(4, 22, by = 2)) +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    values = rescale(c(-vmax, 0, vmax), from = lims),  # forces white -> exactly 0
    limits = lims,
    na.value = "transparent"
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))

base_plot
saveRDS(base_plot, file = here("data", "base_plot_phoneme_FEEMS.rds"))
