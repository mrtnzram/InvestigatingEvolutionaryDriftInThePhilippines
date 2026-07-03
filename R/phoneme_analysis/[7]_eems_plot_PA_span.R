# =============================================================================
# [7] Phoneme Analysis — Render the EEMS migration surface as the base map
# Reads the EEMS MCMC output, builds the log-migration-rate surface, overlays the
# per-language cosine-similarity points, and saves the base map for later overlays.
# Requires GRAMMAR_cossim in the environment (shared color scale across analyses).
# Output: data/base_plot_PA.rds
# =============================================================================

library("reemsplots2")
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)
library(here)

# ---- 1. Load EEMS MCMC output and the similarity scores ----
mcmcpath <- here("matlab", "matlab_eems_phonemes", "datapath-g6x4-simno1")
PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))

# ---- 2. Generate the standard EEMS diagnostic plots ----
plots <- make_eems_plots(mcmcpath, longlat = TRUE)

names(plots)

print(plots$mrates01)
plots$mrates02
plots$qrates02
plots$rdist01
plots$pilogl01

str(plots$mrates02$data)
summary(plots$mrates02$data)

# ---- 3. Quick look: log-migration surface over the Philippines ----
ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")

data_log <- plots$mrates02$data %>%
  mutate(z_log = log(z))

data_pts <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

coords <- st_coordinates(data_pts)
data_pts$x <- coords[, "Y"]
data_pts$y <- coords[, "X"]

ggplot() +
  geom_tile(data = data_pts, aes(x = y, y = x, fill = z_log)) +
  geom_sf(data = ph_shape, fill = NA, color = "black") +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    limits = range(data_pts$z_log, na.rm = TRUE),
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(116, 127), ylim = c(4, 21), expand = FALSE) +
  labs(
    fill = "log(Migration Rate)",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# ---- 4. Base map: migration surface + similarity points (shared color scale) ----
all_vals <- c(PHONEME_cossim$cossim_span,
              GRAMMAR_cossim$cossim_span)

global_lim <- range(all_vals, na.rm = TRUE)

ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")

data_log <- plots$mrates02$data %>%
  mutate(z_log = log(z))

data_pts <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

coords <- st_coordinates(data_pts)
data_pts$x <- coords[, "Y"]
data_pts$y <- coords[, "X"]

world_map <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_tile(data = data_pts, aes(x = y, y = x, fill = z_log), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = PHONEME_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 10, alpha = 0.7) +
  geom_point(data = PHONEME_cossim, aes(x = longitude, y = latitude),
             size = 10, shape = 21, color = "black") +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    limits = range(data_pts$z_log, na.rm = TRUE),
    na.value = "transparent"
  ) +
  scale_color_gradient(low = "white",high = "navy",
                       limits = global_lim) +
  guides(
    fill = guide_colorbar(title = "log(Migration Rate)", title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title = "Cosine Similarity", title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_x_continuous(breaks = seq(115, 130, by = 2)) +
  scale_y_continuous(breaks = seq(4, 22, by = 2)) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

base_plot

saveRDS(base_plot, file = here("data", "base_plot_PA.rds"))
