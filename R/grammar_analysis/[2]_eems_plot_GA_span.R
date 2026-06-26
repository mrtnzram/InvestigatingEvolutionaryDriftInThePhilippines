# =============================================================================
# [2] Grammar Analysis — Render the EEMS migration surface as the base map
# Reads the EEMS MCMC output, builds the (land-masked) log-migration-rate surface,
# overlays per-language cosine-similarity points, and saves the base map. Also
# plots the diversity-rate surface. Requires GRAMFEATURE_match_df and global_lim
# in the environment.
# Output: data/base_plot_GA.rds
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
mcmcpath <- here("matlab", "matlab_eems_grammar", "datapath-g8x8-simno1")
GRAMMAR_cossim <- read.csv(here("data", "GRAMMAR_cossim.csv"))

# ---- 2. Generate the standard EEMS diagnostic plots ----
plots <- make_eems_plots(mcmcpath, longlat = TRUE)

names(plots)

plots$mrates01
plots$qrates02
plots$rdist01
plots$pilogl01

str(plots$mrates02$data)
summary(plots$mrates02$data)

# ---- 3. Quick look: land-masked log-migration surface ----
ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")
ph_union <- st_union(ph_shape)

data_log <- plots$mrates02$data %>%
  mutate(z_log = if_else(z == 0, 0, log(z)))

data_pts <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

data_masked <- data_pts[st_intersects(data_pts, ph_union, sparse = FALSE), ]

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

# ---- 4. Base map: masked migration surface + similarity points ----
data_log <- plots$mrates02$data %>%
  mutate(z_log = if_else(z == 0, 0, log(z)))

tile_sf <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")
ph_union <- st_union(ph_shape)

tiles_masked <- tile_sf[st_intersects(tile_sf, ph_union, sparse = FALSE), ]

coords <- st_coordinates(tiles_masked)
tiles_masked$x <- coords[, "Y"]
tiles_masked$y <- coords[, "X"]

world_map <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

base_plot <- ggplot() +
  geom_tile(data = data_pts, aes(x = y, y = x, fill = z_log), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = GRAMMAR_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 10, alpha = 0.7) +
  geom_point(data = GRAMMAR_cossim, aes(x = longitude, y = latitude),
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

saveRDS(base_plot, file = here("data", "base_plot_GA.rds"))

# ---- 5. Diversity-rate surface with match-score points and top labels ----
data_log <- plots$qrates02$data %>%
  mutate(z_log = log(z))

tile_sf <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")
ph_union <- st_union(ph_shape)

tiles_masked <- tile_sf[st_intersects(tile_sf, ph_union, sparse = FALSE), ]

coords <- st_coordinates(tiles_masked)
tiles_masked$x <- coords[, "Y"]
tiles_masked$y <- coords[, "X"]

world_map <- map_data("world")
map_subset <- world_map %>% filter(region %in% c("Philippines", "Malaysia"))

label_df <- GRAMMAR_cossim %>%
  arrange(desc(cossim_span)) %>%
  distinct(language, .keep_all = TRUE) %>%
  slice_head(n = 10)

ggplot() +
  geom_tile(data = tiles_masked, aes(x = y, y = x, fill = z_log), alpha = 0.6) +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = NA, color = "black") +
  geom_point(data = GRAMMAR_cossim,
             aes(x = longitude, y = latitude, color = cossim_span),
             size = 10, alpha = 0.7) +
  geom_point(data = label_df, aes(x = longitude, y = latitude),
             size = 10, shape = 21, color = "black") +
  geom_text(data = label_df,
            aes(x = longitude, y = latitude, label = language),
            size = 2.5, fontface = "italic", color = "black", check_overlap = TRUE,
            nudge_y = 0.6) +
  scale_fill_gradientn(
    colors = c("orange", "white", "#0096FF"),
    limits = range(tiles_masked$z_log, na.rm = TRUE),
    na.value = "transparent"
  ) +
  scale_color_gradient(low = "white", high = "navy") +
  guides(
    fill = guide_colorbar(title = "log(Diversity Rate)", title.position = "top", title.hjust = 0.5),
    color = guide_colorbar(title = "Cosine Similarity (Spanish)", title.position = "top", title.hjust = 0.5)
  ) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  scale_x_continuous(breaks = seq(115, 130, by = 2)) +
  scale_y_continuous(breaks = seq(4, 22, by = 2)) +
  labs(
    title = "Cosine Similarity Over Diversity Surface",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )
