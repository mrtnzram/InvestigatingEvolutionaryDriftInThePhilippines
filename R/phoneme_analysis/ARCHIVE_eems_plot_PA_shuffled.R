# =============================================================================
# Phoneme Analysis — Averaged shuffled (null-model) EEMS surface
# Reads the 100 shuffled EEMS MCMC runs, averages the log-migration rate per
# location, and plots the averaged null surface for comparison with the real one.
# =============================================================================

library("reemsplots2")
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(purrr)
library(here)

PHONEME_cossim <- read.csv(here("data", "PHONEME_cossim.csv"))

# ---- 1. Process each of the 100 shuffled simulations into log-rate point sets ----
base_path <- here("matlab", "matlab_eems_phonemes_shuffled")

all_sim_data <- vector("list", 100)

for (i in 1:100) {
  sim_path <- file.path(base_path, paste0("datapath-g6x4-simno", i))

  # Process each simulation safely (skip failures)
  try({
    plots <- make_eems_plots(sim_path, longlat = TRUE)

    data_log <- plots$mrates02$data %>%
      mutate(z_log = ifelse(z == 0, 0,log(z)))

    data_pts <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

    coords <- st_coordinates(data_pts)
    data_pts$x <- coords[, "Y"]
    data_pts$y <- coords[, "X"]

    all_sim_data[[i]] <- data_pts
  }, silent = FALSE)
}

# ---- 2. Average the log-migration rate across simulations per location ----
valid_data <- compact(all_sim_data)

master_data <- bind_rows(valid_data)

averaged_data <- master_data %>%
  group_by(x, y) %>%
  summarise(z_log_mean = mean(z_log, na.rm = TRUE), .groups = "drop")

averaged_sf <- st_as_sf(averaged_data, coords = c("y", "x"), crs = 4326)

# ---- 3. Plot the averaged null migration surface ----
ph_shape <- ne_countries(scale = "medium", country = "Philippines", returnclass = "sf")

data_log <- plots$mrates02$data %>%
  mutate(z_log = log(z))

data_pts <- st_as_sf(data_log, coords = c("y", "x"), crs = 4326)

coords <- st_coordinates(data_pts)
data_pts$x <- coords[, "Y"]
data_pts$y <- coords[, "X"]

data_clean <- averaged_data %>%
  filter(is.finite(z_log_mean))

ggplot() +
  geom_tile(data = averaged_data, aes(x = y, y = x, fill = z_log_mean)) +
  geom_sf(data = ph_shape, fill = NA, color = "black") +
  scale_fill_gradientn(
    colors = c("orange", "white", "cyan"),
    limits = range(averaged_data$z_log_mean, na.rm = TRUE),
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
