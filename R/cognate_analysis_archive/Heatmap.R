library(readr)
library(tidyverse)
library(ggplot2)
library(ggmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(maps)
library(akima)
library(patchwork)
library(rlang)
library(here)


PHONEME_cossim_marked_1 <- PHONEME_cossim_marked |> filter(span_influenced == TRUE)

plot_weighted_points_PA <- function(df, value_var, refdf = NULL, lims = NULL) {
  
  world_map <- map_data("world")
  map_subset <- world_map %>% 
    filter(region %in% c("Philippines", "Malaysia"))
  
  max_val <- max(df[[value_var]], na.rm = TRUE)
  min_val <- min(df[[value_var]], na.rm = TRUE)
  
  if (is.null(lims)) {
    lims <- c(min_val, max_val)
  }
  
  label_df <- df %>%
    arrange(desc(.data[[value_var]])) %>%
    distinct(language, .keep_all = TRUE) %>%
    slice_head(n = 10)
  
  p <- ggplot() +
    geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "black")
  
  if (!is.null(refdf)) {
    p <- p + geom_point(data = refdf, aes(x = longitude, y = latitude),
                        color = 'red', size = 2)
  }
  
  p <- p +
    geom_point(data = df,
               aes_string(x = "longitude", y = "latitude", color = value_var),
               size = 10,
               alpha = 0.7) +
    geom_point(data = df, aes(x = longitude, y = latitude),
               size = 10, shape = 21, color = "black") +
    #geom_text(data = label_df,
     #         aes(x = longitude, y = latitude, label = language),
      #        size = 2.5, fontface = "italic", color = "black", check_overlap = TRUE,
       #       nudge_y = 0.6) +
    scale_colour_gradient(
      low = "white",
      high = "navy",
      limits = lims,
      n.breaks = 10
    ) +
    guides(colour = guide_colourbar(
      barheight = unit(15, "cm"),
      barwidth  = unit(0.6, "cm"),
      title.position = "top",
      title.hjust = 0.5
    )) +
    coord_fixed(ratio = 1.3, xlim = c(115, 130), ylim = c(4, 22)) +
    scale_x_continuous(breaks = seq(115, 130, by = 2)) +
    scale_y_continuous(breaks = seq(4, 22, by = 2)) +
    theme_minimal() +
    labs(
      title = paste("Point Plot of", value_var, "over the Philippines"),
      x = "Longitude",
      y = "Latitude",
      color = value_var
    )
  
  return(p)
}


all_vals <- c(similarity_df_h1_span$similarity_score_Spanish,
              similarity_df_h2_span$similarity_score_Spanish)

global_lim <- range(all_vals, na.rm = TRUE)

similarity_df_h1_span <- similarity_df_h1_span %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))
similarity_df_h2_span <- similarity_df_h2_span %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

refdf1 <- data.frame(
  longitude = c(121),
  latitude = c(14.6)
)

refdf2 <- data.frame(
  long = c(125.75),
  lat = c(10.75)
)


PHONEME_cossim_shuffled <- read.csv(here("data", "PHONEME_cossim_shuffled.csv"))

PA_marked <- plot_weighted_points_PA(PHONEME_cossim_marked_1,'cossim_span',refdf1)

PA_p1_cossim <- plot_weighted_points_PA(PHONEME_cossim,"cossim_span",refdf1)
PA_p1_cossim_subset <- plot_weighted_points_PA(PHONEME_cossim_filter_test,"cossim_span",refdf1)
PA_p1_cossim_shuffled <- plot_weighted_points_PA(PHONEME_cossim_shuffled,"cossim_span",refdf1)

PA_p1_cossim + PA_p1_cossim_subset

PA_p2_cossim <- plot_weighted_points_PA(df_jap,"cossim",refdf1)
PA_p3_cossim <- plot_weighted_points_PA(df_eng,"cossim",refdf1)

GA_p1_cossim  <- plot_weighted_points_PA(GRAMMAR_cossim, 'cossim_span',refdf1)
GA_p2_cossim  <- plot_weighted_points_PA(GRAMMAR_cossim, 'cossim_jap',refdf1)
GA_p3_cossim  <- plot_weighted_points_PA(GRAMMAR_cossim, 'cossim_eng',refdf1)


GA_p1_cossim + GA_p2_cossim + GA_p3_cossim

PA_p1_sim <- plot_weighted_points_PA(df_span_ss,"sim",refdf1)
PA_p2_sim <- plot_weighted_points_PA(df_jap_ss,"sim",refdf1)
PA_p3_sim <- plot_weighted_points_PA(df_eng_ss,"sim",refdf1)



PA_p1 <- plot_weighted_points_PA(similarity_df_h1_span,"similarity_score_Spanish",refdf1,global_lim)
PA_p2 <- plot_weighted_points_PA(similarity_df_h2_span,"similarity_score_Spanish",refdf2,global_lim)

match_df_span <- match_df_span %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

similarity_df <- similarity_df %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))

plot_weighted_contour(similarity_df_h1_span, "similarity_score_Spanish")

PA_p1 + PA_p2

plot_weighted_points_PA(match_df_span,"zscore")
plot_weighted_contour(match_df_span,"zscore")


PA_p3 <- plot_weighted_points_PA(similarity_df, "mean_similarity", refdf1)
print(PA_p3)



PHOIBLE_z_score_df <- read_csv(here("data", "PHOIBLE_z_score_df.csv"))

PHOIBLE_z_score_df <- PHOIBLE_z_score_df %>% 
  mutate(latitude = lat.lang(language),
         longitude = long.lang(language))


all_vals <- c(PHOIBLE_z_score_df$zscore_span,
              PHOIBLE_z_score_df$zscore_eng,
              PHOIBLE_z_score_df$zscore_jap)

global_lim <- range(all_vals, na.rm = TRUE)

PA_p1_spm <- plot_weighted_points_PA(PHOIBLE_z_score_df,"zscore_span",lims = global_lim)
PA_p2_spm <- plot_weighted_points_PA(PHOIBLE_z_score_df,"zscore_eng",lims = global_lim)
PA_p3_spm <- plot_weighted_points_PA(PHOIBLE_z_score_df,"zscore_jap",lims = global_lim)

PA_p1_spm + PA_p2_spm + PA_p3_spm


GA_p1_test <- plot_weighted_points_PA(gramfeature_match_df_test_lang,"match_score_span")
GA_p1_test
