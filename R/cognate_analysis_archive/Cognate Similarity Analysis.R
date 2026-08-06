library(readr)
library(tidyverse)
library(tibble)
library(lingtypology)
library(here)

WORDLIST_PH_cognates <- read_csv(here("R", "cognate_analysis", "WORDLIST_PH_cognates_unfiltered.csv"))

COGNATES_df <- WORDLIST_PH_cognates %>%
  drop_na() %>% 
  select(CONCEPT,DOCULECT,COGID) %>% 
  pivot_wider(names_from = "CONCEPT",values_from = "COGID") %>% 
  rename(language = DOCULECT)

interest_languages <- c("SPANISH", "ENGLISH", "JAPANESE")
unrelated_languages <- c("BASQUE","HUNGARIAN","FINNISH","SOUTH SAAMI","GAGAUZ","ESTONIAN","TURKISH","MOKSHA","UGARITIC","LIVE",
                         "KOREAN","FUZHOU CHINESE","MONGOLIAN","SEDANG")

COGNATES_df_unrelated <- COGNATES_df %>% 
  filter(language %in% unrelated_languages)

COGNATES_df_interest <- COGNATES_df %>% 
  filter(language %in% interest_languages)

COGNATES_df_PH <- COGNATES_df %>% 
  filter(!language %in% c(unrelated_languages,interest_languages))


comparison_results <- list()  # to store outputs for each interest language

concept_columns <- setdiff(names(COGNATES_df_interest), "language")

for (lang in interest_languages) {
  # Extract the target interest language row
  target_row <- COGNATES_df_interest %>% filter(language == lang)
  
  # Replicate target row to match PH dataframe size
  target_repeated <- target_row[rep(1, nrow(COGNATES_df_PH)), ]
  
  # Compare concept columns
  comparison_df <- COGNATES_df_PH
  comparison_df[concept_columns] <- mapply(function(ph, ref) ifelse(ph == ref, 1, 0),
                                           comparison_df[concept_columns],
                                           target_repeated[concept_columns])
  # Save result
  comparison_results[[lang]] <- comparison_df
}

# Count number of 1s per row across those columns
comparison_results$SPANISH$similiarity_to_spanish <- apply(
  comparison_results$SPANISH[concept_columns],
  1,
  function(row) {
    match_values <- row[!is.na(row)]         # exclude NAs
    num_ones <- sum(match_values == 1)
    num_comparables <- sum(match_values == 1 | match_values == 0)
    if (num_comparables == 0) NA else (num_ones / num_comparables) * 100
  }
)

comparison_results$ENGLISH$similiarity_to_english <- apply(
  comparison_results$ENGLISH[concept_columns],
  1,
  function(row) {
    match_values <- row[!is.na(row)]         # exclude NAs
    num_ones <- sum(match_values == 1)
    num_comparables <- sum(match_values == 1 | match_values == 0)
    if (num_comparables == 0) NA else (num_ones / num_comparables) * 100
  }
)

comparison_results$JAPANESE$similiarity_to_japanese <- apply(
  comparison_results$JAPANESE[concept_columns],
  1,
  function(row) {
    match_values <- row[!is.na(row)]         # exclude NAs
    num_ones <- sum(match_values == 1)
    num_comparables <- sum(match_values == 1 | match_values == 0)
    if (num_comparables == 0) NA else (num_ones / num_comparables) * 100
  }
)

comparison_results$SPANISH$n_span <- apply(
  comparison_results$SPANISH[concept_columns],
  1,
  function(row) sum(row %in% c(0, 1), na.rm = TRUE)
)

comparison_results$ENGLISH$n_eng <- apply(
  comparison_results$ENGLISH[concept_columns],
  1,
  function(row) sum(row %in% c(0, 1), na.rm = TRUE)
)

comparison_results$JAPANESE$n_jap <- apply(
  comparison_results$JAPANESE[concept_columns],
  1,
  function(row) sum(row %in% c(0, 1), na.rm = TRUE)
)

Spanish_similarity <- comparison_results$SPANISH %>% 
  select(language,similiarity_to_spanish,n_span)
English_similarity <- comparison_results$ENGLISH %>% 
  select(language,similiarity_to_english,n_eng)
Japanese_similarity <- comparison_results$JAPANESE %>% 
  select(language,similiarity_to_japanese,n_jap)

Similarity_df <- Spanish_similarity %>% 
  left_join(English_similarity,by = "language") %>% 
  left_join(Japanese_similarity,by = "language")

Similarity_df <- Similarity_df %>% 
  filter(!language %in% c('ZAMBOANGUENO','CAIJIA'))




# Similarity for Unrelated Languages

# Initialize vector to hold percent match scores
unrelated_match_percentages <- c()

# Initialize list to collect distributions for each PH language
results_list <- list()

for (ph_lang in COGNATES_df_PH$language) {
  # Extract the PH language row
  ph_row <- COGNATES_df_PH %>% filter(language == ph_lang)
  
  # Initialize vector to store match percentages across unrelated languages
  match_percentages <- c()
  
  for (unrel_lang in unrelated_languages) {
    # Get unrelated language row
    ref_row <- COGNATES_df %>% filter(language == unrel_lang)
    if (nrow(ref_row) == 0) next  # skip if missing
    
    # Compare concept columns row-wise
    ph_values <- ph_row[concept_columns]
    ref_values <- ref_row[concept_columns]
    
    binary_vector <- mapply(function(ph, ref) {
      if (is.na(ph) || is.na(ref)) 0 else if (ph == ref) 1 else 0
    }, ph_values, ref_values)
    
    # Compute percentage match
    valid <- binary_vector[!is.na(binary_vector)]
    if (length(valid) == 0) next
    percent_match <- mean(valid) * 100
    match_percentages <- c(match_percentages, percent_match)
  }
  
  # Store results per PH language
  results_list[[ph_lang]] <- list(
    language = ph_lang,
    Mean = mean(match_percentages, na.rm = TRUE),
    Std_dev = sd(match_percentages, na.rm = TRUE),
    Distribution = match_percentages
  )
}

# Convert results list to a dataframe
results_df <- tibble(
  language = names(results_list),
  Unrelated_Mean = sapply(results_list, function(x) x$Mean),
  Unrelated_Std_dev = sapply(results_list, function(x) x$Std_dev),
  Unrelated_Distribution = lapply(results_list, function(x) x$Distribution)
)

Similarity_df <- Similarity_df %>% 
  left_join(results_df,by = "language")



target <- Similarity_df %>% filter(language == "CEBUANO")

# Create a density-friendly dataframe
plot_data <- data.frame(score = unlist(target$Unrelated_Distribution))

# Plot the distribution with reference lines
ref_lines <- data.frame(
  similarity = c(
    target$similiarity_to_spanish,
    target$similiarity_to_english,
    target$similiarity_to_japanese
  ),
  source = factor(c("Spanish", "English", "Japanese"),
                  levels = c("Spanish", "English", "Japanese"))
)

# Plot the density with mapped reference lines
ggplot(plot_data, aes(x = score)) +
  geom_density(fill = "#56B4E9", alpha = 0.6) +
  geom_vline(data = ref_lines, aes(xintercept = similarity, color = source),
             linetype = "dashed", size = 1.2, show.legend = TRUE) +
  scale_color_manual(values = c(
    "Spanish" = "#D55E00",
    "English" = "#009E73",
    "Japanese" = "#F0E442"
  )) +
  labs(
    title = "Cognate Similarity Distribution vs. Unrelated Languages",
    x = "Cognate Similarity (%)",
    y = "Density",
    color = "Interest Language",
    subtitle = "Cebuano"
  ) +
  theme_bw()


long_Similarity_df <- bind_rows(
  Similarity_df %>% select(similarity = Unrelated_Mean) %>% mutate(language = "Unrelated"),
  Similarity_df %>% select(similarity = similiarity_to_spanish) %>% mutate(language = "Spanish"),
  Similarity_df %>% select(similarity = similiarity_to_english) %>% mutate(language = "English"),
  Similarity_df %>% select(similarity = similiarity_to_japanese) %>% mutate(language = "Japanese")
)
means_df <- long_Similarity_df %>%
  group_by(language) %>%
  summarise(mean_similarity = mean(similarity, na.rm = TRUE))

ggplot(long_Similarity_df, aes(x = similarity, fill = language)) +
  geom_density(alpha = 0.25) +
  scale_fill_manual(values = c(
    "Unrelated" = "#D3D3D3",
    "Spanish" = "#D55E00",
    "English" = "#009E73",
    "Japanese" = "#F0E442"
  )) +
  geom_vline(data = means_df, aes(xintercept = mean_similarity, color = language),
             linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c(
    "Unrelated" = "#D3D3D3",
    "Spanish" = "#D55E00",
    "English" = "#009E73",
    "Japanese" = "#F0E442"
  )) +
  labs(
    title = "Cognate Similarity Distribution with Means",
    x = "Cognate Similarity (%)",
    y = "Density",
    fill = "Language",
    color = "Mean"
  ) +
  theme_bw()

Similarity_df <- Similarity_df %>%
  mutate(
    zscore_span = (similiarity_to_spanish - Unrelated_Mean) / Unrelated_Std_dev,
    zscore_eng  = (similiarity_to_english - Unrelated_Mean) / Unrelated_Std_dev,
    zscore_jap  = (similiarity_to_japanese - Unrelated_Mean) / Unrelated_Std_dev,
    ) %>% 
  mutate(
    zscore_span = ifelse(zscore_span < 0, 0, zscore_span),
    zscore_eng = ifelse(zscore_eng < 0, 0, zscore_eng),
    zscore_jap = ifelse(zscore_jap <0, 0, zscore_jap)
  )

write_csv(Similarity_df, here("R", "cognate_analysis", "Similarity_df.csv"))

#map.feature(Similarity_df$language,
#            Similarity_df$zscore_span)
