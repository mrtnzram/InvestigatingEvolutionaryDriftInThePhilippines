# =============================================================================
# Phoneme Analysis — Build EEMS input files (dissimilarity + coordinates)
# Converts the language-to-Spanish cosine-similarity scores into an EEMS
# dissimilarity matrix and coordinate file, plus a shuffled (null-model) variant.
# Outputs: data/datapath.diffs, data/datapath.coord, data/PHONEME_cossim_shuffled.csv
# =============================================================================

library(readr)
library(tidyverse)
library(dplyr)
library(reshape2)
library(here)

# ---- 1. Build the |score_i - score_j| dissimilarity matrix (ordered by longitude) ----
similarity_df <- read.csv(here("data", "PHONEME_cossim.csv"), header = TRUE)

similarity_df <- similarity_df %>%
  select(language, cossim_span, latitude, longitude)

similarity_df <- similarity_df[order(similarity_df$longitude), ]
scores <- similarity_df$cossim_span
names(scores) <- similarity_df$language

dist_matrix <- outer(scores, scores, function(x, y) abs(x - y))
diag(dist_matrix) <- 0
rownames(dist_matrix) <- colnames(dist_matrix) <- similarity_df$language

PHONEME_diss_matrix <- dist_matrix

write.table(dist_matrix, file = here("data", "datapath.diffs"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

# ---- 2. Preview the dissimilarity matrix as a heatmap ----
dist_long <- melt(dist_matrix)
names(dist_long) <- c("Lang1", "Lang2", "Dissimilarity")

dist_long$Lang1 <- factor(dist_long$Lang1, levels = similarity_df$language)
dist_long$Lang2 <- factor(dist_long$Lang2, levels = similarity_df$language)

ggplot(dist_long, aes(x = Lang1, y = Lang2, fill = Dissimilarity)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(name = "Dissimilarity", option = "C") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 3. Write the EEMS coordinate file (latitude, longitude) ----
similarity_df <- similarity_df[order(similarity_df$longitude), ]

coords <- similarity_df[, c("latitude", "longitude")]
write.table(coords, file = here("data", "datapath.coord"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# ---- 4. Shuffled null model: permute scores, rebuild dissimilarity, write outputs ----
similarity_df <- similarity_df %>%
  select(language, cossim_span, latitude, longitude)

similarity_df <- similarity_df[order(similarity_df$longitude), ]

set.seed(2025)
shuffled_scores <- sample(similarity_df$cossim_span)

shuffled_similarity_df <- similarity_df
shuffled_similarity_df$cossim_span <- shuffled_scores

scores <- shuffled_similarity_df$cossim_span
names(scores) <- shuffled_similarity_df$language

dist_matrix <- outer(scores, scores, function(x, y) abs(x - y))
diag(dist_matrix) <- 0
rownames(dist_matrix) <- colnames(dist_matrix) <- shuffled_similarity_df$language

write.table(dist_matrix, file = here("data", "datapath.diffs"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

write.csv(shuffled_similarity_df, file = here("data", "PHONEME_cossim_shuffled.csv"), row.names = FALSE)
