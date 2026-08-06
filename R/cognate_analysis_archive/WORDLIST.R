library(jsonlite)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(glue)
library(readr)
library(tibble)
library(here)

# Set your folder path
folder_path <- here("swadeshlist_jsons")

# Get all JSON files in the folder
json_files <- list.files(folder_path, pattern = "\\.json$", full.names = TRUE)

# Function to extract wordlist from a single JSON file
extract_wordlist <- function(file_path) {
  # Read JSON
  json_data <- fromJSON(file_path)
  
  # Get language name from filename
  language <- tools::file_path_sans_ext(basename(file_path))
  language <- gsub("_", " ", language)
  
  # Extract lines from txt field
  lines <- unlist(strsplit(json_data$txt, "\n"))
  word_lines <- grep("^[0-9]+\\s+\\S+\\t", lines, value = TRUE)
  
  # Parse lines
  parsed <- lapply(word_lines, function(line) {
    parts <- unlist(strsplit(line, "\\t"))
    if (length(parts) >= 2) {
      id_concept <- unlist(strsplit(parts[1], "\\s+", perl = TRUE))
      concept <- id_concept[2]
      words <- gsub(" //", "", parts[2])
      word_list <- unlist(strsplit(words, ",\\s*"))
      # Return one row per word
      data.frame(concept = concept, word = word_list, stringsAsFactors = FALSE)
    }
  })
  
  # Combine and add language column
  df <- bind_rows(parsed)
  df$language <- language
  return(df)
}

# Extract all wordlists
all_wordlists <- map_dfr(json_files, extract_wordlist)

# Choosing only 1 if there are multiple
collapsed <- all_wordlists %>%
  group_by(concept, language) %>%
  summarise(word = first(word), .groups = "drop")

# Pivot to wide format: one row per concept, one column per language
wordlist_wide <- collapsed %>%
  pivot_wider(names_from = language, values_from = word)

# Transpose the dataframe: concepts become column names, languages become rows
wordlist_long <- wordlist_wide %>%
  pivot_longer(-concept, names_to = "language", values_to = "word")

# Now pivot back to wide: one row per language, one column per concept
wordlist_by_language <- wordlist_long %>%
  pivot_wider(names_from = concept, values_from = word)

protected_languages <- c("SPANISH", "ENGLISH", "JAPANESE")
unrelated_languages <- c("BASQUE","HUNGARIAN","FINNISH","SOUTH SAAMI","GAGAUZ","ESTONIAN","TURKISH","MOKSHA","UGARITIC","LIVE",
                         "KOREAN","FUZHOU CHINESE","MONGOLIAN","SEDANG")
always_keep <- unique(c(protected_languages, unrelated_languages))

# Parameters
feature_thresh <- 30   # Minimum number of languages a concept must appear in
language_thresh <- 30  # Minimum number of concepts a language must have

# Normalize row names (language names) for matching
rownames(wordlist_by_language) <- toupper(rownames(wordlist_by_language))

# Initialize
wordlist_snapshots <- list()
iteration <- 1

repeat {
  cat(glue("\n--- Iteration {iteration} ---\n"))
  
  old_dim <- dim(wordlist_by_language)
  
  # Identify concept columns (exclude rownames which are languages)
  concept_cols <- colnames(wordlist_by_language)
  
  # Filter concepts: keep those present in ≥ feature_thresh languages
  keep_concepts <- concept_cols[colSums(!is.na(wordlist_by_language[, concept_cols, drop = FALSE])) >= feature_thresh]
  wordlist_by_language <- wordlist_by_language[, keep_concepts, drop = FALSE]
  cat(glue("Retained {length(keep_concepts)} concepts\n"))
  
  # Filter languages: keep those with ≥ language_thresh non-NA concepts
  lang_filter <- rowSums(!is.na(wordlist_by_language[, keep_concepts, drop = FALSE])) >= language_thresh
  
  # Always retain protected languages
  lang_filter <- lang_filter | rownames(wordlist_by_language) %in% always_keep
  wordlist_by_language <- wordlist_by_language[lang_filter, , drop = FALSE]
  cat(glue("Retained {nrow(wordlist_by_language)} languages\n"))
  
  # Save snapshot
  wordlist_snapshots[[paste0("iter_", iteration)]] <- wordlist_by_language
  
  # Stop if shape is stable
  if (all(dim(wordlist_by_language) == old_dim)) {
    cat(glue("Matrix stabilized after iteration {iteration}.\n"))
    break
  }
  
  # Stop if pruning is too aggressive
  if (nrow(wordlist_by_language) < 2 || length(keep_concepts) < 2) {
    warning("Matrix reduced to nearly nothing — adjust your thresholds.")
    break
  }
  
  iteration <- iteration + 1
}

WORDLISTdf_maximized <- wordlist_snapshots$iter_2

write_csv(WORDLISTdf_maximized, here("R", "cognate_analysis", "WORDLISTdf.csv"))


# Pivot longer — from wide (concepts as columns) to long format
long_format <- WORDLISTdf_maximized %>%
  pivot_longer(cols = -language, names_to = "concept", values_to = "word")

# Pivot wider — from long to wide with concepts as rows and languages as columns
WORDLISTdf_concepts <- long_format %>%
  pivot_wider(names_from = language, values_from = word)

write_csv(WORDLISTdf_concepts, here("R", "cognate_analysis", "WORDLISTdf_forlingpy.csv"))
