# =============================================================================
# Phoneme Analysis — scrub-step sensitivity figure
#
# Three rows, one per [4] scrub variant (signed axis + distance from Manila,
# signed axis only, distance only): the geography fit of phylogenetic Axis 1, and
# Spanish phoneme similarity regressed on the resulting pure phylogeny.
#
# Input:   `df` from [4]_PHONEME_PVR.R
# Outputs: figures/regression/pvr/phoneme_pvr_scrub_sensitivity.png
# =============================================================================

library(tidyverse)
library(here)
source(here("R", "shared", "pure_phylogeny.R"))

stopifnot("Run [4]_PHONEME_PVR.R first: `df` lacks the pure-phylogeny columns." =
            exists("df") && all(c("phylo_axis1", "geo_signed_km", "x_km", "y") %in% names(df)))

FOCUS <- c(Tagalog = "#D62728", Kapampangan = "#FF7F0E", Ilokano = "#1F77B4",
           Cebuano = "#2CA02C", Sangil = "#9467BD")

sens_fig <- plot_scrub_sensitivity(
  df |> mutate(focus = if_else(language %in% names(FOCUS), language, NA_character_)),
  FOCUS, "Normalized Spanish phoneme similarity"
)
print(sens_fig)

dir.create(here("figures", "regression", "pvr"), recursive = TRUE, showWarnings = FALSE)
ggsave(here("figures", "regression", "pvr", "phoneme_pvr_scrub_sensitivity.png"),
       sens_fig, width = 12, height = 13, units = "in", dpi = 300)
