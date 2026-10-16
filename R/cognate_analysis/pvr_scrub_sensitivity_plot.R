# =============================================================================
# Cognate Analysis — scrub-step sensitivity figure
#
# Three rows, one per [4] scrub variant (signed axis + distance from Manila,
# signed axis only, distance only): the geography fit of phylogenetic Axis 1, and
# the normalized Spanish loanword count regressed on the resulting pure phylogeny.
# One figure per zero-loan condition.
#
# Input:   `removed`, `kept` from [4]_COGNATE_PVR.R
# Outputs: figures/regression/pvr/cognate_pvr_scrub_sensitivity[_zeros_kept].png
# =============================================================================

library(tidyverse)
library(here)
source(here("R", "shared", "pure_phylogeny.R"))

stopifnot("Run [4]_COGNATE_PVR.R first: `removed` / `kept` are not defined." =
            exists("removed") && exists("kept"))

# Same five languages and colours as the phoneme figures, matched on glottocode.
FOCUS_GLOTTO <- c(taga1270 = "Tagalog", pamp1243 = "Kapampangan", ilok1237 = "Ilokano",
                  cebu1242 = "Cebuano", sang1337 = "Sangil")
FOCUS <- c(Tagalog = "#D62728", Kapampangan = "#FF7F0E", Ilokano = "#1F77B4",
           Cebuano = "#2CA02C", Sangil = "#9467BD")

dir.create(here("figures", "regression", "pvr"), recursive = TRUE, showWarnings = FALSE)

for (r in list(list(fit = removed, sfx = ""), list(fit = kept, sfx = "_zeros_kept"))) {
  sens_fig <- plot_scrub_sensitivity(
    r$fit$frame |> mutate(focus = unname(FOCUS_GLOTTO[glottocode])),
    FOCUS, "Normalized Spanish loanword count"
  )
  print(sens_fig)
  ggsave(here("figures", "regression", "pvr", paste0("cognate_pvr_scrub_sensitivity", r$sfx, ".png")),
         sens_fig, width = 12, height = 13, units = "in", dpi = 300)
}
