# =============================================================================
# GRAMBANK feature set — the 50 features shared by every grammar analysis
#
# Derived once by the two-pass coverage filter in
# R/grammar_analysis/[0]_GRAMBANKdatabase.R (feature_thresh = 80 languages,
# language_thresh = 50 features, capped at 2 iterations). Those thresholds were
# applied to the Philippine bounding-box set, so the result is specific to that
# starting pool: re-running the loop on a different language set — Austronesian-
# wide, for instance — yields a different 50. The IDs are therefore frozen here
# and sourced, never re-derived downstream.
#
# Of the 50, GB065 and GB130 are 3-state with levels {1, 2, 3}; the other 48 are
# binary {0, 1}. Note that neither multistate feature has a 0 state, so nothing
# downstream may assume states are 0-indexed.
# =============================================================================

GRAMBANK_FEATURE_IDS <- c(
  "GB020", "GB021", "GB022", "GB023", "GB028", "GB030", "GB031", "GB035",
  "GB036", "GB037", "GB042", "GB043", "GB044", "GB051", "GB052", "GB053",
  "GB054", "GB065", "GB070", "GB071", "GB072", "GB073", "GB079", "GB080",
  "GB082", "GB083", "GB084", "GB086", "GB089", "GB090", "GB091", "GB092",
  "GB093", "GB094", "GB107", "GB121", "GB130", "GB131", "GB137", "GB138",
  "GB171", "GB172", "GB186", "GB192", "GB196", "GB197", "GB316", "GB318",
  "GB321", "GB415"
)

stopifnot(
  "the GRAMBANK feature set must hold 50 unique IDs" =
    length(GRAMBANK_FEATURE_IDS) == 50L && !anyDuplicated(GRAMBANK_FEATURE_IDS)
)
