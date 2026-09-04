# =============================================================================
# Shared — Canonical subgroup -> colour palette
#
# The phoneme, grammar, cognate, and genetic pipelines each colour their
# languages/populations by Glottolog family subgroup (the 4th comma-separated
# element of `affiliation`). Previously each domain's [0]_Phylogenetic_Tree.R
# assigned its own Polychrome-36 palette in order of its own clade sizes, so the
# same subgroup got a different colour in every domain (e.g. "Cagayan Valley"
# was #79F287 in phoneme, #F16DEF in grammar). This derives one palette, keyed
# on Glottolog alone rather than on any domain's study set, so it is stable as
# study sets change and every domain's figures agree.
#
# Input:   lingtypology::glottolog (bundled data)
# Outputs: data/shared/subgroup_palette.csv (subgroup, n_glottolog, colour)
# Consumers: R/{phoneme,grammar,cognate,genetic}_analysis/[0]_*.R
# =============================================================================

library(dplyr)
library(stringr)
library(here)
library(lingtypology)

# ── 1. Every Philippine-scoped Glottolog subgroup ───────────────────────────
# `countries` is a 2-letter ISO code list ("PH" not "Philippines"); level is
# left unrestricted (language/dialect/family) to match the widest possible set
# a domain lookup could resolve to.
gl <- lingtypology::glottolog %>%
  filter(str_detect(coalesce(countries, ""), "PH"))

# subgroup = 4th element of the affiliation path; shorter paths fall back to
# their deepest level. Family-level rows are one step shallow relative to a
# language row of the same clade (their own name is the missing tail), so
# append it before slicing — this is what turns the `Subanen` family node into
# subgroup "Subanen" rather than the shallower "Greater Central Philippine".
sub4 <- function(affiliation, name, level) {
  mapply(function(a, n, l) {
    parts <- str_trim(str_split(a, ",")[[1]])
    if (l == "family") parts <- c(parts, n)
    if (length(parts) >= 4) parts[4] else tail(parts, 1)
  }, affiliation, name, level)
}
gl <- gl %>% mutate(subgroup = sub4(affiliation, language, level))

# Non-linguistic Glottolog buckets that show up in the PH-scoped rows but are
# not real subgroups any domain would map a language/population to.
drop_subgroups <- c("Bookkeeping", "ASLic", "Artificial Language", "Unclassifiable")

subgroup_counts <- gl %>%
  filter(!subgroup %in% drop_subgroups) %>%
  count(subgroup, sort = TRUE, name = "n_glottolog") %>%
  arrange(desc(n_glottolog), subgroup)

# ── 2. Palette ────────────────────────────────────────────────────────────
# Softened Polychrome-36, luminance-filtered to drop colours invisible on
# white — identical idiom to the one previously inlined in each domain's
# [0]_Phylogenetic_Tree.R, now computed once here.
soften <- function(cols, s_mult = 0.55, v_mult = 0.95) {
  h <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  grDevices::hsv(h["h", ], h["s", ] * s_mult, pmin(h["v", ] * v_mult, 1))
}
.poly <- grDevices::palette.colors(NULL, "Polychrome 36")
.lum  <- colSums(grDevices::col2rgb(.poly) * c(0.299, 0.587, 0.114))  # 0..255
.pal_slots <- soften(.poly[.lum < 200])

stopifnot(
  "More PH Glottolog subgroups than palette slots — widen the palette source." =
    nrow(subgroup_counts) <= length(.pal_slots)
)

subgroup_palette <- subgroup_counts %>%
  mutate(colour = .pal_slots[seq_len(n())])

write.csv(subgroup_palette, here("data", "shared", "subgroup_palette.csv"), row.names = FALSE)
