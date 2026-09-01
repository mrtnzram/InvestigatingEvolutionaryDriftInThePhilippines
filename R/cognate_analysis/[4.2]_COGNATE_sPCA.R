# =============================================================================
# [4.2] Cognate Analysis — true multivariate sPCA
#
# Builds a one-hot (language x gloss-wordform) presence matrix from the raw
# ABVD wordforms in PH_df.csv, residualizes it against phylogenetic
# eigenvectors, and runs real adegenet::spca() on the result. Runs at language
# grain (94): [0]_Phylogenetic_Tree.R already reduces to one tip per language,
# so there is no dialect-tip collapsing step here unlike phoneme/grammar.
#
# Run order: requires `tree_pruned` and `tip_map` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/cognate/PH_df.csv (raw per-gloss wordforms),
#          data/cognate/COGNATE_final.csv (y/x for eigenvector selection only),
#          data/cognate/COGNATE_dist_matrix.csv
# Outputs: data/cognate/COGNATE_sPCA_results.csv, data/cognate/COGNATE_sPCA_scores.csv,
#          data/cognate/COGNATE_sPCA_loadings.csv,
#          figures/cognate/regression/cognate_sPCA_surface.png,
#          data/cognate/base_plot_cognate_sPCA.rds
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)
library(spdep)
library(adespatial)
library(adegenet)
library(maps)

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tip_map` is not defined." =
    exists("tip_map")
)

PREDICTOR <- "geodist_H1_span"
N_PERM    <- 999   # + observed arrangement = 1000 draws
MIN_LANGS <- 2      # drop (gloss, wordform) pairs attested in fewer languages


# ── 1. Phylogenetic eigenvectors — same PVR() call as [4]/[4.1] ─────────────
COGNATE_final <- read.csv(here("data", "cognate", "network_distance", "COGNATE_final.csv"))

tip_df <- COGNATE_final |>
  dplyr::select(glottocode, y = loans_norm, x_km = all_of(PREDICTOR)) |>
  semi_join(tip_map, by = "glottocode") |>
  filter(!is.na(y), !is.na(x_km)) |>
  as.data.frame()
tip_df <- tip_df[match(tree_pruned$tip.label, tip_df$glottocode), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the tip set." =
    nrow(tip_df) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(tip_df$glottocode) && identical(tip_df$glottocode, tree_pruned$tip.label)
)

pvr_dec <- PVRdecomp(tree_pruned, scale = TRUE)
pvr_fit <- PVR(pvr_dec, phy = tree_pruned, trait = tip_df$y, envVar = tip_df$x_km / 1000,
               method = "moran")
E_sel <- as.matrix(pvr_fit@Selection$Vectors)
rownames(E_sel) <- tip_df$glottocode
k <- ncol(E_sel)
message(k, " phylogenetic eigenvector", if (k == 1) "" else "s",
        " selected (matches [4]/[4.1]_COGNATE_PVR.R).")


# ── 2. Analysis frame: one-hot (gloss, wordform) matrix at language grain ───
META_COLS <- c("language_id", "language", "latitude", "longitude",
               "glottocode", "author", "number_of_entries", "source_url")
PH_df <- read.csv(here("data", "cognate", "initial_datasets", "PH_df.csv"))
gloss_cols <- setdiff(names(PH_df), META_COLS)

ph_sub <- PH_df |>
  dplyr::select(glottocode, language, longitude, latitude, all_of(gloss_cols)) |>
  filter(glottocode %in% rownames(E_sel))

# Cells are ";"-joined wordform variants; split, trim, casefold.
long <- ph_sub |>
  dplyr::select(glottocode, all_of(gloss_cols)) |>
  pivot_longer(-glottocode, names_to = "gloss", values_to = "cell") |>
  filter(!is.na(cell), cell != "") |>
  separate_longer_delim(cell, "; ") |>
  mutate(form = str_trim(tolower(cell))) |>
  filter(form != "") |>
  distinct(glottocode, gloss, form)

# A wordform attested by only one language can't carry spatial signal, so
# singletons are dropped before the matrix is built.
pair_counts <- long |> count(gloss, form, name = "n_langs")
kept_pairs  <- pair_counts |> filter(n_langs >= MIN_LANGS)
long <- long |> semi_join(kept_pairs, by = c("gloss", "form"))

message(nrow(pair_counts), " unique (gloss, wordform) pairs observed; ",
        nrow(kept_pairs), " attested in >= ", MIN_LANGS, " languages, kept.")

feature_wide <- long |>
  mutate(feature = paste(gloss, form, sep = "::"), present = 1L) |>
  distinct(glottocode, feature, present) |>
  pivot_wider(names_from = feature, values_from = present, values_fill = 0L)
feature_cols <- setdiff(names(feature_wide), "glottocode")

df <- ph_sub |>
  dplyr::select(glottocode, language, longitude, latitude) |>
  left_join(feature_wide, by = "glottocode") |>
  mutate(across(all_of(feature_cols), \(x) replace_na(x, 0L))) |>
  arrange(match(glottocode, rownames(E_sel))) |>
  as.data.frame()

stopifnot(
  "One-hot matrix language set disagrees with E_sel's." =
    identical(df$glottocode, rownames(E_sel))
)
N <- nrow(df)
message(N, " languages, ", length(feature_cols), " (gloss, wordform) features.")


# ── 3. Per-column residualization ────────────────────────────────────────────
X <- as.matrix(df[, feature_cols])
X <- X[, apply(X, 2, var) > 0]   # defensive: drop any column constant on this N

# Duplicate presence columns (sister languages sharing a subgroup or a loan)
# push the residual matrix's rank to ade4's zero-eigenvalue threshold and
# crash spca_randtest()'s permutations with mismatched eigenvalue-vector
# lengths; one representative per pattern removes the ambiguity.
X <- X[, !duplicated(t(X)), drop = FALSE]

R <- if (ncol(E_sel) == 0) X else apply(X, 2, \(col) resid(lm(col ~ E_sel)))
message(ncol(X), " of ", length(feature_cols), " features retained (non-invariant, deduplicated); ",
        "residualized against the ", k, "-eigenvector set above.")


# ── 4. Spatial weight matrix: threshold search, maximizing sPCA eigenvalue ──
COGNATE_dist_matrix <- read.csv(here("data", "cognate", "network_distance", "COGNATE_dist_matrix.csv"),
                                row.names = 1, check.names = FALSE) |>
  as.matrix()
Dgeo <- COGNATE_dist_matrix[df$glottocode, df$glottocode]
diag(Dgeo) <- 0

dvals <- Dgeo[upper.tri(Dgeo)]
dvals <- dvals[dvals > 0]
thresholds <- sort(unique(quantile(dvals, probs = seq(0.1, 1, length.out = 20), na.rm = TRUE)))

# PCA step is network-independent, so it runs once outside the threshold loop.
pca_R <- dudi.pca(as.data.frame(R), center = TRUE, scale = FALSE, scannf = FALSE)

best <- list(threshold = NA_real_, eig1 = -Inf, listw = NULL)
for (th in thresholds) {
  W <- 1 / Dgeo^2
  W[!is.finite(W)] <- 0
  W[Dgeo > th] <- 0
  diag(W) <- 0
  # spca()'s matWeight normalization (prop.table) can't tolerate an isolated
  # row, unlike spdep's zero.policy — reject any candidate with one, not just
  # a fully-disconnected graph.
  if (any(rowSums(W) == 0)) next

  lw <- mat2listw(W, style = "W", zero.policy = TRUE)
  ms <- tryCatch(multispati(pca_R, lw, scannf = FALSE, nfposi = 1, nfnega = 0),
                 error = function(e) NULL)
  if (is.null(ms)) next

  if (ms$eig[1] > best$eig1) {
    best <- list(threshold = th, eig1 = ms$eig[1], listw = lw)
  }
}
stopifnot("No candidate threshold produced a usable spatial weights object." =
            !is.null(best$listw))

message(sprintf("Spatial weights: threshold = %.1f km, leading eigenvalue = %.4f.",
                best$threshold, best$eig1))


# ── 5. Real spca() ────────────────────────────────────────────────────────
W_best <- 1 / Dgeo^2
W_best[!is.finite(W_best)] <- 0
W_best[Dgeo > best$threshold] <- 0
diag(W_best) <- 0

spca_fit <- spca(as.data.frame(R), xy = cbind(df$longitude, df$latitude),
                 matWeight = W_best, scannf = FALSE, nfposi = 1, nfnega = 1)
sPC1 <- spca_fit$li[, 1]

perm <- spca_randtest(spca_fit, nperm = N_PERM)
# NOT perm$global$pvalue — that field is broken upstream and always returns 1:
# spca_randtest passes `sim = sum(sims[sims >= 0])`, a scalar summed over every
# permutation, so obs can never exceed it. eigentest uses a proper per-permutation
# vector (`sims[e, ]`); [1, ] is the leading axis.
perm_p <- perm$eigentest[1, "sim_p"]

var_explained <- best$eig1 / sum(abs(pca_R$eig))

message(sprintf("sPCA: variance explained (axis 1) = %.3f, permutation p = %.4f.",
                var_explained, perm_p))


# ── 6. Per-language scores + per-feature loadings ────────────────────────────
scores_df <- tibble(language = df$language, longitude = df$longitude,
                    latitude = df$latitude, sPC1 = sPC1)
write.csv(scores_df, file = here("data", "cognate", "PVR", "COGNATE_sPCA_scores.csv"),
          row.names = FALSE)

loadings_df <- tibble(feature = colnames(R), loading = spca_fit$c1[, 1]) |>
  separate_wider_delim(feature, "::", names = c("gloss", "wordform"), cols_remove = FALSE) |>
  arrange(desc(abs(loading)))
write.csv(loadings_df, file = here("data", "cognate", "PVR", "COGNATE_sPCA_loadings.csv"),
          row.names = FALSE)


# ── 7. Results table ─────────────────────────────────────────────────────────
COGNATE_sPCA_results <- tibble(
  n = N, n_evec_phylo = k, n_features_retained = ncol(R),
  threshold_km = best$threshold, eigenvalue = best$eig1,
  variance_explained = var_explained, perm_p = perm_p, n_perm = N_PERM
)
print(COGNATE_sPCA_results)
write.csv(COGNATE_sPCA_results,
          file = here("data", "cognate", "PVR", "COGNATE_sPCA_results.csv"), row.names = FALSE)


# ── 8. Plot: sPCA point-symbol map (size = |sPC1|, colour = sign) ───────────
dir.create(here("figures", "cognate", "regression"), recursive = TRUE, showWarnings = FALSE)

world_map  <- map_data("world")
map_subset <- world_map |> filter(region %in% c("Philippines", "Malaysia"))

mag_threshold <- median(abs(scores_df$sPC1))
scores_df <- scores_df |>
  mutate(
    sign      = if_else(sPC1 >= 0, "Positive", "Negative"),
    magnitude = if_else(abs(sPC1) >= mag_threshold, "Large", "Small"),
    key       = factor(paste(magnitude, sign),
                       levels = c("Large Positive", "Small Positive",
                                  "Small Negative", "Large Negative"))
  )

key_sizes  <- c("Large Positive" = 10, "Small Positive" = 4,
                "Small Negative" = 4,  "Large Negative" = 10)
key_fills  <- c("Large Positive" = "black", "Small Positive" = "black",
                "Small Negative" = "white", "Large Negative" = "white")
key_labels <- c("Large Positive" = "Large (+)", "Small Positive" = "Small (+)",
                "Small Negative" = "Small (-)", "Large Negative" = "Large (-)")

p_surface <- ggplot() +
  geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
               fill = "gray97", color = "black") +
  geom_point(data = scores_df,
             aes(x = longitude, y = latitude, size = key, fill = key),
             shape = 22, colour = "black", stroke = 0.6) +
  scale_size_manual(values = key_sizes, labels = key_labels, name = "sPC1", drop = FALSE) +
  scale_fill_manual(values = key_fills, labels = key_labels, name = "sPC1", drop = FALSE) +
  coord_fixed(xlim = c(115, 130), ylim = c(4, 22)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text  = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "Cognate sPCA",
       subtitle = sprintf("p = %.3f, r^2 = %.3f", perm_p, var_explained))
print(p_surface)

ggsave(here("figures", "cognate", "regression", "cognate_sPCA_surface.png"),
       p_surface, width = 7.5, height = 6, units = "in", dpi = 300)
saveRDS(p_surface, file = here("data", "cognate", "PVR", "base_plot_cognate_sPCA.rds"))
