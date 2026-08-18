# =============================================================================
# [4.2] Phoneme Analysis — true multivariate sPCA
#
# [4.1]'s companion: instead of a single phylogeny-scrubbed trait fed through
# Moran Eigenvector Maps (a univariate-safe stand-in), this residualizes the
# full 728-column binary phoneme matrix against the same phylogenetic
# eigenvectors and runs real adegenet::spca() on the resulting continuous
# residual matrix — spca() needs >=2 variables to rotate, which a single
# trait never satisfied.
#
# Runs at language grain (58), not tip grain (65): the raw binary matrix is
# one row per language, and duplicating rows per dialect tip (as [4]/[4.1] do
# for PVR's positional matching) creates near-duplicate rows that destabilize
# spca_randtest()'s permutation rank. E_sel is computed at tip grain (PVR
# requires the tree) and then averaged down to one row per language.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/RUHLENdf_PH.csv (full binary matrix), data/PHONEME_final.csv
#          (y/x for eigenvector selection only), data/PHONEME_dist_matrix.csv
# Outputs: data/PHONEME_sPCA_results.csv, data/PHONEME_sPCA_scores.csv,
#          data/PHONEME_sPCA_loadings.csv,
#          figures/phoneme/regression/phoneme_sPCA_surface.png,
#          data/base_plot_phoneme_sPCA.rds
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
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

PREDICTOR <- "geodist_H1_span"
N_PERM    <- 999   # + observed arrangement = 1000 draws

tip_map <- tree_df_matched |> dplyr::select(original, ph)


# ── 1. Phylogenetic eigenvectors at tip grain — same PVR() call as [4]/[4.1] ─
PHONEME_final <- read.csv(here("data", "PHONEME_final.csv"))

tip_df <- PHONEME_final |>
  dplyr::select(language, y = cossim_span_norm, x_km = all_of(PREDICTOR)) |>
  left_join(tip_map, by = c("language" = "ph")) |>
  filter(!is.na(original), !is.na(y), !is.na(x_km)) |>
  as.data.frame()
tip_df <- tip_df[match(tree_pruned$tip.label, tip_df$original), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the tip set." =
    nrow(tip_df) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(tip_df$original) && identical(tip_df$original, tree_pruned$tip.label)
)

pvr_dec <- PVRdecomp(tree_pruned, scale = TRUE)
pvr_fit <- PVR(pvr_dec, phy = tree_pruned, trait = tip_df$y, envVar = tip_df$x_km / 1000,
               method = "moran")
E_sel_tip <- as.matrix(pvr_fit@Selection$Vectors)
k         <- ncol(E_sel_tip)
message(k, " phylogenetic eigenvector", if (k == 1) "" else "s",
        " selected (matches [4]/[4.1]_PHONEME_PVR.R).")

# Collapse to one row per language (mean across a language's tips).
E_sel <- as.data.frame(E_sel_tip) |>
  mutate(language = tip_df$language) |>
  summarise(across(where(is.numeric), mean), .by = language) |>
  column_to_rownames("language") |>
  as.matrix()


# ── 2. Analysis frame: full binary matrix at language grain ─────────────────
phn_cols <- grep("^phoneme_", names(read.csv(here("data", "RUHLENdf_PH.csv"), nrows = 1)),
                 value = TRUE)
df <- read.csv(here("data", "RUHLENdf_PH.csv")) |>
  filter(Language_type == "Philippine Language") |>
  dplyr::select(language, longitude, latitude, all_of(phn_cols)) |>
  filter(language %in% rownames(E_sel)) |>
  arrange(match(language, rownames(E_sel)))

stopifnot(
  "Binary-matrix language set disagrees with E_sel's." =
    identical(df$language, rownames(E_sel))
)
N <- nrow(df)
message(N, " languages, ", length(phn_cols), " phoneme columns.")


# ── 3. Per-column residualization ────────────────────────────────────────────
X <- as.matrix(df[, phn_cols])
X <- X[, apply(X, 2, var) > 0]   # drop invariant columns (no signal to scrub or rotate)
# k = 0 (lm() rejects a zero-column matrix term) means nothing to scrub.
R <- if (ncol(E_sel) == 0) X else apply(X, 2, \(col) resid(lm(col ~ E_sel)))
message(ncol(X), " of ", length(phn_cols), " phoneme columns retained (non-invariant); ",
        "residualized against the ", k, "-eigenvector set above.")


# ── 4. Spatial weight matrix: threshold search, maximizing sPCA eigenvalue ──
PHONEME_dist_matrix <- read.csv(here("data", "PHONEME_dist_matrix.csv"),
                                row.names = 1, check.names = FALSE) |>
  as.matrix()
Dgeo <- PHONEME_dist_matrix[df$language, df$language]
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
write.csv(scores_df, file = here("data", "PHONEME_sPCA_scores.csv"), row.names = FALSE)

loadings_df <- tibble(phoneme = colnames(R), loading = spca_fit$c1[, 1]) |>
  arrange(desc(abs(loading)))
write.csv(loadings_df, file = here("data", "PHONEME_sPCA_loadings.csv"), row.names = FALSE)


# ── 7. Results table ─────────────────────────────────────────────────────────
PHONEME_sPCA_results <- tibble(
  n = N, n_evec_phylo = k, n_features_retained = ncol(R),
  threshold_km = best$threshold, eigenvalue = best$eig1,
  variance_explained = var_explained, perm_p = perm_p, n_perm = N_PERM
)
print(PHONEME_sPCA_results)
write.csv(PHONEME_sPCA_results, file = here("data", "PHONEME_sPCA_results.csv"), row.names = FALSE)


# ── 8. Plot: sPCA point-symbol map (size = |sPC1|, colour = sign) ───────────
dir.create(here("figures", "phoneme", "regression"), recursive = TRUE, showWarnings = FALSE)

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
  labs(title = "Phoneme sPCA: phylogeny-scrubbed multivariate spatial structure",
       subtitle = sprintf("permutation p = %.3f (%d permutations) | variance explained = %.3f",
                          perm_p, N_PERM + 1, var_explained))
print(p_surface)

ggsave(here("figures", "phoneme", "regression", "phoneme_sPCA_surface.png"),
       p_surface, width = 7.5, height = 6, units = "in", dpi = 300)
saveRDS(p_surface, file = here("data", "base_plot_phoneme_sPCA.rds"))
