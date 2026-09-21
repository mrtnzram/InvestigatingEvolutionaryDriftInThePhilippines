# =============================================================================
# [4.2] Grammar Analysis — true multivariate sPCA
#
# Residualizes the full 50-column binary Grambank matrix against phylogenetic
# eigenvectors and runs real adegenet::spca() on the result, recovering a true
# multivariate spatial component rather than a single-trait MEM approximation.
# Runs at language grain (26, not tip grain): E_sel is computed at tip grain
# (PVR requires the tree) then averaged down to one row per language, since
# duplicating rows per dialect tip would destabilize spca_randtest()'s
# permutation rank.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/grammar/GRAMBANKdf_full.csv (full binary matrix)
#          data/network_distance/GRAMMAR_final.csv (y/x for eigenvector selection only)
#          data/network_distance/GRAMMAR_dist_matrix.csv
# Outputs: data/spca/GRAMMAR_sPCA_results.csv, data/spca/GRAMMAR_sPCA_scores.csv
#          data/spca/GRAMMAR_sPCA_loadings.csv
#          data/spca/base_plot_grammar_sPCA.rds
#          figures/regression/grammar_sPCA_surface.png
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)
library(spdep)
library(adespatial)
library(adegenet)
library(maps)

source(here("R", "shared", "spatial_threshold.R"))

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

PREDICTOR <- "geodist_H1_span"
N_PERM    <- 999   # + observed arrangement = 1000 draws

tip_map <- tree_df_matched |> dplyr::select(original, gram)


# ── 1. Phylogenetic eigenvectors at tip grain — same PVR() call as [4]/[4.1] ─
GRAMMAR_final <- read.csv(here("data", "network_distance", "GRAMMAR_final.csv"))

tip_df <- GRAMMAR_final |>
  dplyr::select(language, y = cossim_span_norm, x_km = all_of(PREDICTOR)) |>
  left_join(tip_map, by = c("language" = "gram")) |>
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
        " selected (matches [4]/[4.1]_GRAMMAR_PVR.R).")

# Collapse to one row per language (mean across a language's tips).
E_sel <- as.data.frame(E_sel_tip) |>
  mutate(language = tip_df$language) |>
  summarise(across(where(is.numeric), mean), .by = language) |>
  column_to_rownames("language") |>
  as.matrix()


# ── 2. Analysis frame: full binary matrix at language grain ─────────────────
gb_cols <- grep("^GB", names(read.csv(here("data", "grammar", "GRAMBANKdf_full.csv"), nrows = 1)),
                value = TRUE)
df <- read.csv(here("data", "grammar", "GRAMBANKdf_full.csv")) |>
  filter(Language_Type == "Philippine Language") |>
  dplyr::select(language, longitude = Longitude, latitude = Latitude, all_of(gb_cols)) |>
  filter(language %in% rownames(E_sel)) |>
  arrange(match(language, rownames(E_sel)))

stopifnot(
  "Binary-matrix language set disagrees with E_sel's." =
    identical(df$language, rownames(E_sel))
)
N <- nrow(df)
message(N, " languages, ", length(gb_cols), " Grambank columns.")


# ── 3. Per-column residualization ────────────────────────────────────────────
X <- as.matrix(df[, gb_cols])
X <- X[, apply(X, 2, var) > 0]   # drop invariant columns (no signal to scrub or rotate)
# k = 0 (lm() rejects a zero-column matrix term) means nothing to scrub.
R <- if (ncol(E_sel) == 0) X else apply(X, 2, \(col) resid(lm(col ~ E_sel)))
message(ncol(X), " of ", length(gb_cols), " Grambank columns retained (non-invariant); ",
        "residualized against the ", k, "-eigenvector set above.")


# ── 4. Spatial weight matrix: fixed geographic threshold ────────────────────
GRAMMAR_dist_matrix <- read.csv(here("data", "network_distance", "GRAMMAR_dist_matrix.csv"),
                                row.names = 1, check.names = FALSE) |>
  as.matrix()
Dgeo <- GRAMMAR_dist_matrix[df$language, df$language]
diag(Dgeo) <- 0

# Geography-only tau (longest MST edge; see R/shared/spatial_threshold.R). It
# guarantees every unit a neighbour, which spca()'s matWeight normalization
# (prop.table) requires — an isolated row would divide by zero.
tau    <- mst_threshold(Dgeo)
W_best <- threshold_weights(Dgeo, tau)
stopifnot("A unit has no neighbour within tau." = all(rowSums(W_best) > 0))

pca_R <- dudi.pca(as.data.frame(R), center = TRUE, scale = FALSE, scannf = FALSE)
ms    <- multispati(pca_R, mat2listw(W_best, style = "W", zero.policy = TRUE),
                    scannf = FALSE, nfposi = 1, nfnega = 0)
best  <- list(threshold = tau, eig1 = ms$eig[1])

message(sprintf("Spatial weights: threshold = %.1f km, leading eigenvalue = %.4f.",
                best$threshold, best$eig1))


# ── 5. Real spca() ────────────────────────────────────────────────────────

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
write.csv(scores_df, file = here("data", "spca", "GRAMMAR_sPCA_scores.csv"), row.names = FALSE)

loadings_df <- tibble(feature = colnames(R), loading = spca_fit$c1[, 1]) |>
  arrange(desc(abs(loading)))
write.csv(loadings_df, file = here("data", "spca", "GRAMMAR_sPCA_loadings.csv"), row.names = FALSE)


# ── 7. Results table ─────────────────────────────────────────────────────────
GRAMMAR_sPCA_results <- tibble(
  n = N, n_evec_phylo = k, n_features_retained = ncol(R),
  threshold_km = best$threshold, eigenvalue = best$eig1,
  variance_explained = var_explained, perm_p = perm_p, n_perm = N_PERM
)
print(GRAMMAR_sPCA_results)
write.csv(GRAMMAR_sPCA_results, file = here("data", "spca", "GRAMMAR_sPCA_results.csv"), row.names = FALSE)


# ── 8. Plot: sPCA point-symbol map (size = |sPC1|, colour = sign) ───────────
dir.create(here("figures", "regression"), recursive = TRUE, showWarnings = FALSE)

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
  labs(title = "Grammar sPCA",
       subtitle = sprintf("p = %.3f, r^2 = %.3f", perm_p, var_explained))
print(p_surface)

ggsave(here("figures", "regression", "grammar_sPCA_surface.png"),
       p_surface, width = 7.5, height = 6, units = "in", dpi = 300)
saveRDS(p_surface, file = here("data", "spca", "base_plot_grammar_sPCA.rds"))
