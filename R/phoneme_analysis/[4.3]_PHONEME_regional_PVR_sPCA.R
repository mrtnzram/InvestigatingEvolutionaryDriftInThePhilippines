# =============================================================================
# [4.3] Phoneme Analysis — regional PVR-MEM
#
# Re-runs [4.1]'s analysis independently within Luzon, Visayas and Mindanao, to
# test whether localized spatial structure is significant where the global fit
# is not. Everything is re-derived per region — the tree is re-pruned and the
# phylogenetic eigenvectors recomputed on it, not sliced from the global fit.
#
# Regions below MIN_N are skipped, not run with a reduced permutation count:
# spdep::moran.mc(nsim = 999) errors when 999 > N!, i.e. below N = 7.
#
# Run order: requires `tree_pruned` and `tree_df_matched` from [0]_Phylogenetic_Tree.R.
#
# Input:   data/PHONEME_final.csv, data/PHONEME_dist_matrix.csv (both from [3]),
#          data/PHONEME_subgroup_lookup.csv (region column, from [0])
# Outputs: data/PHONEME_regional_pvr_mem_results.csv,
#          data/PHONEME_regional_pvr_mem_scores.csv,
#          figures/phoneme/regression/phoneme_regional_pvr_mem_<region>.png
#          (one standalone file per region — not combined into a single figure,
#          since each region's map has its own aspect ratio)
# =============================================================================

library(PVR)
library(ape)
library(tidyverse)
library(here)
library(spdep)
library(adespatial)
library(maps)

stopifnot(
  "Run [0]_Phylogenetic_Tree.R first: `tree_pruned` is not defined." =
    exists("tree_pruned"),
  "Run [0]_Phylogenetic_Tree.R first: `tree_df_matched` is not defined." =
    exists("tree_df_matched")
)

PREDICTOR <- "geodist_H1_span"
N_PERM    <- 999
MIN_N     <- 7      # moran.mc(nsim = 999) needs N! > 999
REGIONS   <- c("Luzon", "Visayas", "Mindanao")

PHONEME_final  <- read.csv(here("data", "phoneme", "network_distance", "PHONEME_final.csv"))
region_lookup  <- read.csv(here("data", "phoneme", "initial_datasets", "PHONEME_subgroup_lookup.csv")) |>
  dplyr::select(language, region)
tip_map        <- tree_df_matched |> dplyr::select(original, ph)
dist_matrix    <- read.csv(here("data", "phoneme", "network_distance", "PHONEME_dist_matrix.csv"),
                           row.names = 1, check.names = FALSE) |> as.matrix()


# ── 1. Tip-grain frame carrying region (same construction as [4.1]) ─────────
df_all <- PHONEME_final |>
  dplyr::select(language, y = cossim_span_norm, x_km = all_of(PREDICTOR),
                longitude, latitude) |>
  left_join(tip_map, by = c("language" = "ph")) |>
  left_join(region_lookup, by = "language") |>
  filter(!is.na(original), !is.na(y), !is.na(x_km)) |>
  as.data.frame()
df_all <- df_all[match(tree_pruned$tip.label, df_all$original), ]

stopifnot(
  "Analysis frame and pruned tree disagree on the tip set." =
    nrow(df_all) == length(tree_pruned$tip.label),
  "Analysis frame is not in tip.label order — PVR matches by position." =
    !anyNA(df_all$original) && identical(df_all$original, tree_pruned$tip.label),
  "Some tips have no region — rerun [0]_Phylogenetic_Tree.R." = !anyNA(df_all$region)
)


# ── 2. Per-region fit ────────────────────────────────────────────────────────
results <- list()
scores  <- list()

for (rg in REGIONS) {
  df <- df_all[df_all$region == rg, ]
  N  <- nrow(df)
  row <- tibble(region = rg, n = N, n_langs = n_distinct(df$language),
                status = NA_character_, n_evec_phylo = NA_integer_,
                threshold_km = NA_real_, moran_I = NA_real_,
                moran_perm_p = NA_real_, n_perm = N_PERM,
                n_mem_selected = NA_integer_, mem_selection = NA_character_,
                mem_r2 = NA_real_)

  if (N < MIN_N) {
    row$status <- "skipped_small_n"
    message(sprintf("%-9s N=%2d — skipped (need >= %d for a %d-permutation test).",
                    rg, N, MIN_N, N_PERM))
    results[[rg]] <- row; next
  }

  # Re-prune. drop.tip collapses singleton paths and merges branch lengths, so
  # the regional patristic distances are genuinely re-derived, not a submatrix.
  tree_rg <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, df$original))
  df <- df[match(tree_rg$tip.label, df$original), ]
  stopifnot(
    "Regional frame is not in regional tip.label order." =
      identical(df$original, tree_rg$tip.label),
    "Regional subtree is not ultrametric." = is.ultrametric(tree_rg)
  )

  y <- df$y; x <- df$x_km / 1000

  pvr_dec <- PVRdecomp(tree_rg, scale = TRUE)
  pvr_fit <- PVR(pvr_dec, phy = tree_rg, trait = y, envVar = x, method = "moran")
  E_sel   <- as.matrix(pvr_fit@Selection$Vectors)
  k       <- ncol(E_sel)
  row$n_evec_phylo <- k

  # A saturated design (k = N-1) makes resid_y identically zero, which would
  # poison every statistic below without erroring.
  if (k >= N - 1) {
    row$status <- "saturated_evec"
    message(sprintf("%-9s N=%2d — skipped (PVR selected k=%d >= N-1; fit saturated).",
                    rg, N, k))
    results[[rg]] <- row; next
  }

  resid_y <- resid(lm(y ~ E_sel))

  # Spatial weights: same quantile grid / cost^-2 / row-standardized search as
  # [4.1], minimizing |Moran's I|.
  Dgeo <- dist_matrix[df$language, df$language]
  dimnames(Dgeo) <- list(df$original, df$original)
  diag(Dgeo) <- 0
  dvals <- Dgeo[upper.tri(Dgeo)]; dvals <- dvals[dvals > 0]
  thresholds <- sort(unique(quantile(dvals, probs = seq(0.1, 1, length.out = 20),
                                     na.rm = TRUE)))

  best <- list(threshold = NA_real_, moran_I = Inf, listw = NULL)
  for (th in thresholds) {
    W <- 1 / Dgeo^2
    W[!is.finite(W)] <- 0
    W[Dgeo > th] <- 0
    diag(W) <- 0
    if (all(rowSums(W) == 0)) next
    lw <- mat2listw(W, style = "W", zero.policy = TRUE)
    mt <- tryCatch(moran.test(resid_y, lw, zero.policy = TRUE), error = function(e) NULL)
    if (is.null(mt)) next
    mI <- unname(mt$estimate[["Moran I statistic"]])
    if (abs(mI) < abs(best$moran_I)) best <- list(threshold = th, moran_I = mI, listw = lw)
  }
  if (is.null(best$listw)) {
    row$status <- "no_usable_threshold"
    message(sprintf("%-9s N=%2d — skipped (no candidate threshold gave usable weights).", rg, N))
    results[[rg]] <- row; next
  }

  perm <- moran.mc(resid_y, best$listw, nsim = N_PERM, zero.policy = TRUE)

  mems  <- mem(best$listw, MEM.autocor = "positive")
  mem_p <- vapply(as.data.frame(mems),
                  function(v) summary(lm(resid_y ~ v))$coefficients[2, 4], numeric(1))
  sel   <- which(mem_p < 0.05)
  # Record which branch was taken: a single max-correlation MEM at small N is an
  # overfit, and that should be visible in the results table.
  branch <- if (length(sel) == 0) "fallback_max_cor" else "p_filter"
  if (length(sel) == 0) sel <- which.max(abs(cor(resid_y, as.matrix(mems))))

  MEM_sel  <- as.matrix(mems)[, sel, drop = FALSE]
  mem_fit  <- lm(resid_y ~ MEM_sel)
  sPC1     <- fitted(mem_fit)

  row$status         <- "ok"
  row$threshold_km   <- best$threshold
  row$moran_I        <- best$moran_I
  row$moran_perm_p   <- perm$p.value
  row$n_mem_selected <- length(sel)
  row$mem_selection  <- branch
  row$mem_r2         <- summary(mem_fit)$r.squared
  results[[rg]] <- row

  # Tip -> language, matching [4.1]'s convention.
  scores[[rg]] <- df |>
    mutate(sPC1 = sPC1) |>
    summarise(sPC1 = mean(sPC1), .by = c(language, longitude, latitude)) |>
    mutate(region = rg)

  message(sprintf("%-9s N=%2d | k=%d | th=%.1f km | I=%.4f | p=%.4f | %d MEM (%s) | R2=%.3f",
                  rg, N, k, best$threshold, best$moran_I, perm$p.value,
                  length(sel), branch, summary(mem_fit)$r.squared))
}

PHONEME_regional_results <- bind_rows(results)
print(PHONEME_regional_results)
write.csv(PHONEME_regional_results,
          here("data", "phoneme", "PVR", "PHONEME_regional_pvr_mem_results.csv"), row.names = FALSE)

scores_all <- bind_rows(scores)
write.csv(scores_all,
          here("data", "phoneme", "PVR", "PHONEME_regional_pvr_mem_scores.csv"), row.names = FALSE)
stopifnot("No region produced scores — nothing to plot." = nrow(scores_all) > 0)


# ── 3. Plot: one panel per region, legend normalized across all of them ─────
dir.create(here("figures", "phoneme", "regression"), recursive = TRUE, showWarnings = FALSE)

map_subset <- map_data("world") |> filter(region %in% c("Philippines", "Malaysia"))

# ONE median split over the pooled scores, so square sizes are comparable
# between panels. Each region's scores come from an independently-estimated
# model, so this is a visual aid, not a claim the axes share a scale.
mag_threshold <- median(abs(scores_all$sPC1))

key_sizes  <- c("Large Positive" = 10, "Small Positive" = 4,
                "Small Negative" = 4,  "Large Negative" = 10)
key_fills  <- c("Large Positive" = "black", "Small Positive" = "black",
                "Small Negative" = "white", "Large Negative" = "white")
key_labels <- c("Large Positive" = "Large (+)", "Small Positive" = "Small (+)",
                "Small Negative" = "Small (-)", "Large Negative" = "Large (-)")

PAD <- 0.5
panel <- function(rg) {
  s <- scores_all |>
    filter(region == rg) |>
    mutate(sign = if_else(sPC1 >= 0, "Positive", "Negative"),
           magnitude = if_else(abs(sPC1) >= mag_threshold, "Large", "Small"),
           key = factor(paste(magnitude, sign),
                        levels = c("Large Positive", "Small Positive",
                                   "Small Negative", "Large Negative")))
  st <- PHONEME_regional_results |> filter(region == rg)
  ggplot() +
    geom_polygon(data = map_subset, aes(x = long, y = lat, group = group),
                 fill = "gray97", color = "black") +
    geom_point(data = s, aes(x = longitude, y = latitude, size = key, fill = key),
               shape = 22, colour = "black", stroke = 0.6) +
    scale_size_manual(values = key_sizes, labels = key_labels, name = "sPC1", drop = FALSE) +
    scale_fill_manual(values = key_fills, labels = key_labels, name = "sPC1", drop = FALSE) +
    coord_fixed(xlim = range(s$longitude) + c(-PAD, PAD),
                ylim = range(s$latitude)  + c(-PAD, PAD)) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.title = element_blank(), axis.ticks = element_blank()) +
    labs(title = rg,
         subtitle = sprintf("N=%d | p = %.3f | R^2 = %.3f",
                            st$n, st$moran_perm_p, st$mem_r2))
}

# One standalone file per region (not a combined multi-panel figure — each
# region's coord_fixed() box has its own aspect ratio, e.g. Mindanao's is much
# wider than Luzon's, so cramming them into one row via patchwork stretched
# panels unevenly). mag_threshold/key_* above are still shared across regions,
# so the legend in every file reads on the same scale even though the files
# are separate.
run_regions <- PHONEME_regional_results |> filter(status == "ok") |> pull(region)
panels <- lapply(run_regions, panel)
names(panels) <- run_regions

for (rg in run_regions) {
  print(panels[[rg]])
  ggsave(here("figures", "phoneme", "regression",
              paste0("phoneme_regional_pvr_mem_", tolower(rg), ".png")),
         panels[[rg]], width = 5.5, height = 6, units = "in", dpi = 300)
}

skipped <- PHONEME_regional_results |> filter(status != "ok")
if (nrow(skipped) > 0) {
  message("Not run: ", paste(sprintf("%s (N=%d, %s)", skipped$region, skipped$n, skipped$status),
                              collapse = "; "))
}
