# =============================================================================
# [4.1] Genetic Analysis — PVR-sPCA: population-structure-scrubbed spatial
#       structure
#
# Second-stage companion to [4]_GENETIC_PVR.R. Regresses population-structure-
# scrubbed Spanish admixture on Moran eigenvector maps to recover its dominant
# spatial component (sPC1), plotted as sPCA point symbols. Uses population-
# structure PCA eigenvectors (data/pca_results_phil_only.eigenvec) in place of
# phylogenetic eigenvectors — see [4]_GENETIC_PVR.R's header for why.
#
# Run order: no [0]_Phylogenetic_Tree.R dependency (unlike phoneme/grammar/
#            cognate) — reads data/pca_results_phil_only.eigenvec directly,
#            same as [4]_GENETIC_PVR.R.
#
# Input:   data/GENETIC_final.csv, data/GENETIC_dist_matrix.csv (both from [3]),
#          data/pca_results_phil_only.eigenvec
# Outputs: data/GENETIC_pvr_spca_results.csv, data/GENETIC_spca_scores.csv,
#          figures/genetic/regression/genetic_pvr_spca_surface.png,
#          data/base_plot_genetic_PVR_sPCA.rds
# =============================================================================

library(tidyverse)
library(here)
library(spdep)
library(adespatial)
library(maps)

source(here("R", "shared", "select_moran_eigenvectors.R"))

PREDICTOR <- "geodist_H1_span"
N_PERM    <- 999   # + observed arrangement = 1000 draws

GENETIC_final <- read.csv(here("data", "GENETIC_final.csv"))


# ── 1. Population-structure eigenvectors (same construction as [4]) ─────────
eigenvec_raw <- read.table(here("data", "pca_results_phil_only.eigenvec"),
                           header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eigenvec_raw) - 2
names(eigenvec_raw) <- c("population", "sample", paste0("PC", seq_len(n_pc)))

# Same population-label reconciliation [4]/[5] apply to the PLINK .mdist.id table.
eigenvec_raw$population <- ifelse(eigenvec_raw$population == "ManoboRajahKabunsuwan",
                                  "ManoboRK", eigenvec_raw$population)

eigenvec_pop <- eigenvec_raw |>
  group_by(population) |>
  summarise(across(starts_with("PC"), mean), .groups = "drop")

stopifnot(
  "eigenvec population set != GENETIC_final population set" =
    setequal(eigenvec_pop$population, GENETIC_final$population)
)


# ── 2. Analysis frame ────────────────────────────────────────────────────────
# Matched by name (no phylo object, so no position-match needed). qpAdm source
# populations carry no estimate and drop out, same as [4].
df <- GENETIC_final |>
  dplyr::select(population, span_admx, span_w, x_km = all_of(PREDICTOR),
                longitude, latitude) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  mutate(y = span_admx / max(span_admx)) |>
  as.data.frame()

y <- df$y
x <- df$x_km / 1000
E <- as.matrix(df[, paste0("PC", seq_len(n_pc)), drop = FALSE])
N <- nrow(df)

message(N, " populations with both an admixture estimate and PCA coverage.")


# ── 3. Eigenvector selection (same Moran's-I stopping rule as [4]) ──────────
E_sel <- select_moran_eigenvectors(y, x, E)
k     <- ncol(E_sel)
message(k, " population-structure PC", if (k == 1) "" else "s",
        " selected (matches [4]_GENETIC_PVR.R).")

# lm() rejects a zero-column matrix term, so k = 0 drops E_sel from the formula.
fit_e_only <- function(y, E) if (ncol(E) == 0) lm(y ~ 1) else lm(y ~ E)


# ── 4. Population-structure-scrubbed trait — no geographic covariate ────────
resid_y <- resid(fit_e_only(y, E_sel))


# ── 5. Spatial weight matrix: threshold search ──────────────────────────────
GENETIC_dist_matrix <- read.csv(here("data", "GENETIC_dist_matrix.csv"),
                                row.names = 1, check.names = FALSE) |>
  as.matrix()
Dgeo <- GENETIC_dist_matrix[df$population, df$population]
diag(Dgeo) <- 0

# cost^-2 weights, row-standardized; threshold minimizing |Moran's I| wins.
dvals <- Dgeo[upper.tri(Dgeo)]
dvals <- dvals[dvals > 0]
thresholds <- sort(unique(quantile(dvals, probs = seq(0.1, 1, length.out = 20), na.rm = TRUE)))

best <- list(threshold = NA_real_, moran_I = Inf, listw = NULL)
for (th in thresholds) {
  W <- 1 / Dgeo^2
  W[!is.finite(W)] <- 0   # zero-distance pairs -> 0, not Inf
  W[Dgeo > th] <- 0
  diag(W) <- 0
  if (all(rowSums(W) == 0)) next   # degenerate: no candidate has any neighbor

  lw <- mat2listw(W, style = "W", zero.policy = TRUE)
  mt <- tryCatch(moran.test(resid_y, lw, zero.policy = TRUE), error = function(e) NULL)
  if (is.null(mt)) next

  moran_I <- unname(mt$estimate[["Moran I statistic"]])
  if (abs(moran_I) < abs(best$moran_I)) {
    best <- list(threshold = th, moran_I = moran_I, listw = lw)
  }
}
stopifnot("No candidate threshold produced a usable spatial weights object." =
            !is.null(best$listw))

# Permutation test on the chosen weights: shuffles resid_y N_PERM times.
perm <- moran.mc(resid_y, best$listw, nsim = N_PERM, zero.policy = TRUE)
best$perm_p <- perm$p.value

message(sprintf("Spatial weights: threshold = %.1f km, Moran's I = %.4f, permutation p = %.4f.",
                best$threshold, best$moran_I, best$perm_p))


# ── 6. Spatial component: Moran Eigenvector Maps + regression ──────────────
# adegenet::spca() needs >=2 variables to rotate; mem() + regression is the
# univariate-safe equivalent. Fitted values = sPC1.
mems <- mem(best$listw, MEM.autocor = "positive")

mem_p <- vapply(as.data.frame(mems), function(v) {
  summary(lm(resid_y ~ v))$coefficients[2, 4]
}, numeric(1))
mem_sel_idx <- which(mem_p < 0.05)
if (length(mem_sel_idx) == 0) {
  mem_sel_idx <- which.max(abs(cor(resid_y, as.matrix(mems))))
  message("No MEM cleared p < 0.05 individually; falling back to the single ",
          "best-correlated MEM (", colnames(mems)[mem_sel_idx], ").")
}

MEM_sel  <- as.matrix(mems)[, mem_sel_idx, drop = FALSE]
spca_fit <- lm(resid_y ~ MEM_sel)
sPC1     <- fitted(spca_fit)

message(sprintf("sPC1: %d MEM%s selected, R^2 = %.3f.",
                length(mem_sel_idx), if (length(mem_sel_idx) == 1) "" else "s",
                summary(spca_fit)$r.squared))


# ── 7. Per-population sPC1 scores ────────────────────────────────────────────
scores_df <- df |>
  mutate(sPC1 = sPC1) |>
  summarise(sPC1 = mean(sPC1), .by = c(population, longitude, latitude))

write.csv(scores_df, file = here("data", "GENETIC_spca_scores.csv"), row.names = FALSE)


# ── 8. Results table ─────────────────────────────────────────────────────────
GENETIC_pvr_spca_results <- tibble(
  n = N, n_evec_popstruct = k,
  threshold_km = best$threshold, moran_I = best$moran_I, moran_perm_p = best$perm_p,
  n_perm = N_PERM, n_mem_selected = length(mem_sel_idx), spca_r2 = summary(spca_fit)$r.squared
)
print(GENETIC_pvr_spca_results)
write.csv(GENETIC_pvr_spca_results,
          file = here("data", "GENETIC_pvr_spca_results.csv"), row.names = FALSE)


# ── 9. Plot: sPCA point-symbol map (size = |sPC1|, colour = sign) ───────────
dir.create(here("figures", "genetic", "regression"), recursive = TRUE, showWarnings = FALSE)

world_map  <- map_data("world")
map_subset <- world_map |> filter(region %in% c("Philippines", "Malaysia"))

# Legend: 4 fixed keys (median |sPC1| split x sign), merged via shared factor.
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
  labs(title = "Genetic PVR-sPCA: population-structure-scrubbed spatial structure",
       subtitle = sprintf("Moran's I permutation p = %.3f (%d permutations) | eigenvector R^2 = %.3f",
                          best$perm_p, N_PERM + 1, summary(spca_fit)$r.squared))
print(p_surface)

ggsave(here("figures", "genetic", "regression", "genetic_pvr_spca_surface.png"),
       p_surface, width = 7.5, height = 6, units = "in", dpi = 300)
saveRDS(p_surface, file = here("data", "base_plot_genetic_PVR_sPCA.rds"))
