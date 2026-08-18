#!/usr/bin/env Rscript
# =============================================================================
# [4.2] Genetic Analysis — true multivariate sPCA (REMOTE — run on the server
# where phil_only_pruned lives; see README_GENETIC_sPCA_remote.md)
#
# Genetic analogue of [4.2]_PHONEME_sPCA.R/[4.2]_GRAMMAR_sPCA.R: residualizes
# a genotype matrix (aggregated to population allele frequency via bigsnpr/
# bigstatsr's file-backed matrix) against the same population-structure PCs
# [4]_GENETIC_PVR.R uses, then runs real adegenet::spca(). Population is read
# from the .fam FID directly (matches pca_results_phil_only.eigenvec's
# convention, NOT genetic_feems.py's PHIL_<population>_<code> scheme) —
# --id-prefix/--alias are a fallback if phil_only_pruned differs.
#
# Usage:
#   Rscript "[4.2]_GENETIC_sPCA_remote.R" --bfile /path/to/phil_only_pruned \
#     --eigenvec pca_results_phil_only.eigenvec --genetic-final GENETIC_final.csv \
#     --dist-matrix GENETIC_dist_matrix.csv --out ./out
#
# Input:   PLINK --bfile, pca_results_phil_only.eigenvec, GENETIC_final.csv,
#          GENETIC_dist_matrix.csv (files to scp up: see README)
# Outputs: GENETIC_sPCA_results.csv, GENETIC_sPCA_scores.csv,
#          GENETIC_sPCA_loadings.csv (written to --out; scp back into data/)
# =============================================================================

suppressPackageStartupMessages({
  library(bigsnpr)
  library(bigstatsr)
  library(spdep)
  library(adespatial)
  library(adegenet)
  library(RSpectra)
  library(dplyr)
  library(tibble)
})

N_PERM <- 999   # + observed arrangement = 1000 draws

# ── CLI args (base commandArgs — no new parsing dependency) ─────────────────
parse_args <- function(argv) {
  defaults <- list(bfile = NULL, eigenvec = NULL, `genetic-final` = NULL,
                   `dist-matrix` = NULL, out = "./out",
                   `id-prefix` = "", alias = character(0), `max-snps` = 0)
  i <- 1
  while (i <= length(argv)) {
    flag <- sub("^--", "", argv[i])
    if (!flag %in% names(defaults)) stop("unknown flag: --", flag)
    if (flag == "alias") {
      defaults$alias <- c(defaults$alias, argv[i + 1])
    } else {
      defaults[[flag]] <- argv[i + 1]
    }
    i <- i + 2
  }
  defaults$`max-snps` <- as.integer(defaults$`max-snps`)
  required <- c("bfile", "eigenvec", "genetic-final", "dist-matrix")
  missing <- required[vapply(defaults[required], is.null, logical(1))]
  if (length(missing) > 0) stop("missing required flag(s): --", paste(missing, collapse = ", --"))
  defaults
}
args <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)


# ── 1. Population aliases — same reconciliation [4]_GENETIC_PVR.R/genetic_feems.py use ──
DEFAULT_ALIASES <- c(manoborajahkabunsuwan = "ManoboRK")
aliases <- DEFAULT_ALIASES
for (spec in args$alias) {
  kv <- strsplit(spec, "=", fixed = TRUE)[[1]]
  aliases[tolower(kv[1])] <- kv[2]
}

resolve_population <- function(fid, id_prefix, aliases, known_pops) {
  stripped <- if (nzchar(id_prefix) && startsWith(tolower(fid), tolower(id_prefix))) {
    substring(fid, nchar(id_prefix) + 1)
  } else fid
  candidates <- unique(c(stripped, strsplit(stripped, "_", fixed = TRUE)[[1]][1]))
  for (cand in candidates) {
    key <- tolower(cand)
    if (key %in% tolower(known_pops)) return(known_pops[tolower(known_pops) == key][1])
    if (key %in% names(aliases)) {
      aliased <- aliases[[key]]
      if (tolower(aliased) %in% tolower(known_pops)) return(aliased)
    }
  }
  NA_character_
}


# ── 2. Efficient load: bigsnpr FBM, not a dense in-memory matrix ────────────
message("Loading ", args$bfile, ".bed via bigsnpr ...")
rds_path <- paste0(args$bfile, ".rds")
rds <- if (file.exists(rds_path)) rds_path else snp_readBed(paste0(args$bfile, ".bed"))
obj <- snp_attach(rds)
G   <- obj$genotypes
message(nrow(G), " individuals x ", ncol(G), " SNPs (post-pruning).")

GENETIC_final <- read.csv(args$`genetic-final`)
population <- vapply(obj$fam$family.ID, resolve_population, character(1),
                     id_prefix = args$`id-prefix`, aliases = aliases,
                     known_pops = GENETIC_final$population)
n_unmatched <- sum(is.na(population))
if (n_unmatched > 0) {
  message(n_unmatched, " of ", length(population),
          " individuals had no population match — dropped. Check --id-prefix/--alias",
          " if this is more than a handful.")
}
keep <- !is.na(population)


# ── 3. Population-level allele frequency, chunked over the FBM ─────────────
# The only step that needs bigsnpr's efficient loading — everything after
# this operates on a plain 115 x p_snps dense matrix.
pop_idx <- split(which(keep), population[keep])
message("Aggregating to ", length(pop_idx), " populations ...")

freq <- big_apply(G, a.FUN = function(X, ind, pop_idx) {
  block <- X[, ind, drop = FALSE]
  sapply(pop_idx, function(rows) colMeans(block[rows, , drop = FALSE], na.rm = TRUE) / 2)
}, a.combine = "rbind", pop_idx = pop_idx, block.size = 2000)
freq <- t(freq)
rownames(freq) <- names(pop_idx)
colnames(freq) <- obj$map$marker.ID   # SNP identity, else lost by the loadings output

if (args$`max-snps` > 0 && ncol(freq) > args$`max-snps`) {
  set.seed(1)
  freq <- freq[, sample(ncol(freq), args$`max-snps`)]
  message("Subsampled to ", ncol(freq), " SNPs (--max-snps).")
}

freq <- freq[, apply(freq, 2, var, na.rm = TRUE) > 0]   # drop invariant columns
freq[is.nan(freq)] <- NA
freq <- freq[, colSums(is.na(freq)) == 0]                # drop columns any population missed entirely
message(ncol(freq), " SNPs retained (non-invariant, fully observed).")


# ── 4. Population-structure eigenvectors — same E_sel as [4]_GENETIC_PVR.R ──
# select_moran_eigenvectors.R is scp'd flat next to this script (see README).
script_dir <- {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1) dirname(sub("^--file=", "", file_arg)) else "."
}
source(file.path(script_dir, "select_moran_eigenvectors.R"))

eigenvec_raw <- read.table(args$eigenvec, header = FALSE, stringsAsFactors = FALSE)
n_pc <- ncol(eigenvec_raw) - 2
names(eigenvec_raw) <- c("population", "sample", paste0("PC", seq_len(n_pc)))
eigenvec_raw$population <- ifelse(eigenvec_raw$population == "ManoboRajahKabunsuwan",
                                  "ManoboRK", eigenvec_raw$population)
eigenvec_pop <- eigenvec_raw |>
  group_by(population) |>
  summarise(across(starts_with("PC"), mean), .groups = "drop")

PREDICTOR <- "geodist_H1_span"
gf <- GENETIC_final |>
  dplyr::select(population, span_admx, x_km = all_of(PREDICTOR)) |>
  inner_join(eigenvec_pop, by = "population") |>
  filter(!is.na(span_admx), !is.na(x_km)) |>
  mutate(y = span_admx / max(span_admx)) |>
  as.data.frame()

y_pc <- gf$y
x_pc <- gf$x_km / 1000
E_pc <- as.matrix(gf[, paste0("PC", seq_len(n_pc)), drop = FALSE])
E_sel <- select_moran_eigenvectors(y_pc, x_pc, E_pc)
k <- ncol(E_sel)
rownames(E_sel) <- gf$population
message(k, " population-structure PC", if (k == 1) "" else "s",
        " selected (matches [4]_GENETIC_PVR.R).")


# ── 5. Align + residualize (hat-matrix projection, cheap at N=115) ─────────
common <- intersect(rownames(freq), rownames(E_sel))
freq  <- freq[common, ]
E_sel <- E_sel[common, , drop = FALSE]
N <- length(common)
message(N, " populations with both genotype and population-structure coverage.")

# k = 0 (no PC clears select_moran_eigenvectors()'s bar) means nothing to
# project out — solve() on a 0-column crossproduct errors, so short-circuit.
if (ncol(E_sel) == 0) {
  R <- freq
} else {
  H <- E_sel %*% solve(t(E_sel) %*% E_sel) %*% t(E_sel)
  R <- freq - H %*% freq
}


# ── 6. Reduce to <=N-1 PC scores before any multispati/spca call ────────────
# adespatial::multispati() builds a p x p covariance matrix and calls eigen()
# on it regardless of how many axes are requested — cubic in SNP count
# (measured: 71s for a single call at p=4542; ~9 days extrapolated at
# p=100,000). Centering an N-row matrix caps its rank at N-1, so re-expressing
# R in "sample space" via truncated SVD first is lossless (verified: identical
# eigenvalues/row-scores up to an arbitrary sign flip) and turns every
# downstream multispati/spca call into an r x r problem instead (measured:
# 0.002s at the same p=4542 — 35,000x faster). Loadings come back in
# score-space and are projected back to per-SNP space via V at the end.
Rc <- scale(R, center = TRUE, scale = FALSE)
r  <- min(nrow(Rc) - 1, ncol(Rc))
sv <- svds(Rc, k = r)
scores <- sv$u %*% diag(sv$d, r, r)
V <- sv$v
rownames(scores) <- rownames(R)
message("Reduced ", ncol(R), " SNPs to ", r, " PC scores (lossless — full rank of the centered residual).")


# ── 7. Spatial weight matrix: threshold search, maximizing sPCA eigenvalue ──
GENETIC_dist_matrix <- read.csv(args$`dist-matrix`, row.names = 1, check.names = FALSE) |>
  as.matrix()
Dgeo <- GENETIC_dist_matrix[common, common]
diag(Dgeo) <- 0

dvals <- Dgeo[upper.tri(Dgeo)]
dvals <- dvals[dvals > 0]
thresholds <- sort(unique(quantile(dvals, probs = seq(0.1, 1, length.out = 20), na.rm = TRUE)))

pca_scores <- dudi.pca(as.data.frame(scores), center = FALSE, scale = FALSE, scannf = FALSE)

best <- list(threshold = NA_real_, eig1 = -Inf, listw = NULL)
for (th in thresholds) {
  W <- 1 / Dgeo^2
  W[!is.finite(W)] <- 0
  W[Dgeo > th] <- 0
  diag(W) <- 0
  if (any(rowSums(W) == 0)) next   # spca()'s matWeight normalization needs no isolates

  lw <- mat2listw(W, style = "W", zero.policy = TRUE)
  ms <- tryCatch(multispati(pca_scores, lw, scannf = FALSE, nfposi = 1, nfnega = 0),
                 error = function(e) NULL)
  if (is.null(ms)) next
  if (ms$eig[1] > best$eig1) best <- list(threshold = th, eig1 = ms$eig[1], listw = lw)
}
stopifnot("No candidate threshold produced a usable spatial weights object." =
            !is.null(best$listw))
message(sprintf("Spatial weights: threshold = %.1f km, leading eigenvalue = %.4f.",
                best$threshold, best$eig1))


# ── 8. Real spca() ────────────────────────────────────────────────────────
coords <- GENETIC_final[match(common, GENETIC_final$population), c("longitude", "latitude")]

W_best <- 1 / Dgeo^2
W_best[!is.finite(W_best)] <- 0
W_best[Dgeo > best$threshold] <- 0
diag(W_best) <- 0

spca_fit <- spca(as.data.frame(scores), xy = as.matrix(coords),
                 matWeight = W_best, scannf = FALSE, nfposi = 1, nfnega = 1)
sPC1 <- spca_fit$li[, 1]

perm <- spca_randtest(spca_fit, nperm = N_PERM)
# NOT perm$global$pvalue — that field is broken upstream and always returns 1:
# spca_randtest passes `sim = sum(sims[sims >= 0])`, a scalar summed over every
# permutation, so obs can never exceed it. eigentest uses a proper per-permutation
# vector (`sims[e, ]`); [1, ] is the leading axis.
perm_p <- perm$eigentest[1, "sim_p"]
var_explained <- best$eig1 / sum(abs(pca_scores$eig))
message(sprintf("sPCA: variance explained (axis 1) = %.3f, permutation p = %.4f.",
                var_explained, perm_p))


# ── 9. Write outputs ──────────────────────────────────────────────────────
scores_df <- tibble(population = common, longitude = coords$longitude,
                    latitude = coords$latitude, sPC1 = sPC1)
write.csv(scores_df, file.path(args$out, "GENETIC_sPCA_scores.csv"), row.names = FALSE)

# Per-SNP loading = V %*% (loading in score-space) — exact, not approximate,
# because `scores = Rc %*% V` is itself lossless (see §6).
loadings_df <- tibble(snp = colnames(R), loading = as.numeric(V %*% as.matrix(spca_fit$c1)[, 1])) |>
  arrange(desc(abs(loading)))
write.csv(loadings_df, file.path(args$out, "GENETIC_sPCA_loadings.csv"), row.names = FALSE)

results_df <- tibble(
  n = N, n_evec_popstruct = k, n_snps_retained = ncol(R), n_pc_scores = r,
  threshold_km = best$threshold, eigenvalue = best$eig1,
  variance_explained = var_explained, perm_p = perm_p, n_perm = N_PERM
)
print(results_df)
write.csv(results_df, file.path(args$out, "GENETIC_sPCA_results.csv"), row.names = FALSE)

message("Done. scp ", args$out, "/GENETIC_sPCA_*.csv back into data/.")
