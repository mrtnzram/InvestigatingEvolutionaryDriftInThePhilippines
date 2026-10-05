# =============================================================================
# Shared — pure phylogeny on the geographic axis
#
# Fitters and figures for the [4]_*_PVR.R scripts: orient a phylogenetic axis
# with distance from Manila, residualize it on that distance, and regress the
# Spanish signal on geography and on the residual.
#
# Inputs:  none (sourced)
# Outputs: none (functions only)
# =============================================================================

library(tidyverse)
library(patchwork)

DISPLAY_UNIT_KM <- 100   # geography slopes and axes are reported per this many km


#' First principal coordinate of a pairwise distance matrix.
#'
#' @param D Square distance matrix with matching row and column names.
#' @return list(axis1 = named numeric vector, var_share = Axis-1 eigenvalue share).
pcoa_axis1 <- function(D) {
  stopifnot("Distance matrix row and column names disagree." =
              identical(rownames(D), colnames(D)))
  p <- cmdscale(as.dist(D), k = 1, eig = TRUE)
  list(axis1 = p$points[, 1], var_share = p$eig[1] / sum(p$eig[p$eig > 0]))
}


#' One results row from a single-predictor lm.
lm_row <- function(m, term, model, response, predictor, scale = 1) {
  s  <- summary(m)
  cf <- s$coefficients[term, ]
  fs <- s$fstatistic
  tibble(model, response, predictor, n = nobs(m),
         slope         = cf[["Estimate"]]   * scale,
         slope_se      = cf[["Std. Error"]] * scale,
         f.statistic   = unname(fs[1]),
         p.value       = unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE)),
         r.squared     = s$r.squared,
         adj.r.squared = s$adj.r.squared)
}


#' Decompose a phylogenetic axis against geography and regress y on both parts.
#'
#' @param frame     Data frame with `y` and `x_km`, rows aligned to `axis1`.
#' @param axis1     Phylogenetic Axis 1, one value per row of `frame`.
#' @param var_share Axis-1 variance share, carried into the results table.
#' @param predictor Name of the geography column, for the results table.
#' @param response  Name of the y column, for the results table.
#' @param axis_name Name of the phylogenetic axis, for the results table.
#' @param anc_label What Axis 1 measures, for the variance-partition labels.
#' @param zeros     Zero-handling label ("removed"/"kept"), or NULL to omit the column.
#' @return list(frame, pvr_geo_lm, vp, r2_full, rows).
pure_phylo_fit <- function(frame, axis1, var_share, predictor, response,
                           axis_name = "phylo_axis1", anc_label = "phylogeny",
                           zeros = NULL) {
  stopifnot("axis1 and frame differ in length." = length(axis1) == nrow(frame))
  x <- frame$x_km / 1000                             # fit on Mm
  axis1 <- unname(axis1)

  # Orient Axis 1 to rise with distance from Manila.
  flipped <- cor(axis1, x) < 0
  if (flipped) axis1 <- -axis1

  pvr_geo_lm <- lm(axis1 ~ x)
  r <- unname(resid(pvr_geo_lm))
  frame$phylo_axis1        <- axis1
  frame$pure_phylogeny     <- r
  frame$pure_phylogeny_pos <- r - min(r)            # 0 -> + like geography; a constant shift, so no fit changes
  stopifnot("pure_phylogeny is not orthogonal to geography." = abs(cor(r, x)) < 1e-8)

  d <- data.frame(y = frame$y, x = x, axis1 = axis1, pure = frame$pure_phylogeny_pos)
  m_geo  <- lm(y ~ x,         data = d)
  m_pure <- lm(y ~ pure,      data = d)
  m_phy  <- lm(y ~ axis1,     data = d)
  m_full <- lm(y ~ x + axis1, data = d)

  R2_geo  <- summary(m_geo)$r.squared
  R2_pure <- summary(m_pure)$r.squared
  R2_phy  <- summary(m_phy)$r.squared
  R2_full <- summary(m_full)$r.squared
  stopifnot("geography R^2 + pure-phylogeny R^2 != R^2(y ~ geography + axis 1)." =
              isTRUE(all.equal(R2_geo + R2_pure, R2_full, tolerance = 1e-8)))

  # b is a difference of R^2, not a variance, and goes negative under suppression.
  vp <- tibble(
    component = c("a", "b", "c", "d"),
    meaning   = c("geography only", "shared", paste(anc_label, "only"), "unexplained"),
    fraction  = c(R2_full - R2_phy, R2_geo + R2_phy - R2_full, R2_pure, 1 - R2_full)
  )
  stopifnot("Variance partition does not sum to 1." =
              isTRUE(all.equal(sum(vp$fraction), 1, tolerance = 1e-8)))

  geo_lab <- paste0(predictor, " (per ", DISPLAY_UNIT_KM, " km)")
  scl     <- DISPLAY_UNIT_KM / 1000
  rows <- bind_rows(
    lm_row(pvr_geo_lm, "x",    "phylogeny_vs_geography", axis_name, geo_lab, scl),
    lm_row(m_geo,      "x",    "A_y_vs_geography",       response,  geo_lab, scl),
    lm_row(m_pure,     "pure", "B_y_vs_pure_phylogeny",  response,  "pure_phylogeny_pos")
  ) |>
    mutate(axis1_var_share = var_share, axis1_flipped = flipped, .after = n)
  if (!is.null(zeros)) rows <- mutate(rows, zeros = zeros, .after = model)

  list(frame = frame, pvr_geo_lm = pvr_geo_lm, vp = vp, r2_full = R2_full, rows = rows)
}


#' Figure 1: geography and pure phylogeny profiles, languages ordered by distance from Manila.
plot_profile <- function(fit, title, anc_label) {
  d <- fit$frame |>
    arrange(x_km, pure_phylogeny_pos) |>
    mutate(rank = row_number())
  expect <- -min(fit$frame$pure_phylogeny)           # level of pure_phylogeny_pos where pure_phylogeny = 0

  p_geo <- ggplot(d, aes(rank, x_km / DISPLAY_UNIT_KM)) +
    geom_line(colour = "grey60", linewidth = 0.4) +
    geom_point(size = 1.6) +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = title, x = NULL,
         y = sprintf("Distance from Manila (%d km)", DISPLAY_UNIT_KM))

  p_pure <- ggplot(d, aes(rank, pure_phylogeny_pos)) +
    geom_hline(yintercept = expect, linetype = "dashed", colour = "grey40") +
    geom_segment(aes(xend = rank, y = expect, yend = pure_phylogeny_pos),
                 colour = "grey60", linewidth = 0.4) +
    geom_point(size = 1.6) +
    theme_bw() +
    labs(x = "Rank by distance from Manila", y = paste("Pure", anc_label))

  p_geo / p_pure
}


#' One scatter panel with an lm fit; subtitle carries that fit's p and r^2.
lm_panel <- function(d, xvar, yvar, xlab, ylab, ptitle) {
  s  <- summary(lm(d[[yvar]] ~ d[[xvar]]))
  fs <- s$fstatistic
  ggplot(d, aes(.data[[xvar]], .data[[yvar]])) +
    geom_point(size = 2, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x,
                colour = "steelblue", fill = "steelblue", alpha = 0.25, linewidth = 1.2) +
    theme_bw() +
    labs(title = ptitle,
         subtitle = sprintf("p = %.3g, r^2 = %.3f",
                            unname(pf(fs[1], fs[2], fs[3], lower.tail = FALSE)), s$r.squared),
         x = xlab, y = ylab)
}


#' Figure 2: oriented Axis 1 and pure phylogeny against geography, with lm fits.
plot_scatter <- function(fit, title, anc_label, axis_label) {
  d    <- fit$frame |> mutate(x_disp = x_km / DISPLAY_UNIT_KM)
  xlab <- sprintf("Distance from Manila (%d km)", DISPLAY_UNIT_KM)
  (lm_panel(d, "x_disp", "phylo_axis1", xlab, paste(axis_label, "(oriented)"), "Raw") |
     lm_panel(d, "x_disp", "pure_phylogeny_pos", xlab, paste("Pure", anc_label), "Pure")) +
    plot_annotation(title = title)
}


#' Figure 3: the Spanish signal against geography and against pure phylogeny, on a shared y-axis.
plot_y_scatter <- function(fit, title, anc_label, y_label) {
  d <- fit$frame |> mutate(x_disp = x_km / DISPLAY_UNIT_KM)
  (lm_panel(d, "x_disp", "y", sprintf("Distance from Manila (%d km)", DISPLAY_UNIT_KM),
            y_label, "Geography") |
     lm_panel(d, "pure_phylogeny_pos", "y", paste("Pure", anc_label), y_label,
              paste("Pure", anc_label))) +
    plot_layout(axes = "collect_y") +
    plot_annotation(title = title) &
    coord_cartesian(ylim = range(d$y))
}


#' Variance-partition diagram: total variance is the full-width line, geography
#' the bar above it (a + b) and ancestry the bar below (b + c), overlapping on b.
vp_diagram <- function(vp, title, subtitle, anc_label) {
  f   <- setNames(vp$fraction, vp$component)[c("a", "b", "c", "d")]
  x   <- cumsum(c(0, f[["a"]], f[["b"]], f[["c"]]))
  mid <- x + f / 2

  # Narrow spans cannot sit under their own label, so labels are pushed apart to
  # a minimum gap and tied back to their span by a leader.
  gap <- 0.085
  o <- order(mid); p <- mid[o]
  for (i in seq_along(p)[-1]) p[i] <- max(p[i], p[i - 1] + gap)
  if (max(p) > 1) {
    p <- p - (max(p) - 1)
    for (i in rev(seq_along(p))[-1]) p[i] <- min(p[i], p[i + 1] - gap)
  }
  lx <- numeric(length(mid)); lx[o] <- p

  box <- data.frame(xmin = c(min(x[1], x[3]), min(x[2], x[4])),
                    xmax = c(max(x[1], x[3]), max(x[2], x[4])),
                    ymin = c(0.500, 0.375), ymax = c(0.625, 0.500))
  box <- box[box$xmax > box$xmin, ]   # a zero-width bar would draw as a stray tick

  ggplot() +
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", colour = "black", linewidth = 0.4) +
    geom_segment(data = data.frame(x = x), aes(x = x, xend = x, y = 0.5, yend = 0.68),
                 linetype = "dashed", linewidth = 0.3) +
    geom_segment(data = data.frame(mid = mid, lx = lx),
                 aes(x = mid, xend = lx, y = 0.68, yend = 0.755),
                 linewidth = 0.25, colour = "grey45") +
    annotate("segment", x = 0, xend = 1, y = 0.5, yend = 0.5, linewidth = 0.9) +
    geom_text(data = data.frame(lx = lx, lab = sprintf("%s\n%.3f", names(f), f)),
              aes(lx, 0.80, label = lab), size = 3.4, lineheight = 0.95, vjust = 0) +
    geom_text(data = data.frame(y = c(0.5625, 0.4375), lab = c("geography", anc_label)),
              aes(-0.012, y, label = lab), hjust = 1, size = 3.4) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0.34, 0.94), clip = "off") +
    theme_void(base_size = 11) +
    theme(plot.margin = margin(10, 14, 10, 96), plot.title.position = "plot") +
    labs(title = title, subtitle = subtitle)
}
