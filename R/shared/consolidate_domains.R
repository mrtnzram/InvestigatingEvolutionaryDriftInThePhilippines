# =============================================================================
# Shared — run a method across every domain and consolidate one results table
#
# Each analysis method ([4] PVR, [5] MMRR, [9] Procrustes) is implemented once
# per domain, and each copy used to write its own one-or-two-row summary CSV —
# 15 near-duplicate files holding a single domain's row of the same statistic.
# The `[N]_ALL_*.R` runners in this folder replace that: they source each
# domain's scripts, harvest the results object left in the environment, stack
# the rows under a `domain` key and write one consolidated table.
#
# The domain scripts no longer persist their own row, so a method's results
# table only exists after its runner has been run.
#
# Consumers: R/shared/[4]_ALL_PVR.R, [5]_ALL_MMRR.R, [9]_ALL_PROCRUSTES.R
# =============================================================================

library(dplyr)
library(here)

# Canonical domain order — every consolidated table is sorted by this, so rows
# line up across methods regardless of the order the runners ran them in.
DOMAIN_ORDER <- c("phoneme", "grammar", "cognate", "genetic")

# Marker the two database scripts use to separate the language-set definition
# (PART A) from the post-tree feature build (PART B). PART B asserts that the
# tree script has already run, so a runner that needs only `Ph_Languages` must
# stop at this line.
PART_A_MARKER <- ">>> END OF PART A <<<"

#' Mark a script so run_domain() sources only its PART A.
#'
#' [0]_CREANZA_RUHLENdatabase.R and [0]_GRAMBANKdatabase.R are split around
#' PART_A_MARKER: PART A defines `Ph_Languages` (and, for grammar,
#' `GRAMBANKdf_PH_maximized`) which [0]_Phylogenetic_Tree.R consumes, and PART B
#' then consumes the tree's `Ph_Languages_pruned`. Only PART A is a prerequisite
#' of the [4] scripts, and sourcing the whole file would fail on PART B's guard.
part_a <- function(path) structure(path, class = "part_a_script")

#' Source a domain's scripts into one isolated environment.
#'
#' A fresh environment per domain is what makes this safe: all four [9] scripts
#' name their table `results_df`, and every domain's [0]_Phylogenetic_Tree.R
#' defines `tree_pruned`. The parent is globalenv() so `library()` attachments
#' and `here()` resolve normally, and so the `stopifnot(exists("tree_pruned"))`
#' guards in the [4] scripts see objects the tree script left in the same env.
#'
#' @param scripts Path, or list/vector of paths, sourced in order. Prerequisites
#'                come first. Wrap an entry in part_a() to stop at PART_A_MARKER.
#' @return The environment the scripts were evaluated in.
run_domain <- function(scripts) {
  env <- new.env(parent = globalenv())
  for (s in as.list(scripts)) {
    if (!file.exists(s)) stop("Script not found: ", s, call. = FALSE)
    if (inherits(s, "part_a_script")) {
      message("  sourcing ", basename(s), " (PART A only)")
      lines <- readLines(s, warn = FALSE)
      cut <- grep(PART_A_MARKER, lines, fixed = TRUE)
      if (length(cut) != 1L) {
        stop("Expected exactly one '", PART_A_MARKER, "' line in ", basename(s),
             "; found ", length(cut), ".", call. = FALSE)
      }
      eval(parse(text = paste(lines[seq_len(cut - 1L)], collapse = "\n")), envir = env)
    } else {
      message("  sourcing ", basename(s))
      source(s, local = env, echo = FALSE)
    }
  }
  env
}

#' Pull a named results object out of a run_domain() environment.
#'
#' @param env    Environment from run_domain().
#' @param object Name of the results data frame the script builds.
#' @param domain Lowercase domain label, prepended as the `domain` column.
#' @return The results data frame with `domain` as its first column.
harvest <- function(env, object, domain) {
  if (!exists(object, envir = env, inherits = FALSE)) {
    stop("`", object, "` was not defined for domain '", domain,
         "'. Did the script's results block change name?", call. = FALSE)
  }
  df <- get(object, envir = env)
  stopifnot("harvested object is not a data frame" = is.data.frame(df))
  # Domains that already carry their own key (Procrustes' `domain`, MMRR's
  # uppercase `dataset`) get it replaced by the canonical lowercase label.
  df %>%
    select(-any_of(c("domain", "dataset"))) %>%
    mutate(domain = domain, .before = 1)
}

#' Stack per-domain rows and write the consolidated table.
#'
#' bind_rows() takes the union of columns, so a domain that lacks one gets NA —
#' that is how `zeros` (cognate/genetic only), `beta_phylo`/`p_phylo` (absent for
#' genetic, which has no tree) and `n_individual` (genetic only) are handled.
#'
#' @param rows List of data frames from harvest().
#' @param path Output path, from here().
#' @return The consolidated data frame, invisibly.
write_consolidated <- function(rows, path) {
  out <- bind_rows(rows) %>%
    mutate(domain = factor(domain, levels = DOMAIN_ORDER)) %>%
    arrange(domain) %>%
    mutate(domain = as.character(domain))

  stopifnot("unrecognised domain label" = !any(is.na(out$domain)))

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(out, file = path, row.names = FALSE)
  message("\nwrote ", nrow(out), " rows -> ",
          sub(here(), "", path, fixed = TRUE))
  print(out)
  invisible(out)
}
