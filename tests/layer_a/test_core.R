#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
root <- normalizePath(file.path(
  dirname(sub("^--file=", "", arg)), "..", ".."
), mustWork = TRUE)

e <- new.env(parent = globalenv())
e$logmsg <- function(...) invisible(NULL)
source(file.path(root, "workflow/layer_a/R/01_helpers_qc.R"), local = e)
source(file.path(root, "workflow/layer_a/R/04_integration.R"), local = e)

toy <- data.table(
  candidate_id = c("TOY_A", "TOY_B"), FILTER = "PASS",
  DP = c(20, 10), GQ = 35, AD_ref = 12, AD_alt = 8,
  MAX_AF = .001, IMPACT = "HIGH"
)
qc <- list(dp_gt = 10, min_gq = 20, min_alt_reads = 4, max_af = .01)
x <- e$apply_core_qc(toy, qc)
stopifnot(
  identical(x$pass_core_qc, c(TRUE, FALSE)),
  all(abs(x$VAF - .4) < 1e-12),
  !"pass_core_qc" %in% names(toy)
)
cat("PASS: QC boundaries, informative VAF and copy semantics\n")

scores <- rbind(
  c(genomic = .8, hpo = .2),
  c(.8, NA_real_),
  c(NA_real_, NA_real_)
)
x <- e$weighted_available_mean(scores, c(hpo = .75, genomic = .25))
stopifnot(
  abs(x[1] - .35) < 1e-12,
  abs(x[2] - .8) < 1e-12,
  is.na(x[3]),
  all(attr(x, "n_views") == c(2, 1, 0))
)
cat("PASS: name-aligned weights and missing evidence\n")

toy <- data.table(
  candidate_id = c("TOY_A", "TOY_B"),
  gene = c("TOY_G1", "TOY_G2"),
  genomic = c(.8, .3), hpo = c(.6, .2)
)
x <- e$integrate_tier1(toy, c(genomic = .5, hpo = .5), character())
stopifnot(
  isTRUE(all.equal(x$posterior_informed, x$posterior_agnostic)),
  abs(sum(x$posterior_agnostic) - 1) < 1e-12
)
cat("PASS: empty-panel informed/agnostic agreement\n")

x <- e$shapley_attribution(
  c(genomic = .8, hpo = .4, network = NA_real_),
  c(genomic = .5, hpo = .3, network = .2)
)
stopifnot(abs(sum(x) - attr(x, "full_minus_empty")) < 1e-10)
cat("PASS: Shapley efficiency with missing evidence\n")
