#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "stan_runner.R"))
# B11: M2 variants on the TRAIN split (analytic-gradient conditional softmax). Optional --backend stan runs grouped OOF with the Stan models.
a <- parse_args(c("config", "features", "outdir"), c("backend")); cfg <- read_config(a$config); feat <- read_gz(a$features)
tr <- feat[split == "train"]; variants <- list(M2_independent = c("G", "H", "N"), M2_interaction = c("G", "H", "N", "GxH"))
fits <- lapply(variants, function(cols) fit_conditional_softmax(tr[, ..cols], tr$case_id, tr$label, lambda = 1e-3))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); saveRDS(list(fits = fits, variants = variants), file.path(od, "m2_fits.rds"))
write_table(rbindlist(lapply(names(fits), function(n) data.table(model = n, term = names(fits[[n]]$beta), beta = as.numeric(fits[[n]]$beta), converged = fits[[n]]$converged))), file.path(od, "m2_coefficients.tsv"))
outs <- file.path(od, c("m2_fits.rds", "m2_coefficients.tsv"))
if (identical(a$backend, "stan")) {                                  # NOT exercised by the toy test; see stan_runner.R and the Stan README
  Sys.setenv(MVA_STAN_DIR = normalizePath(file.path(dirname(normalizePath(.sf)), "..", "stan")))
  dev <- copy(tr); dev[, `:=`(group = "all", hpo_score = H, network_score = N)]
  gs <- unique(dev[, .(duplicate_group = case_id)]); f <- rep_len(seq_len(3), nrow(gs)); dev[, inner_fold := f[match(case_id, gs$duplicate_group)]]
  cm <- compile_models(c("m2a", "m2b", "m3a", "m3b")); oof <- run_oof(dev, c("G", "H", "N"), "all", cm, profile = "quick", out_dir = file.path(od, "stan_oof"), seed = cfg$run$seed)
  outs <- c(outs, file.path(od, "stan_oof", c("oof_predictions.tsv.gz", "oof_diagnostics.tsv", "oof_param_summaries.tsv.gz")))
}
write_node_manifest(od, "b11", cfg, inputs = a$features, outputs = outs, extra = list(backend = a$backend %||% "clogit_optim", fit_split = "train", converged = as.list(vapply(fits, function(f) f$converged, TRUE))))
