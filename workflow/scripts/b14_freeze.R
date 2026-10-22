#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B14: freeze. Model, transform state, conformal threshold, seeds, checksums; read-only; reload test; no private strings.
a <- parse_args(c("config", "selected", "state", "conformal", "outdir")); cfg <- read_config(a$config)
sel <- readRDS(a$selected); st <- readRDS(a$state); cf <- jsonlite::read_json(a$conformal)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); rds <- file.path(od, "frozen_M2.rds")
if (file.exists(rds)) Sys.chmod(rds, "644")
fr <- list(model = sel$model, schema = sel$schema, beta = sel$beta, vcov = sel$vcov, state = st, conformal_qhat = cf$qhat, conformal_alpha = cf$alpha,
           train_candidates_per_case = st$train_candidates_per_case, seed = cfg$external$split_seed, refit_allowed = FALSE)
saveRDS(fr, rds); chk <- readRDS(rds); stopifnot(identical(chk$schema, fr$schema), isTRUE(all.equal(chk$beta, fr$beta)))
mf <- list(run_id = cfg$run$run_id, model = fr$model, sha256 = sha256_file(rds), conformal_qhat = cf$qhat, locked_test_coverage = cf$locked_test_coverage, locked_test_top1 = cf$locked_test_top1,
           train_candidates_per_case = fr$train_candidates_per_case, r_version = R.version.string, time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
fm <- file.path(od, "freeze_manifest.json"); if (file.exists(fm)) Sys.chmod(fm, "644")
assert_no_private(mf); jsonlite::write_json(mf, file.path(od, "freeze_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
write_bundle_checksums(c(rds, file.path(od, "freeze_manifest.json")), file.path(od, "checksums.sha256"), root = od)
Sys.chmod(c(rds, file.path(od, "freeze_manifest.json")), "444"); write_node_manifest(od, "b14", cfg, inputs = c(a$selected, a$state, a$conformal), outputs = rds, extra = list(frozen = TRUE, refit_allowed = FALSE))
logmsg("B14: frozen -> ", rds)
