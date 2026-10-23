#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B12.5: conformal threshold on validation positives; coverage on the locked test, evaluated once
a <- parse_args(c("config", "features", "selected", "outdir")); cfg <- read_config(a$config); feat <- read_gz(a$features); sel <- readRDS(a$selected); cols <- sel$schema
va <- feat[split == "validation"]; te <- feat[split == "locked_test"]
pv <- predict_softmax(as.matrix(va[, ..cols]), sel$beta, va$case_id); qhat <- conformal_qhat(1 - pv[va$label == 1L], cfg$external$conformal_alpha)
pt <- predict_softmax(as.matrix(te[, ..cols]), sel$beta, te$case_id); cover <- mean(1 - pt[te$label == 1L] <= qhat); tm <- case_metrics(pt, te$case_id, te$label)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(list(qhat = qhat, alpha = cfg$external$conformal_alpha, locked_test_coverage = cover, locked_test_top1 = tm$top1, locked_test_mrr = tm$mrr, n_validation_cases = uniqueN(va$case_id), n_test_cases = uniqueN(te$case_id)), file.path(od, "conformal.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
write_node_manifest(od, "b12_5", cfg, inputs = c(a$features, a$selected), outputs = file.path(od, "conformal.json"), extra = list(locked_test_used_once = TRUE))
logmsg("B12.5: qhat=", signif(qhat, 3), " locked-test coverage=", round(cover, 3))
