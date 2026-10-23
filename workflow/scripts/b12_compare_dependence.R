#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B12: validation-only comparison; the simplest variant unless the larger one wins by more than the declared margin. Locked test is not used.
a <- parse_args(c("config", "features", "fits", "outdir")); cfg <- read_config(a$config); feat <- read_gz(a$features); ft <- readRDS(a$fits); va <- feat[split == "validation"]
met <- rbindlist(lapply(names(ft$fits), function(nm) { cols <- ft$variants[[nm]]; p <- predict_softmax(as.matrix(va[, ..cols]), ft$fits[[nm]]$beta, va$case_id); m <- case_metrics(p, va$case_id, va$label)
  data.table(model = nm, n_params = length(cols), mrr = m$mrr, top1 = m$top1, top5 = m$top5, nll = m$nll, converged = ft$fits[[nm]]$converged) }))
sel <- if (met[model == "M2_interaction", mrr] > met[model == "M2_independent", mrr] + cfg$external$mrr_margin) "M2_interaction" else "M2_independent"
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(met, file.path(od, "dependence_report.tsv"))
saveRDS(list(model = sel, schema = ft$variants[[sel]], beta = ft$fits[[sel]]$beta, vcov = ft$fits[[sel]]$vcov), file.path(od, "selected_M2.rds"))
write_node_manifest(od, "b12", cfg, inputs = c(a$features, a$fits), outputs = file.path(od, c("dependence_report.tsv", "selected_M2.rds")), extra = list(selected = sel, margin = cfg$external$mrr_margin, used_locked_test = FALSE))
logmsg("B12: selected ", sel)
