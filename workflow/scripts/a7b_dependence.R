#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A7b: Gaussian-copula dependence between HPO and network channels
a <- parse_args(c("config", "sample", "integrated", "outdir")); cfg <- read_config(a$config); sl <- read_gz(a$integrated)
dep <- rbindlist(lapply(c("string", "multiplex"), function(v) { d <- estimate_dependence(sl$hpo_lin_bma_pct, sl[[paste0(v, "_agnostic_pct")]]); data.table(network = v, rho = d$rho, spearman = d$spearman, n = d$n) }))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(dep, file.path(od, "dependence.tsv"))
write_node_manifest(od, "a7b", cfg, inputs = a$integrated, outputs = file.path(od, "dependence.tsv"), file = "manifest_a7b.json")
