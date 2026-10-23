#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D3.1: direction gate on the normalized action vocabulary (uncertain never passes)
a <- parse_args(c("config", "sample", "targets", "outdir")); cfg <- read_config(a$config); tu <- read_gz(a$targets)
out <- tu[, .(scenario_id, drug_id, gene, direction, action_normalized)][, direction_gate := direction_gate(direction, action_normalized)]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(out, file.path(od, "method1_target_query.tsv.gz"))
write_node_manifest(od, "d3_1", cfg, inputs = a$targets, outputs = file.path(od, "method1_target_query.tsv.gz"), file = "manifest_d3_1.json", extra = list(gate_counts = as.list(table(out$direction_gate))))
