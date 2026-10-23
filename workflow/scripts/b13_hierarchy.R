#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B13: hierarchical (M3) blocks. In the compressed schedule M3 is NOT run; this node records that explicitly (never silent).
a <- parse_args(c("config", "outdir")); cfg <- read_config(a$config); od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
h <- data.table(block = c("disease_group", "phenotype_group", "mechanism_group"), status = "not_run", reason = "cut-list item applied in plan v3.2; use stan/m3a.stan, m3b.stan with grouped OOF if time allows")
write_table(h, file.path(od, "hierarchy_report.tsv")); write_node_manifest(od, "b13", cfg, outputs = file.path(od, "hierarchy_report.tsv"), extra = list(m3_run = FALSE))
