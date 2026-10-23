#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# C2: private-only baselines (M0 = integrated score, M1 = informed posterior), normalized within the private set
a <- parse_args(c("config", "sample", "a10", "outdir")); cfg <- read_config(a$config); a10 <- read_gz(a$a10)
m0 <- a10$best_integrated_score; m1 <- if ("posterior_informed" %in% names(a10)) a10$posterior_informed else a10$best_integrated_score
out <- data.table(candidate_id = a10$best_candidate_id, gene = a10$gene_symbol, m0 = m0 / sum(m0, na.rm = TRUE), m1 = m1 / sum(m1, na.rm = TRUE)); out[, `:=`(rank_m0 = frank(-m0, ties.method = "min"), rank_m1 = frank(-m1, ties.method = "min"))]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(out, file.path(od, "m0_m1.tsv.gz")); write_node_manifest(od, "c2", cfg, inputs = a$a10, outputs = file.path(od, "m0_m1.tsv.gz"))
