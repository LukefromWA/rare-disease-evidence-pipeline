#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A8: informative-depth beta-binomial VAF posterior and germline-het vs mosaic log Bayes factor
a <- parse_args(c("config", "sample", "shortlist", "outdir")); cfg <- read_config(a$config); sl <- read_gz(a$shortlist)
for (j in c("DP", "GQ", "AD_ref", "AD_alt")) set(sl, j = j, value = as.numeric(sl[[j]]))
qc <- list(dp_gt = cfg$qc$min_dp, min_gq = cfg$qc$min_gq, mosaic_min_alt_reads = cfg$qc$mosaic_min_alt_reads, mosaic_vaf_range = unlist(cfg$qc$mosaic_vaf_range))
t9 <- run_tier9(sl, qc)[, .(candidate_id, gene_symbol, VAF_observed, VAF_post_mean, VAF_95CrI_lower, VAF_95CrI_upper, log_bf_mosaic_vs_het, mosaic_review)]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(t9, file.path(od, "mosaicism.tsv"))
write_node_manifest(od, "a8", cfg, inputs = a$shortlist, outputs = file.path(od, "mosaicism.tsv"), extra = list(n_mosaic_review = sum(t9$mosaic_review, na.rm = TRUE)))
