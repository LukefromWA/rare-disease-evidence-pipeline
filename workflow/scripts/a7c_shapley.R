#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A7c: exact Shapley attribution of the top informed candidate's channels
a <- parse_args(c("config", "sample", "integrated", "outdir")); cfg <- read_config(a$config); sl <- read_gz(a$integrated)
W <- unlist(cfg$integration$channel_weights); top <- which.min(sl$rank_informed)
sc <- c(genomic = sl$genomic[top], hpo = sl$hpo_lin_bma_pct[top], string = sl$string_informed_pct[top], multiplex = sl$multiplex_informed_pct[top]); sc[!is.finite(sc)] <- NA
sh <- shapley_attribution(sc, W[names(sc)]); stopifnot(abs(sum(sh) - attr(sh, "full_minus_empty")) < 1e-9)       # efficiency axiom
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(data.table(candidate_id = sl$candidate_id[top], channel = names(sh), shapley = as.numeric(sh)), file.path(od, "shapley_top.tsv"))
write_node_manifest(od, "a7c", cfg, inputs = a$integrated, outputs = file.path(od, "shapley_top.tsv"), file = "manifest_a7c.json")
