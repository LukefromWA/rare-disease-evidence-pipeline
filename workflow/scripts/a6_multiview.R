#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A6: view agreement and coverage (diagnostic only; never fed back into ranking)
a <- parse_args(c("config", "sample", "shortlist", "hpo", "network", "outdir")); cfg <- read_config(a$config)
sl <- read_gz(a$shortlist); hpo <- read_gz(a$hpo); nw <- read_gz(a$network)
g <- unique(merge(unique(sl[, .(gene_symbol)]), merge(hpo[, .(gene_symbol, hpo = hpo_lin_bma_pct)], nw[, .(gene_symbol, string = string_agnostic_pct, multiplex = multiplex_agnostic_pct)], by = "gene_symbol"), by = "gene_symbol"))
pairs <- list(c("hpo", "string"), c("hpo", "multiplex"), c("string", "multiplex"))
ag <- rbindlist(lapply(pairs, function(p) data.table(view_a = p[1], view_b = p[2], spearman = suppressWarnings(cor(g[[p[1]]], g[[p[2]]], method = "spearman", use = "complete.obs")), n = sum(complete.cases(g[[p[1]]], g[[p[2]]])),
  top10_overlap = { k <- 10; a_ <- g[order(-get(p[1]))][seq_len(min(k, .N)), gene_symbol]; b_ <- g[order(-get(p[2]))][seq_len(min(k, .N)), gene_symbol]; length(intersect(a_, b_)) })))
cov <- g[, .(n_views = rowSums(!is.na(cbind(hpo, string, multiplex))))][, .N, by = n_views][order(n_views)]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(ag, file.path(od, "view_agreement.tsv")); write_table(cov, file.path(od, "coverage_counts.tsv"))
write_node_manifest(od, "a6", cfg, inputs = c(a$shortlist, a$hpo, a$network), outputs = file.path(od, c("view_agreement.tsv", "coverage_counts.tsv")))
