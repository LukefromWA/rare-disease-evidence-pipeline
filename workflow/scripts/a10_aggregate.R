#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A10: gene-level ranked table (best variant per gene), coverage and stability classes, review tier. Same column names as the recorded A10.
a <- parse_args(c("config", "sample", "integrated", "mosaicism", "stability", "outdir")); cfg <- read_config(a$config)
sl <- read_gz(a$integrated); mo <- read_gz(a$mosaicism); fr <- read_gz(a$stability)
sl <- merge(sl, mo[, .(candidate_id, VAF_post_mean, VAF_95CrI_lower, VAF_95CrI_upper, log_bf_mosaic_vs_het, mosaic_review)], by = "candidate_id", all.x = TRUE, sort = FALSE)
for (j in c("GQ", "MAX_AF", "VAF", "DP", "AD_alt")) set(sl, j = j, value = as.numeric(sl[[j]]))
a10 <- aggregate_ranked_candidates(sl, fr); stopifnot(!anyDuplicated(a10$gene_symbol), all(c("posterior_informed", "posterior_agnostic") %in% names(a10)))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(a10, file.path(od, "ranked_candidates.tsv.gz"))
write_node_manifest(od, "a10", cfg, inputs = c(a$integrated, a$mosaicism, a$stability), outputs = file.path(od, "ranked_candidates.tsv.gz"),
  extra = list(n_genes = nrow(a10), top_genes_informed = as.list(head(a10$gene_symbol, 5)), both_posteriors_present = TRUE))
logmsg("A10: ", nrow(a10), " genes; top: ", paste(head(a10$gene_symbol, 3), collapse = ", "))
