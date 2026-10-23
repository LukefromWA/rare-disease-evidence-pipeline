#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D3.3: Jaccard overlap between the drug's target set and the genes of the mutated gene's pathways (Reactome release recorded)
a <- parse_args(c("config", "sample", "targets", "outdir")); cfg <- read_config(a$config); b <- require_bundle(cfg); tu <- read_gz(a$targets); pw <- read_gz(b$pathways); dr <- read_gz(b$drugs)
need_cols(pw, c("pathway_id", "gene_symbol"), "reactome_gene_sets"); genes_by_pw <- split(pw$gene_symbol, pw$pathway_id); pw_of_gene <- split(pw$pathway_id, pw$gene_symbol)
tgt <- dr[mapping_status == "mapped", .(targets = list(unique(target_gene))), by = drug_id]
out <- unique(tu[, .(scenario_id, drug_id, gene)])
out[, `:=`(pathway_jaccard = mapply(function(g, d) { P_ <- unique(unlist(genes_by_pw[pw_of_gene[[g]]], use.names = FALSE)); T_ <- tgt[drug_id == d, targets][[1]]; if (!length(P_)) NA_real_ else length(intersect(P_, T_)) / length(union(P_, T_)) }, gene, drug_id),
           n_shared_pathways = mapply(function(g, d) length(intersect(pw_of_gene[[g]], unique(unlist(pw_of_gene[tgt[drug_id == d, targets][[1]]])))), gene, drug_id))]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(out, file.path(od, "method3_pathway.tsv.gz"))
write_node_manifest(od, "d3_3", cfg, inputs = c(a$targets, b$pathways, b$drugs), outputs = file.path(od, "method3_pathway.tsv.gz"), file = "manifest_d3_3.json", extra = list(multiple_testing = "none (descriptive Jaccard; no p-values)"))
