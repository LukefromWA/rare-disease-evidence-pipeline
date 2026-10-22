#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A4: gene-level HPO evidence for the shortlist genes (proband-only terms)
a <- parse_args(c("config", "sample", "shortlist", "outdir")); cfg <- read_config(a$config); sl <- read_gz(a$shortlist)
genes <- unique(sl$gene_symbol[!is.na(sl$gene_symbol)])
terms <- read_proband_terms(rp(cfg, cfg$private$proband_hpo), rp(cfg, cfg$private$family_history_hpo))
ont <- build_ontology(parse_obo(need_files(rp(cfg, cfg$resources$hpo_obo), "HPO obo")))
g2p <- data.table::fread(need_files(rp(cfg, cfg$resources$gene_to_hpo), "gene_to_hpo"), sep = "\t", quote = "")
hpo <- compute_hpo_evidence(genes, terms, ont, gene_term_sets(g2p, ont))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(hpo, file.path(od, "hpo_gene_scores.tsv"))
write_node_manifest(od, "a4", cfg, inputs = c(a$shortlist, rp(cfg, c(cfg$private$proband_hpo, cfg$resources$hpo_obo, cfg$resources$gene_to_hpo))), outputs = file.path(od, "hpo_gene_scores.tsv"),
  extra = list(n_patient_terms_provided = length(terms), n_patient_terms_active = length(clean_terms(terms, ont)), patient_terms = as.list(terms), n_genes = length(genes), n_scored = sum(is.finite(hpo$hpo_lin_bma)),
               family_history_terms_used = FALSE))
