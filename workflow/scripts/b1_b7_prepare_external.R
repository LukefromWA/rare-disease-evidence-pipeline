#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B1-B7 (consume): external cases + public resource index. HARD BOUNDARY: reads no private path; asserts none appear in outputs.
a <- parse_args(c("config", "outdir")); cfg <- read_config(a$config)
sg <- load_string_graph(need_files(rp(cfg, cfg$resources$string_links), "STRING links"), need_files(rp(cfg, cfg$resources$string_info), "STRING info"), cfg$ppr$string_min_score)
er <- load_external_resources(rp(cfg, cfg$resources$gene_to_hpo), sg$edges)
cases <- read_external_cases(rp(cfg, cfg$external$dedup_cases), cfg$external$min_observed_hpo); cases <- cases[true_gene %in% er$universe]
logmsg("B1-B7: cases with causal gene in the HPO+STRING universe: ", nrow(cases), " | universe genes: ", length(er$universe))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(cases, file.path(od, "cases.tsv.gz")); saveRDS(er, file.path(od, "resources.rds"), compress = "xz")
write_node_manifest(od, "b1_b7", cfg, inputs = rp(cfg, c(cfg$external$dedup_cases, cfg$resources$gene_to_hpo, cfg$resources$string_links, cfg$resources$string_info)), outputs = file.path(od, c("cases.tsv.gz", "resources.rds")),
  extra = list(n_cases = nrow(cases), universe_size = length(er$universe), phenopacket_version = cfg$external$phenopacket_store_version %||% NA))
