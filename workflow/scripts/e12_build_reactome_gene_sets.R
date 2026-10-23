#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("registry", "raw", "lookup", "out")); reg <- yaml::read_yaml(a$registry)$track2_resources; s <- reg$sources$reactome; assert_frozen_resource(a$raw, "reactome_raw")
raw <- data.table::fread(a$raw, header = FALSE, sep = "\t", quote = "", colClasses = "character", col.names = c("uniprot_id", "pathway_id", "url", "pathway_name", "evidence", "species"))
write_table(build_reactome_gene_sets(raw, read_gz(a$lookup), s$version %||% "unknown"), a$out); write_resource_manifest(a$out, "reactome_gene_sets", s$url, s$version, s$license)
