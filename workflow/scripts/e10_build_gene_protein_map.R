#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("registry", "raw", "out_map", "out_lookup")); reg <- yaml::read_yaml(a$registry)$track2_resources; s <- reg$sources$uniprot
assert_frozen_resource(a$raw, "uniprot_raw"); r <- build_gene_protein_map(read_gz(a$raw, colClasses = "character"))
write_table(r$gene_protein_map, a$out_map); write_table(r$uniprot_lookup, a$out_lookup)
for (p in c(a$out_map, a$out_lookup)) write_resource_manifest(p, basename(p), s$url, s$version, s$license)
