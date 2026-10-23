#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("registry", "raw", "lookup", "out_edges", "out_conflicts", "out_audit")); reg <- yaml::read_yaml(a$registry)$track2_resources; s <- reg$sources$signor
assert_frozen_resource(a$raw, "signor_raw"); r <- build_signed_edges(read_gz(a$raw, colClasses = "character"), read_gz(a$lookup), s$version %||% "unknown"); sp <- split_edge_conflicts(r$edges)
stopifnot(!anyDuplicated(sp$usable_edges[, .(source_gene, target_gene)]))
write_table(sp$usable_edges, a$out_edges); write_table(sp$conflict_audit, a$out_conflicts); write_table(r$drop_audit, a$out_audit)
for (p in c(a$out_edges, a$out_conflicts)) write_resource_manifest(p, basename(p), s$url, s$version, s$license)
