#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# E20: graph + audit + metadata from ONE edge snapshot
a <- parse_args(c("edges", "drugs", "out_graph", "out_audit", "out_meta")); e <- read_gz(a$edges); dr <- read_gz(a$drugs)
e[, `:=`(sign = as.integer(sign), edge_weight = as.numeric(edge_weight), edge_probability = as.numeric(edge_probability))]
g <- build_signed_graph(e); saveRDS(g, a$out_graph); au <- audit_signed_graph(g, e, dr); write_table(au, a$out_audit)
jsonlite::write_json(list(n_edges_table = nrow(e), n_edges_graph = igraph::ecount(g), edges_sha256 = sha256_file(a$edges), audit = as.list(au)), a$out_meta, auto_unbox = TRUE, pretty = TRUE)
stopifnot(igraph::ecount(g) == nrow(e)); for (p in c(a$out_graph, a$out_audit, a$out_meta)) write_resource_manifest(p, basename(p))
