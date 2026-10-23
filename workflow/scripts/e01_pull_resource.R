#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# E01: acquisition run only (allow_live_downloads must be true in the registry). Main scoring uses frozen local snapshots.
a <- parse_args(c("registry", "root", "source")); reg <- yaml::read_yaml(a$registry)$track2_resources; s <- reg$sources[[a$source]]; if (is.null(s)) stop("unknown source: ", a$source)
p <- download_resource(s$url, file.path(a$root, s$raw_path), isTRUE(reg$allow_live_downloads), expected_sha256 = s$expected_sha256)
write_resource_manifest(p, a$source, s$url, s$version, s$license); logmsg("pulled ", a$source)
