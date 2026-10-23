#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("registry", "raw", "out")); reg <- yaml::read_yaml(a$registry)$track2_resources; s <- reg$sources$pediatric_safety; assert_frozen_resource(a$raw, "safety_raw")
write_table(build_pediatric_safety_map(read_gz(a$raw, colClasses = "character"), s$version %||% "unknown"), a$out); write_resource_manifest(a$out, "pediatric_safety_map", s$url, s$version, s$license)
