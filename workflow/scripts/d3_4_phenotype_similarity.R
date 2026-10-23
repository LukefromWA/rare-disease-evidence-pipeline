#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("config", "sample", "outdir")); cfg <- read_config(a$config)
write_method_status(a$outdir, "d3_4_phenotype_similarity", "unavailable_data", "no drug-indication HPO resource configured", cfg)
