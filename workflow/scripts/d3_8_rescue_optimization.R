#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("config", "sample", "outdir")); cfg <- read_config(a$config)
write_method_status(a$outdir, "d3_8_rescue_optimization", "not_run_cut_list", "cut-list item applied in plan v3.2; no rescue ranking available while D3.6 is disabled", cfg)
