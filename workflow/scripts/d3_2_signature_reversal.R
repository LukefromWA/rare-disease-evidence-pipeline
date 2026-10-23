#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D3.2: needs a pinned local LINCS/CLUE snapshot, which is not configured -> explicit status
a <- parse_args(c("config", "sample", "outdir")); cfg <- read_config(a$config)
write_method_status(a$outdir, "d3_2_signature_reversal", "unavailable_data", "no pinned LINCS/CLUE snapshot configured", cfg)
