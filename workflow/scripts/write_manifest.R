#!/usr/bin/env Rscript
# Generic manifest writer for bash/python nodes: write_manifest.R --config C --node a3 --outdir D --inputs "f1,f2" --outputs "g1,g2"
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("config", "node", "outdir", "inputs", "outputs")); cfg <- read_config(a$config)
write_node_manifest(a$outdir, a$node, cfg, inputs = strsplit(a$inputs, ",")[[1]], outputs = strsplit(a$outputs, ",")[[1]])
