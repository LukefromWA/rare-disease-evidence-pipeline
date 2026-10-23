#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B10 (runs BEFORE B9): whole duplicate groups assigned together; disjointness asserted inside define_casewise_splits
a <- parse_args(c("config", "candidates", "outdir")); cfg <- read_config(a$config); cand <- read_gz(a$candidates)
man <- define_casewise_splits(unique(cand[, .(case_id, duplicate_group)]), seed = cfg$external$split_seed)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(man, file.path(od, "split_manifest.tsv"))
write_node_manifest(od, "b10", cfg, inputs = a$candidates, outputs = file.path(od, "split_manifest.tsv"), extra = list(split_counts = as.list(table(man$split)), split_seed = cfg$external$split_seed))
logmsg("B10: ", paste(names(table(man$split)), table(man$split), collapse = " / "), " cases")
