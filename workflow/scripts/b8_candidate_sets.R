#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B8: candidate sets per case; background = HPO+STRING-mapped gene universe minus the positive (policy recorded)
a <- parse_args(c("config", "cases", "resources", "outdir")); cfg <- read_config(a$config); cases <- read_gz(a$cases); er <- readRDS(a$resources)
cand <- build_candidate_sets(cases, er$universe, cfg$external$negatives_per_case, cfg$external$split_seed)
stopifnot(all(cand[, sum(label), by = case_id]$V1 == 1L), !anyDuplicated(cand$candidate_id)); assert_no_private(cand)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(cand, file.path(od, "candidate_sets.tsv.gz"))
write_node_manifest(od, "b8", cfg, inputs = c(a$cases, a$resources), outputs = file.path(od, "candidate_sets.tsv.gz"),
  extra = list(negatives_per_case = cfg$external$negatives_per_case, background_policy = "HPO+STRING-mapped gene universe minus the case's positive gene", n_rows = nrow(cand)))
