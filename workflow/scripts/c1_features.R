#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# C1: private features with the SAME compute_features() as B9, using the FROZEN state. Nothing is estimated from private data.
a <- parse_args(c("config", "sample", "a10", "frozen", "resources", "outdir")); cfg <- read_config(a$config)
a10 <- read_gz(a$a10); need_cols(a10, c("gene_symbol", "best_candidate_id", "best_integrated_score", "best_IMPACT"), "A10"); stopifnot(!anyDuplicated(a10$best_candidate_id))
fz <- readRDS(a$frozen); er <- readRDS(a$resources); terms <- read_proband_terms(rp(cfg, cfg$private$proband_hpo), rp(cfg, cfg$private$family_history_hpo))
cand <- data.table(case_id = "MVA", candidate_id = a10$best_candidate_id, gene = a10$gene_symbol)
f <- compute_features(cand, list(MVA = terms), er$res, fz$state); stopifnot(all(fz$schema %in% names(f)), nrow(f) == nrow(cand))     # full A10 set, never a top-k subset
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(f, file.path(od, "mva_feature_matrix.tsv.gz"))
write_node_manifest(od, "c1", cfg, inputs = c(a$a10, a$frozen, a$resources), outputs = file.path(od, "mva_feature_matrix.tsv.gz"), extra = list(frozen_state_used = TRUE, n_candidates = nrow(f), parity_function = "compute_features (shared with B9)"))
