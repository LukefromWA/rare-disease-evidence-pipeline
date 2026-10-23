#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# B9: features. Transform state is fit on the TRAIN split only and applied unchanged to validation and locked test.
a <- parse_args(c("config", "cases", "candidates", "splits", "resources", "outdir")); cfg <- read_config(a$config)
cases <- read_gz(a$cases); cand <- read_gz(a$candidates); man <- read_gz(a$splits); er <- readRDS(a$resources)
cand <- merge(cand, man[, .(case_id, split)], by = "case_id"); case_terms <- stats::setNames(strsplit(cases$observed_hpo_ids_normalized, "|", fixed = TRUE), cases$case_id)
state <- fit_feature_state(cand[split == "train"], case_terms, er$res)
feat <- compute_features(cand, case_terms, er$res, state); feat <- merge(feat, cand[, .(candidate_id, split, label)], by = "candidate_id")
tr_cases <- unique(cand[split == "train", case_id]); stopifnot(!any(unique(cand[split != "train", case_id]) %in% tr_cases))     # no validation/test case informed the state
state$train_candidates_per_case <- cfg$external$negatives_per_case + 1L; state$train_case_ids_n <- length(tr_cases)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(feat, file.path(od, "features.tsv.gz")); saveRDS(state, file.path(od, "feature_state.rds"))
write_node_manifest(od, "b9", cfg, inputs = c(a$cases, a$candidates, a$splits, a$resources), outputs = file.path(od, c("features.tsv.gz", "feature_state.rds")),
  extra = list(state_fit_on = "train split only", schema = as.list(state$schema), gene_prior_note = "G = log1p(leave-one-case-out causal-gene count in TRAIN); rename to gene_prior in reports"))
