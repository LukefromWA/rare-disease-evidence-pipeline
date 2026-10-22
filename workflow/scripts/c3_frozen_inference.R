#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# C3: frozen inference. STOP if B14 is absent or altered, if the schema mismatches, or if a refit is attempted (this script has no fitting call).
a <- parse_args(c("config", "sample", "features", "frozen", "freeze_manifest", "outdir"), c("refit")); cfg <- read_config(a$config)
if (!is.null(a$refit) && tolower(a$refit) %in% c("true", "1", "yes")) stop("C3 STOP: refit is forbidden on a frozen model", call. = FALSE)
fm <- jsonlite::read_json(need_files(a$freeze_manifest, "B14 freeze manifest")); assert_frozen(a$frozen, fm$sha256); fz <- readRDS(a$frozen); stopifnot(identical(fz$refit_allowed, FALSE))
f <- read_gz(a$features); sch <- fz$schema; need_cols(f, sch, "C1 features")
if (nrow(f) > 2 * fz$train_candidates_per_case) warning(sprintf("private set has %d candidates vs %d per training case: softmax probabilities are not comparable across set sizes; read ranks and q_j as relative", nrow(f), fz$train_candidates_per_case), call. = FALSE)
p <- predict_softmax(as.matrix(f[, ..sch]), fz$beta, f$case_id)
out <- data.table(candidate_id = f$candidate_id, gene = f$gene, p_frozen = p, rank_frozen = frank(-p, ties.method = "min"), in_conformal_set = (1 - p) <= fz$conformal_qhat)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(out, file.path(od, "frozen_scores.tsv.gz"))
write_node_manifest(od, "c3", cfg, inputs = c(a$features, a$frozen, a$freeze_manifest), outputs = file.path(od, "frozen_scores.tsv.gz"), extra = list(refit = FALSE, frozen_sha256 = fm$sha256, candidate_set_size = nrow(f), train_candidates_per_case = fz$train_candidates_per_case))
