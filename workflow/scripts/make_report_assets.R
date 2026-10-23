#!/usr/bin/env Rscript
# Writes reports/values.tex and reports/generated/*.tex from frozen outputs and manifests. Numbers never reach LaTeX any other way.
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("config", "sample", "outroot", "reports")); cfg <- read_config(a$config); o <- a$outroot; gen <- file.path(a$reports, "generated"); dir.create(gen, recursive = TRUE, showWarnings = FALSE)
qs <- fread(file.path(o, "a1", "qc_summary.tsv")); a10 <- read_gz(file.path(o, "a10", "ranked_candidates.tsv.gz")); c4 <- read_gz(file.path(o, "c4", "hypothesis_weights.tsv.gz"))
ex <- fread(file.path(o, "a2", "exclusion_summary.tsv")); d7 <- fread(file.path(o, "d7", "final_report.tsv")); m10 <- read_manifest(file.path(o, "a10", "manifest.json"))
v <- c(as.list(setNames(qs$value, paste0("a1_", qs$metric))), list(a10_n_genes = nrow(a10), a10_top_gene_informed = a10$gene_symbol[1], a10_top_gene_agnostic = a10[order(rank_agnostic)]$gene_symbol[1],
        c4_n_hypotheses = nrow(c4), c4_sum_q = round(sum(c4$q_j), 6), run_id = cfg$run$run_id, a10_manifest_sha8 = substr(m10$outputs[[1]]$sha256, 1, 8)))
write_values_tex(v, file.path(a$reports, "values.tex"))
write_latex_table(qs, file.path(gen, "a1_qc_summary.tex"), "A1 QC summary"); write_latex_table(ex, file.path(gen, "a2_filter_audit.tex"), "A2 exclusion summary")
write_latex_table(a10[, .(gene_priority_rank, gene_symbol, best_candidate_id, best_IMPACT, posterior_informed, posterior_agnostic, coverage_class, review_tier)], file.path(gen, "a10_ranked_candidates.tex"), "A10 ranked candidates (top rows)", max_rows = 15)
write_latex_table(c4[order(-q_j), .(gene, candidate_id, q_j, in_conformal_set)], file.path(gen, "c4_hypothesis_set.tex"), "C4 hypothesis weights (top rows)", max_rows = 15)
write_latex_table(d7, file.path(gen, "d7_track2_shortlist.tex"), "D7 research hypotheses (not treatment recommendations)", max_rows = 20)
logmsg("report assets written")
