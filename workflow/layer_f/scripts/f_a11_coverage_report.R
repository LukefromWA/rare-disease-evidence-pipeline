#!/usr/bin/env Rscript

# Creates a report-ready ranked candidate table with raw evidence, Monte Carlo
# rank-stability metrics, coverage classes, and a top-25 TSV for review.
# This script uses only existing private outputs; it does not download data.

options(stringsAsFactors = FALSE)

.source_file <- normalizePath(sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
.runtime_file <- normalizePath(file.path(dirname(.source_file),"..","..","lib","a11_runtime.R"))
source(.runtime_file)
.cli <- a11_cli(c("input","consensus","copula","a11-manifest","outdir"),list(`focus-gene`=""))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Required preinstalled package: data.table")
}
suppressPackageStartupMessages(library(data.table))

INPUT_FILE <- .cli[["input"]]
CONSENSUS_FILE <- .cli[["consensus"]]
COPULA_FILE <- .cli[["copula"]]
OUT_DIR <- .cli[["outdir"]]
FOCUS_GENE <- .cli[["focus-gene"]]

for (f in c(INPUT_FILE, CONSENSUS_FILE, COPULA_FILE)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}
.producer <- a11_verify(.cli[["a11-manifest"]],INPUT_FILE,CONSENSUS_FILE,COPULA_FILE)
.ctx <- a11_begin(OUT_DIR,c(evidence=INPUT_FILE,consensus=CONSENSUS_FILE,copula=COPULA_FILE,producer=.cli[["a11-manifest"]]),.source_file,.runtime_file)

input <- fread(INPUT_FILE)
consensus <- fread(CONSENSUS_FILE)
dependence <- fread(COPULA_FILE)

required_input <- c(
  "candidate_id", "gene_symbol", "genomic", "hpo", "string", "multiplex",
  "has_genomic", "has_hpo", "has_string", "has_multiplex"
)
required_consensus <- c(
  "candidate_id", "gene_symbol", "posterior_rank", "posterior_score_mean",
  "posterior_rank_mean", "p_rank1", "p_top3", "p_top5", "p_top10",
  "score_lower_envelope", "score_upper_envelope",
  "rank_lower_envelope", "rank_upper_envelope"
)

missing_input <- setdiff(required_input, names(input))
missing_consensus <- setdiff(required_consensus, names(consensus))
if (length(missing_input)) stop("Input missing columns: ", paste(missing_input, collapse = ", "))
if (length(missing_consensus)) stop("Consensus missing columns: ", paste(missing_consensus, collapse = ", "))

input <- unique(input[, ..required_input])
consensus <- unique(consensus[, ..required_consensus])

if (anyDuplicated(input$candidate_id)) stop("Duplicate candidate_id values in evidence input.")
if (anyDuplicated(consensus$candidate_id)) stop("Duplicate candidate_id values in consensus input.")
.mi <- match(consensus$candidate_id,input$candidate_id)
if (anyNA(.mi) || !identical(as.character(consensus$gene_symbol),as.character(input$gene_symbol[.mi]))) stop("Consensus/evidence gene-key mismatch")

out <- merge(
  consensus,
  input[, .(candidate_id, genomic, hpo, string, multiplex, has_genomic, has_hpo, has_string, has_multiplex)],
  by = "candidate_id",
  all.x = TRUE,
  suffixes = c("", "_input")
)

if ("gene_symbol_input" %in% names(out)) {
  out[is.na(gene_symbol) | gene_symbol == "", gene_symbol := gene_symbol_input]
  out[, gene_symbol_input := NULL]
}

out[, evidence_view_count := as.integer(has_genomic) + as.integer(has_hpo) + as.integer(has_string | has_multiplex)]
out[, network_view_count := as.integer(has_string) + as.integer(has_multiplex)]
out[, network_agreement := abs(string - multiplex)]

out[, coverage_class := fcase(
  has_genomic & has_hpo & has_string & has_multiplex &
    genomic >= 0.70 & hpo >= 0.70 & pmin(string, multiplex) >= 0.70,
  "three_view_concordant",

  has_genomic & has_hpo & has_string & has_multiplex,
  "three_view_mixed",

  has_genomic & has_hpo & !has_string & !has_multiplex,
  "two_view_genomic_hpo",

  has_genomic & has_hpo & (has_string | has_multiplex),
  "three_view_partial_network",

  has_genomic & !has_hpo & (has_string | has_multiplex),
  "two_view_genomic_network",

  !has_genomic & has_hpo & (has_string | has_multiplex),
  "two_view_hpo_network",

  evidence_view_count == 1L,
  "single_view_only",

  default = "no_usable_evidence"
)]

out[, review_tier := fcase(
  posterior_rank == 1L,
  "rank_1",

  posterior_rank == 2L,
  "rank_2",

  posterior_rank == 3L,
  "rank_3",

  posterior_rank <= 10L,
  "ranks_4_to_10",

  posterior_rank <= 25L,
  "ranks_11_to_25",

  default = "ranks_26_plus"
)]

out[, rank_stability_class := fcase(
  p_rank1 >= 0.90,
  "rank1_robust",

  p_top3 >= 0.90,
  "top3_robust",

  p_top5 >= 0.90,
  "top5_robust",

  p_top10 >= 0.90,
  "top10_robust",

  default = "scenario_sensitive"
)]

setorder(out, posterior_rank)

full_columns <- c(
  "posterior_rank", "candidate_id", "gene_symbol", "review_tier",
  "coverage_class", "rank_stability_class", "evidence_view_count",
  "network_view_count", "genomic", "hpo", "string", "multiplex",
  "network_agreement", "posterior_score_mean", "posterior_rank_mean",
  "p_rank1", "p_top3", "p_top5", "p_top10",
  "score_lower_envelope", "score_upper_envelope",
  "rank_lower_envelope", "rank_upper_envelope",
  "has_genomic", "has_hpo", "has_string", "has_multiplex"
)

full <- out[, ..full_columns]
top25 <- full[posterior_rank <= 25L]
top10 <- full[posterior_rank <= 10L]
focus <- if (nzchar(FOCUS_GENE)) full[toupper(gene_symbol)==toupper(FOCUS_GENE)] else full[0L]

summary_table <- data.table(
  metric = c(
    "n_ranked_candidates",
    "n_three_view_concordant",
    "n_three_view_mixed",
    "n_two_view_genomic_hpo",
    "n_rank1_robust",
    "n_top3_robust",
    "focus_gene_p_rank1",
    "focus_gene_p_top3",
    "focus_gene_p_top5",
    "focus_gene_p_top10",
    "focus_gene_mean_rank",
    "focus_gene_rank_envelope"
  ),
  value = c(
    nrow(full),
    sum(full$coverage_class == "three_view_concordant"),
    sum(full$coverage_class == "three_view_mixed"),
    sum(full$coverage_class == "two_view_genomic_hpo"),
    sum(full$rank_stability_class == "rank1_robust"),
    sum(full$rank_stability_class == "top3_robust"),
    if (nrow(focus)) sprintf("%.6f", focus$p_rank1[1]) else NA_character_,
    if (nrow(focus)) sprintf("%.6f", focus$p_top3[1]) else NA_character_,
    if (nrow(focus)) sprintf("%.6f", focus$p_top5[1]) else NA_character_,
    if (nrow(focus)) sprintf("%.6f", focus$p_top10[1]) else NA_character_,
    if (nrow(focus)) sprintf("%.6f", focus$posterior_rank_mean[1]) else NA_character_,
    if (nrow(focus)) paste0(focus$rank_lower_envelope[1], "-", focus$rank_upper_envelope[1]) else NA_character_
  )
)

fwrite(full, file.path(OUT_DIR, "bayesian_dependence_ranked_candidates_with_coverage.csv"))
fwrite(top25, file.path(OUT_DIR, "bayesian_dependence_top25_evidence.tsv"), sep = "\t")
fwrite(top10, file.path(OUT_DIR, "bayesian_dependence_top10_evidence.tsv"), sep = "\t")
fwrite(focus, file.path(OUT_DIR, "focus_gene_integrated_evidence_summary.tsv"), sep = "\t")
fwrite(summary_table, file.path(OUT_DIR, "bayesian_dependence_report_summary.tsv"), sep = "\t")
fwrite(dependence, file.path(OUT_DIR, "copula_dependence_summary.tsv"), sep = "\t")

cat("\n============================================================\n")
cat("BAYESIAN DEPENDENCE REPORT CREATED\n")
cat("============================================================\n\n")
cat("Output directory:\n", normalizePath(OUT_DIR), "\n\n", sep = "")
cat("Top 10 candidates:\n")
print(top10[, .(
  posterior_rank, gene_symbol, candidate_id, review_tier, coverage_class,
  posterior_rank_mean, p_rank1, p_top3, p_top5, rank_lower_envelope, rank_upper_envelope
)])
cat("\nfocus_gene summary:\n")
print(focus)

a11_finish(.ctx,OUT_DIR,"a11_coverage_report",list(producer_run_id=.producer$run_id,focus_gene=FOCUS_GENE))
