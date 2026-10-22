#!/usr/bin/env Rscript
# Evidence-preserving handoff. Does not fit, score, reweight, or modify frozen artifacts.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) %% 2L) stop("Use --a1 PATH --a2 PATH --a10 PATH --features PATH --predictions PATH --outdir NEW_DIRECTORY")
opt <- as.list(args[seq(2L, length(args), 2L)])
names(opt) <- sub("^--", "", args[seq(1L, length(args), 2L)])
required <- c("a1", "a2", "a10", "features", "predictions", "outdir")
if (!all(required %in% names(opt))) stop("Missing arguments: ", paste(setdiff(required, names(opt)), collapse = ", "))
for (pkg in c("data.table", "digest", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
}
suppressPackageStartupMessages(library(data.table))
sha <- function(p) digest::digest(file = p, algo = "sha256")
inputs <- unlist(opt[c("a1", "a2", "a10", "features", "predictions")], use.names = TRUE)
if (any(!file.exists(inputs))) stop("Missing input: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
if (file.exists(opt$outdir)) stop("Refusing to reuse output directory: ", opt$outdir)
source_sha <- vapply(inputs, sha, character(1))
sidecar <- paste0(opt$predictions, ".sha256")
if (!file.exists(sidecar) || !identical(trimws(readLines(sidecar, n = 1L)), sha(opt$predictions))) {
  stop("Missing or mismatched prediction checksum")
}
read_tab <- function(p) {
  if (grepl("[.]gz$", p, ignore.case = TRUE)) {
    fread(cmd = paste("gzip -cd", shQuote(p)), colClasses = "character", na.strings = NULL)
  } else fread(p, colClasses = "character", na.strings = NULL)
}
need <- function(x, cols, label) {
  if (!all(cols %in% names(x))) stop(label, " missing: ", paste(setdiff(cols, names(x)), collapse = ", "))
}
blank <- function(x) is.na(x) | trimws(x) %in% c("", ".", "NA", "NaN")
truth <- function(x) toupper(trimws(x)) %in% c("TRUE", "T", "1")
gene_tokens <- function(x) {
  if (is.na(x)) return(character())
  z <- trimws(unlist(strsplit(x, "[;,|&]")))
  unique(z[!blank(z)])
}
is_het <- function(x) {
  if (is.na(x) || !grepl("^[0-9]+[/|][0-9]+$", x)) return(FALSE)
  z <- as.integer(strsplit(x, "[/|]")[[1]])
  z[1] != z[2] && any(z > 0L)
}
a1 <- read_tab(opt$a1)
a2 <- read_tab(opt$a2)
a10 <- read_tab(opt$a10)
f <- read_tab(opt$features)
pred <- readRDS(opt$predictions)
s <- as.data.table(pred$summary)
need(a1, c("candidate_id", "gene_symbol", "GT", "pass_core_qc", "MAX_AF", "CLIN_SIG"), "A1")
need(a2, "candidate_id", "A2")
need(a10, c("gene_symbol", "best_candidate_id", "best_integrated_score", "best_IMPACT"), "A10")
need(f, c("case_id", "candidate_id", "gene", "G", "H", "N"), "C1 features")
need(s, c("case_id", "candidate_id", "probability", "rank_mean_probability"), "Prediction summary")
if (!nrow(a1) || !nrow(a10)) stop("Empty A1 or A10")
if (any(blank(a1$candidate_id)) || anyDuplicated(a1$candidate_id)) stop("A1 must have unique nonblank allele IDs; inspect duplicate/transcript rows first")
if (any(blank(a10$gene_symbol)) || any(blank(a10$best_candidate_id)) ||
    anyDuplicated(a10$gene_symbol) || anyDuplicated(a10$best_candidate_id)) stop("A10 must have unique genes and representatives")
if (anyDuplicated(f$candidate_id) || anyDuplicated(s$candidate_id)) stop("Expected one private case and one row per representative")
if (any(blank(f$case_id)) || anyNA(s$case_id) || uniqueN(f$case_id) != 1L || uniqueN(s$case_id) != 1L ||
    !identical(unique(f$case_id), as.character(unique(s$case_id)))) stop("C1/prediction case mismatch")
for (ids in list(f$candidate_id, as.character(s$candidate_id))) {
  if (length(ids) != nrow(a10) || !setequal(ids, a10$best_candidate_id)) stop("A10/C1/prediction candidate set mismatch")
}
mi <- match(a10$best_candidate_id, a1$candidate_id)
if (anyNA(mi)) stop("Some A10 representatives are missing from this A1 matrix")
if (!all(vapply(seq_len(nrow(a10)), function(i) a10$gene_symbol[i] %in% gene_tokens(a1$gene_symbol[mi[i]]), logical(1)))) stop("A10/A1 representative gene mismatch")
fi <- match(a10$best_candidate_id, f$candidate_id)
if (!identical(a10$gene_symbol, f$gene[fi])) stop("A10/C1 gene mismatch")
if (!all(a2$candidate_id %in% a1$candidate_id)) stop("A2 contains alleles missing from A1; source run mismatch")

# Preserve every A1 row, including failed-QC rows and genes outside the A10 set.
variants <- copy(a1)
variants[, handoff_source_row := .I]
variants[, handoff_in_A2 := candidate_id %in% a2$candidate_id]
variants[, handoff_A10_representative := candidate_id %in% a10$best_candidate_id]
variants[, handoff_heterozygous_GT := vapply(GT, is_het, logical(1))]
variants[, handoff_core_qc_pass := truth(pass_core_qc)]
variants[, handoff_core_qc_unknown := blank(pass_core_qc)]
variants[, handoff_AF_missing := blank(MAX_AF)]
variants[, handoff_clinical_annotation_missing := blank(CLIN_SIG)]
variants[, handoff_gene_tokens := lapply(gene_symbol, gene_tokens)]
links <- variants[, .(handoff_gene = unlist(handoff_gene_tokens)), by = handoff_source_row]
links <- unique(links[handoff_gene %in% a10$gene_symbol])
variants[, handoff_gene_tokens := NULL]
linked <- merge(links, variants, by = "handoff_source_row", all.x = TRUE, sort = FALSE)
if (anyNA(linked$candidate_id)) stop("Broken gene/variant link")
setorder(linked, handoff_gene, candidate_id)

counts <- linked[, .(
  n_A1_variants = uniqueN(candidate_id),
  n_A2_variants = uniqueN(candidate_id[handoff_in_A2]),
  n_core_qc_pass = uniqueN(candidate_id[handoff_core_qc_pass]),
  n_core_qc_unknown = uniqueN(candidate_id[handoff_core_qc_unknown]),
  n_heterozygous_calls = uniqueN(candidate_id[handoff_heterozygous_GT]),
  n_A2_qc_heterozygous_calls = uniqueN(candidate_id[handoff_in_A2 & handoff_core_qc_pass & handoff_heterozygous_GT]),
  n_AF_missing = uniqueN(candidate_id[handoff_AF_missing]),
  n_clinical_annotation_missing = uniqueN(candidate_id[handoff_clinical_annotation_missing])
), by = .(gene_symbol = handoff_gene)]
genes <- merge(copy(a10), counts, by = "gene_symbol", all.x = TRUE, sort = FALSE)
if (anyNA(genes$n_A1_variants)) stop("An A10 gene has no linked A1 evidence")
si <- match(genes$best_candidate_id, as.character(s$candidate_id))
for (nm in names(s)) genes[[paste0("stan_", nm)]] <- s[[nm]][si]
fi <- match(genes$best_candidate_id, f$candidate_id)
for (nm in intersect(c("G", "H", "N", "GxH"), names(f))) {
  genes[[paste0("c1_", nm)]] <- f[[nm]][fi]
}
genes[, multiple_het_calls_for_review := n_A2_qc_heterozygous_calls >= 2L]
genes[, phase_assessment := "not_assessed_by_this_adapter"]
genes[, genomic_evidence_in_frozen_Stan := FALSE]
setorder(genes, stan_rank_mean_probability, gene_symbol)
if (nrow(genes) != nrow(a10)) stop("Gene count changed")
if (!identical(genes$stan_probability, s$probability[match(genes$best_candidate_id, as.character(s$candidate_id))])) stop("Prediction values changed")
if (nrow(variants) != nrow(a1) || !identical(variants$candidate_id, a1$candidate_id)) stop("A1 evidence lost")
if (!identical(vapply(inputs, sha, character(1)), source_sha)) stop("An input changed while being read")

if (!dir.create(opt$outdir, recursive = TRUE)) stop("Could not create output directory")
outputs <- c(
  variants = file.path(opt$outdir, "variant_evidence_all.tsv"),
  gene_variant_links = file.path(opt$outdir, "candidate_gene_variant_evidence.tsv"),
  genes = file.path(opt$outdir, "gene_evidence_handoff.tsv")
)
fwrite(variants, outputs[["variants"]], sep = "\t", na = "")
fwrite(linked, outputs[["gene_variant_links"]], sep = "\t", na = "")
fwrite(genes, outputs[["genes"]], sep = "\t", na = "")
contract <- c(
  "Evidence-preserving handoff v1; NOT a revised statistical predictor.",
  "All A1 rows retained. All selected A10 genes retained. No top-k or new QC filter applied.",
  "A10 scores, C1 features and frozen Stan summaries remain separate fields.",
  "c1_G is the original gene prior, not variant-level pathogenicity evidence.",
  "Multiple heterozygous calls flag a review queue, not confirmed compound heterozygosity.",
  "Phase is not inferred, including from allele balance or GT order.",
  "Missing AF is not zero AF. Missing clinical annotation is not benign or pathogenic.",
  "Clinical annotation strings are preserved, not automatically adjudicated.",
  "No fit, new probabilities, revised ranking, calibration, or source edits were performed.",
  "The upstream best_integrated_score is preserved, not interpreted as a clinical probability.",
  "Input hashes document this run; they do not establish upstream biological correctness."
)
writeLines(contract, file.path(opt$outdir, "HANDOFF_CONTRACT.txt"))
manifest <- list(
  adapter_version = "evidence_handoff_v1",
  time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = lapply(names(inputs), function(k) list(role = k, path = normalizePath(inputs[[k]]), sha256 = unname(source_sha[[k]]))),
  outputs = lapply(names(outputs), function(k) list(role = k, path = normalizePath(outputs[[k]]), sha256 = sha(outputs[[k]]))),
  counts = list(A1_rows = nrow(a1), A10_genes = nrow(a10), linked_gene_variant_rows = nrow(linked), multiple_het_review_genes = sum(genes$multiple_het_calls_for_review)),
  predictor_id = pred$predictor_id,
  candidate_protocol_id = pred$candidate_protocol_id,
  original_prediction_interpretation = pred$interpretation,
  provisional_divergences = pred$provisional_divergences,
  source_files_unchanged = TRUE,
  predictor_modified = FALSE,
  phase_inferred = FALSE,
  new_probabilities_computed = FALSE,
  contract = contract,
  session_info = capture.output(sessionInfo())
)
jsonlite::write_json(manifest, file.path(opt$outdir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("\nHANDOFF COMPLETE\n")
cat("A1 rows preserved:", nrow(a1), "\n")
cat("A10 genes preserved:", nrow(genes), "\n")
cat("Gene/variant evidence rows:", nrow(linked), "\n")
cat("Genes with >=2 shortlisted QC-passing heterozygous calls:", sum(genes$multiple_het_calls_for_review), "\n")
cat("Frozen probabilities unchanged; phase not inferred.\n")
cat("Output directory:", opt$outdir, "\n")
