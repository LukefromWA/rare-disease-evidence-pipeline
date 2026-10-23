#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: make_a11_input.R <A10.tsv.gz> <output.csv>")
}

x <- fread(args[1])
for (nm in intersect(c(
  "best_genomic_score", "best_HPO_semantic_score",
  "best_PPI_network_score", "best_multiplex_score"
), names(x))) {
  raw <- x[[nm]]
  value <- suppressWarnings(as.numeric(raw))
  if (any(!is.na(raw) & is.na(value)) ||
      any(is.nan(value) | is.infinite(value))) {
    stop(nm, " contains invalid numeric or nonfinite values.")
  }
}

required <- c(
  "best_candidate_id", "gene_symbol",
  "best_genomic_score", "best_HPO_semantic_score",
  "best_PPI_network_score", "best_multiplex_score"
)

missing <- setdiff(required, names(x))
if (length(missing)) {
  stop("Missing A10 columns: ", paste(missing, collapse = ", "))
}

out <- x[, .(
  candidate_id = best_candidate_id,
  gene_symbol,
  genomic = as.numeric(best_genomic_score),
  hpo = as.numeric(best_HPO_semantic_score),
  string = as.numeric(best_PPI_network_score),
  multiplex = as.numeric(best_multiplex_score)
)]

if (!nrow(out)) stop("No candidates in A10.")
if (anyNA(out$candidate_id) ||
    any(!nzchar(trimws(as.character(out$candidate_id)))) ||
    anyDuplicated(out$candidate_id)) {
  stop("Candidate IDs must be present and unique.")
}

for (nm in c("genomic", "hpo", "string", "multiplex")) {
  v <- out[[nm]]
  if (any(is.finite(v) & (v < 0 | v > 1))) {
    stop(nm, " contains values outside [0,1].")
  }
}

if (anyNA(out$gene_symbol) ||
    any(!nzchar(trimws(as.character(out$gene_symbol))))) {
  stop("Gene symbols must be present and nonblank.")
}
for (nm in c("genomic","hpo","string","multiplex")) {
  data.table::set(out, j=paste0("has_",nm),
                  value=is.finite(out[[nm]]))
}
if (file.exists(args[2])) stop("Output already exists; refusing to overwrite.")
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
fwrite(out, args[2])
cat("A11 input:", nrow(out), "best-per-gene candidates written.\n")
