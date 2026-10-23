#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# E02: convert ONE native UniProtKB TSV export into the prepared table E10 expects.
# No gene recovery from synonyms, no identifier guessing, no edits to the native file.
a <- parse_args(c("native", "native_sha256", "out", "audit", "source_url", "source_version", "license"))
if (file.exists(a$out) || file.exists(paste0(a$out, ".resource_manifest.json")) || file.exists(a$audit)) stop("Refusing to overwrite prepared output, its manifest, or the audit")
if (!identical(tolower(trimws(a$native_sha256)), sha256_file(a$native))) stop("Native snapshot checksum differs from the supplied value")
x <- data.table::fread(a$native, sep = "\t", colClasses = "character", na.strings = NULL, quote = "")
want <- c("Entry", "Gene Names (primary)", "Gene Names", "Reviewed", "Organism (ID)", "Date of last modification")
if (!identical(names(x), want)) stop("Unexpected native header: ", paste(names(x), collapse = " | "))
data.table::setnames(x, c("entry", "gene_primary", "gene_names", "reviewed_raw", "organism_id", "date_modified"))
for (j in names(x)) data.table::set(x, j = j, value = trimws(x[[j]]))
if (!nrow(x)) stop("Empty native table")
if (any(!nzchar(x$entry)) || anyDuplicated(x$entry)) stop("Entry must be unique and nonblank")
if (!all(x$organism_id == "9606")) stop("Non-human rows present")
if (!all(tolower(x$reviewed_raw) == "reviewed")) stop("Unreviewed rows present")
status <- ifelse(!nzchar(x$gene_primary), "no_primary_gene_name",
          ifelse(grepl("[;,[:space:]]", x$gene_primary), "ambiguous_primary_gene_name", "kept"))
data.table::set(x, j = "status", value = status)
kept <- x[x$status == "kept"]
if (!nrow(kept)) stop("No usable rows after normalization")
prepared <- data.table::data.table(uniprot_id = kept$entry, gene_symbol = kept$gene_primary, reviewed = "TRUE")
sym_n <- table(kept$gene_primary)
audit <- data.table::data.table(entry = x$entry, gene_primary = x$gene_primary, gene_names = x$gene_names, status = x$status,
  kept_rows_sharing_symbol = ifelse(x$status == "kept", as.integer(sym_n[x$gene_primary]), NA_integer_))
write_table(prepared, a$out); write_table(audit, a$audit)
chk <- read_gz(a$out, colClasses = "character")
if (nrow(chk) != nrow(prepared) || !identical(names(chk), c("uniprot_id", "gene_symbol", "reviewed")) || !all(chk$reviewed == "TRUE")) stop("Prepared table failed read-back check")
n_none <- sum(x$status == "no_primary_gene_name"); n_amb <- sum(x$status == "ambiguous_primary_gene_name")
shared <- sym_n[sym_n > 1L]
notes <- sprintf("Prepared from native UniProtKB TSV sha256=%s. Kept rows with one nonblank primary gene name; reviewed set to TRUE for reviewed entries; no synonym-based gene recovery. Rows: native=%d kept=%d dropped_no_primary=%d dropped_ambiguous=%d. Symbols shared by >1 kept entry: %d.",
  sha256_file(a$native), nrow(x), nrow(kept), n_none, n_amb, length(shared))
write_resource_manifest(a$out, "uniprot_reviewed_human_prepared", a$source_url, a$source_version, a$license, notes = notes)
assert_frozen_resource(a$out, "uniprot_reviewed_human_prepared")
cat("\nUNIPROT NORMALIZATION COMPLETE\n")
cat("Native rows:", nrow(x), "\nKept:", nrow(kept), "\nDropped, no primary gene name:", n_none, "\nDropped, ambiguous primary gene name:", n_amb, "\n")
cat("Symbols shared by more than one kept entry:", length(shared), "(entries involved:", sum(shared), ")\n")
cat("Prepared:", a$out, "\nAudit:", a$audit, "\n")
