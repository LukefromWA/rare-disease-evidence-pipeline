# workflow/lib/csq.R -- VCF/CSQ parsing, FORMAT fields (incl. PGT/PID/PS), transcript selection, AF policy, A1 evidence matrix
suppressPackageStartupMessages(library(data.table))

read_vcf_header <- function(path) {
  need_files(path, "VCF"); con <- gzfile(path, "rt"); on.exit(close(con)); hdr <- character()
  repeat { l <- readLines(con, n = 1L); if (!length(l)) stop("no #CHROM line in ", path); hdr <- c(hdr, l); if (startsWith(l, "#CHROM")) break }
  hdr
}
read_vcf_records <- function(path) {
  dt <- data.table::fread(cmd = paste("gzip -cd", shQuote(path), "| grep -v '^##'"), sep = "\t", header = TRUE, colClasses = "character", quote = "", na.strings = "")
  data.table::setnames(dt, 1L, "CHROM"); dt
}
read_csq_schema_line <- function(line) {
  if (!grepl("Format: ", line, fixed = TRUE)) stop("Invalid CSQ header")
  make.unique(strsplit(sub('.*Format: ([^"]*).*', "\\1", line), "|", fixed = TRUE)[[1]])
}
csq_schema_from_header <- function(hdr) { l <- grep("ID=CSQ", hdr, value = TRUE); if (!length(l)) stop("no CSQ INFO header: was the VCF annotated with VEP?", call. = FALSE); read_csq_schema_line(l[1]) }
info_csq <- function(info) { out <- rep(NA_character_, length(info)); has <- grepl("(^|;)CSQ=", info); out[has] <- sub("^;?CSQ=", "", regmatches(info, regexpr("(^|;)CSQ=[^;]*", info))); out }

parse_ad <- function(x) {                       # "21,25" -> ref 21, alt 25 (needs bcftools norm -m -any first)
  pieces <- data.table::tstrsplit(x, ",", fixed = TRUE, fill = NA_character_)
  if (length(pieces) != 2L) stop("AD has != 2 values: run bcftools norm -m -any before A1", call. = FALSE)
  list(ref = as.numeric(pieces[[1L]]), alt = as.numeric(pieces[[2L]]))
}
FORMAT_KEYS <- c("GT", "DP", "GQ", "AD", "PGT", "PID", "PS")
parse_format_fields <- function(format, sample) {
  n <- length(format); out <- data.table(GT = rep(NA_character_, n), DP = NA_real_, GQ = NA_real_, AD = NA_character_, PGT = NA_character_, PID = NA_character_, PS = NA_character_)
  for (f in unique(format)) { keys <- strsplit(f, ":", fixed = TRUE)[[1]]; idx <- which(format == f)
    parts <- data.table::tstrsplit(sample[idx], ":", fixed = TRUE, fill = NA_character_)
    for (k in intersect(FORMAT_KEYS, keys)) { j <- match(k, keys); if (j > length(parts)) next; v <- parts[[j]]; if (k %in% c("DP", "GQ")) v <- suppressWarnings(as.numeric(v)); data.table::set(out, idx, k, v) } }
  out
}
# One row per (candidate, transcript), CSQ fields as columns
parse_csq_annotations <- function(dt, schema) {
  assert_columns(dt, c("candidate_id", "CSQ"), "VCF table")
  long <- dt[, .(csq = unlist(strsplit(CSQ, ",", fixed = TRUE))), by = candidate_id][!is.na(csq)]
  parts <- data.table::tstrsplit(long$csq, "|", fixed = TRUE, fill = "")
  if (length(parts) < length(schema)) parts <- c(parts, rep(list(rep("", nrow(long))), length(schema) - length(parts)))
  if (length(parts) > length(schema)) stop("CSQ has more fields than the header schema")
  names(parts) <- schema; out <- cbind(long[, .(candidate_id)], data.table::as.data.table(parts))
  for (col in schema) data.table::set(out, which(missing_value(out[[col]])), col, NA_character_)
  out <- merge(out, out[, .(n_transcripts = .N), by = candidate_id], by = "candidate_id")
  merge(out, dt[, !"CSQ"], by = "candidate_id", all.x = TRUE)
}
IMPACT_LEVELS <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
# MANE Select, then canonical, then worst impact, then transcript id (deterministic)
select_representative_transcript <- function(long) {
  long <- data.table::copy(long)
  long[, impact_rank := fifelse(is.na(match(IMPACT, IMPACT_LEVELS)), 5L, match(IMPACT, IMPACT_LEVELS))]
  long[, `:=`(mane_flag = as.integer(!is.na(MANE_SELECT)), canon_flag = as.integer(CANONICAL %in% "YES"))]
  long[, best_impact_any := IMPACT_LEVELS[min(impact_rank, 4L)], by = candidate_id]
  data.table::setorder(long, candidate_id, -mane_flag, -canon_flag, impact_rank, Feature, na.last = TRUE)
  out <- long[, .SD[1L], by = candidate_id]; if (anyDuplicated(out$candidate_id)) stop("transcript selection left duplicates"); out[]
}
# Named AF policy: use only fields listed in the config; record the source and parse status; never invent AF
max_af_from_csq <- function(rep_dt, schema, fields = NULL) {
  cand <- if (is.null(fields)) schema[grepl("(^|_)AF$|^MAX_AF$", schema) & schema != "AF" | schema == "MAX_AF"] else fields
  af_cols <- intersect(unique(cand), names(rep_dt)); rep_dt
  if (!length(af_cols)) return(list(value = rep(NA_real_, nrow(rep_dt)), source = rep(NA_character_, nrow(rep_dt)), cols = character()))
  m <- vapply(af_cols, function(cn) vapply(strsplit(as.character(rep_dt[[cn]]), "&", fixed = TRUE), function(z) { z <- suppressWarnings(as.numeric(z)); if (all(is.na(z))) NA_real_ else max(z, na.rm = TRUE) }, 0), numeric(nrow(rep_dt)))
  m <- matrix(m, nrow = nrow(rep_dt), dimnames = list(NULL, af_cols))
  list(value = apply(m, 1, function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)),
       source = apply(m, 1, function(r) if (all(is.na(r))) NA_character_ else af_cols[which.max(ifelse(is.na(r), -Inf, r))]), cols = af_cols)
}

# A1 -- evidence matrix from the A0 handoff VCF (biallelic, VEP-annotated)
build_evidence_matrix <- function(vcf, qc, sample_id = NULL, sample_index = 1L, af_fields = NULL, require_format = c("GT", "DP", "GQ", "AD")) {
  hdr <- read_vcf_header(vcf); schema <- csq_schema_from_header(hdr); rec <- read_vcf_records(vcf)
  need_cols(rec, c("CHROM", "POS", "REF", "ALT", "FILTER", "INFO", "FORMAT"), "VCF records")
  if (any(grepl(",", rec$ALT, fixed = TRUE))) stop("multi-allelic ALT present: A0 must hand off a biallelic VCF", call. = FALSE)
  samp_col <- if (!is.null(sample_id) && sample_id %in% names(rec)) sample_id else names(rec)[9L + sample_index]
  rec[, candidate_id := make_candidate_id(CHROM, POS, REF, ALT)]; if (anyDuplicated(rec$candidate_id)) stop("duplicate candidate ids in VCF", call. = FALSE)
  keys_present <- unique(unlist(strsplit(rec$FORMAT, ":", fixed = TRUE))); miss <- setdiff(require_format, keys_present)
  if (length(miss)) stop("required FORMAT field(s) absent: ", paste(miss, collapse = ", "), call. = FALSE)
  gt <- parse_format_fields(rec$FORMAT, rec[[samp_col]]); ad <- parse_ad(fifelse(is.na(gt$AD) | gt$AD == ".", "NA,NA", gt$AD))
  base <- cbind(rec[, .(candidate_id, CHROM, POS = as.integer(POS), REF, ALT, FILTER, CSQ = info_csq(INFO))], gt[, .(GT, DP, GQ, PGT, PID, PS)], AD_ref = ad$ref, AD_alt = ad$alt)
  rep <- select_representative_transcript(parse_csq_annotations(base, schema))
  need_cols(rep, c("SYMBOL", "Consequence", "IMPACT", "Feature"), "CSQ fields (check VEP flags)")
  af <- max_af_from_csq(rep, schema, af_fields); rep[, `:=`(MAX_AF = af$value, af_source = af$source)]
  if (!length(af$cols)) warning("no allele-frequency fields found in CSQ: every row is flagged missing_af_kept. Re-run VEP with --max_af / gnomAD fields.", call. = FALSE)
  data.table::setnames(rep, c("SYMBOL", "Consequence"), c("gene_symbol", "consequence"))
  keep <- intersect(c("candidate_id", "CHROM", "POS", "REF", "ALT", "FILTER", "GT", "DP", "GQ", "PGT", "PID", "PS", "AD_ref", "AD_alt", "MAX_AF", "af_source", "gene_symbol", "consequence",
                      "IMPACT", "Feature", "HGVSc", "HGVSp", "CLIN_SIG", "n_transcripts", "best_impact_any"), names(rep))
  em <- apply_core_qc(rep[, ..keep], qc); stopifnot(nrow(em) == nrow(base))
  list(matrix = em[], schema = schema, af_fields_used = af$cols, n_vcf_records = nrow(rec), format_keys = keys_present)
}
