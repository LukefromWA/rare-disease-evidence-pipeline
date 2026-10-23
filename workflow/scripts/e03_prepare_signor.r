#!/usr/bin/env Rscript
.sf <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]), fixed = TRUE)
source(file.path(dirname(normalizePath(.sf, mustWork = TRUE)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("native", "native_sha256", "out", "audit", "summary", "source_url", "source_version", "license"))
for (p in c(a$out, paste0(a$out, ".resource_manifest.json"), a$audit, a$summary)) {
  if (file.exists(p)) stop("Refusing to overwrite: ", p)
}
if (!identical(tolower(trimws(a$native_sha256)), sha256_file(a$native))) stop("Native SIGNOR checksum mismatch")
x <- data.table::fread(a$native, sep = "\t", quote = "", colClasses = "character", na.strings = NULL, fill = FALSE, check.names = FALSE)
need_cols(x, c("ENTITYA", "TYPEA", "IDA", "DATABASEA", "ENTITYB", "TYPEB", "IDB", "DATABASEB", "EFFECT", "MECHANISM", "TAX_ID", "PMID", "DIRECT", "SCORE", "SIGNOR_ID"), "Native SIGNOR")
if (!nrow(x)) stop("Empty native SIGNOR table")
if (anyDuplicated(names(x))) stop("Duplicate column names")
native_columns <- names(x)
exact_duplicate <- duplicated(x[, ..native_columns])
trim <- function(v) trimws(as.character(v))
missing <- function(v) is.na(v) | trim(v) %in% c("", ".", "NA", "NaN")
ids_a <- trim(x$IDA); ids_b <- trim(x$IDB)
fx <- tolower(trim(x$EFFECT))
sign <- rep(NA_integer_, nrow(x))
sign[!is.na(fx) & grepl("^up-regulates($| )", fx)] <- 1L
sign[!is.na(fx) & grepl("^down-regulates($| )", fx)] <- -1L
effect_class <- data.table::fcase(
  !is.na(fx) & grepl("quantity", fx, fixed = TRUE), "quantity",
  !is.na(fx) & grepl("activity", fx, fixed = TRUE), "activity",
  !is.na(sign), "signed_unspecified",
  default = "unsigned"
)
human <- !missing(x$TAX_ID) & trim(x$TAX_ID) == "9606"
proteins <- tolower(trim(x$TYPEA)) == "protein" & tolower(trim(x$TYPEB)) == "protein"
uniprot <- toupper(trim(x$DATABASEA)) == "UNIPROT" & toupper(trim(x$DATABASEB)) == "UNIPROT"
valid_ids <- !missing(ids_a) & !missing(ids_b)
direct <- toupper(trim(x$DIRECT)); direct[missing(direct)] <- "UNK"
if (any(!direct %in% c("YES", "NO", "UNK"))) stop("Unexpected DIRECT values; review native vocabulary")
reason <- data.table::fcase(
  !human, "excluded_not_explicit_human_9606",
  !proteins, "excluded_not_protein_to_protein",
  !uniprot, "excluded_not_both_UNIPROT_identifiers",
  !valid_ids, "excluded_missing_identifier",
  is.na(sign), "excluded_no_directional_sign",
  exact_duplicate, "excluded_exact_native_duplicate",
  default = "kept"
)
aud <- data.table::data.table(
  native_row = seq_len(nrow(x)), native_signor_id = x$SIGNOR_ID,
  source_label = x$ENTITYA, target_label = x$ENTITYB,
  source_entity = ids_a, target_entity = ids_b,
  source_type = x$TYPEA, target_type = x$TYPEB,
  source_database = x$DATABASEA, target_database = x$DATABASEB,
  tax_id = x$TAX_ID, native_effect = x$EFFECT, sign = sign,
  effect_class = effect_class, mechanism = x$MECHANISM,
  direct = direct, native_score = x$SCORE, pmid = x$PMID,
  exact_native_duplicate = exact_duplicate, status = reason
)
keep <- which(reason == "kept")
if (!length(keep)) stop("No signed human protein-to-protein records survived; no outputs written")
prepared <- data.table::data.table(
  source_entity = ids_a[keep], target_entity = ids_b[keep],
  effect = ifelse(sign[keep] == 1L, "up-regulates", "down-regulates"),
  direct = direct[keep], confidence = rep(NA_real_, length(keep)), evidence_count = rep(NA_integer_, length(keep)),
  native_effect = x$EFFECT[keep], effect_class = effect_class[keep],
  mechanism = x$MECHANISM[keep], native_score = x$SCORE[keep],
  pmid = x$PMID[keep], native_signor_id = x$SIGNOR_ID[keep],
  native_source_row = keep, tax_id = rep("9606", length(keep)),
  source_label = x$ENTITYA[keep], target_label = x$ENTITYB[keep]
)
if (!all(prepared$effect %in% c("up-regulates", "down-regulates")) || !all(prepared$tax_id == "9606")) stop("Prepared contract failed")
summary <- aud[, .N, by = status]
if (sum(summary$N) != nrow(x)) stop("Audit row accounting failed")
score_text <- trim(x$SCORE)
score_num <- suppressWarnings(as.numeric(score_text))
score_diag <- list(finite_rows = sum(is.finite(score_num)), missing_rows = sum(missing(score_text)),
  nonmissing_nonfinite_rows = sum(!missing(score_text) & !is.finite(score_num)),
  minimum = if (any(is.finite(score_num))) min(score_num[is.finite(score_num)]) else NA_real_,
  maximum = if (any(is.finite(score_num))) max(score_num[is.finite(score_num)]) else NA_real_)
notes <- paste0(
  "Native SIGNOR sha256=", sha256_file(a$native), ". Explicit TAX_ID=9606, TYPEA/TYPEB=protein, DATABASEA/DATABASEB=UNIPROT; signed effects only; exact whole-native-row duplicates excluded with audit. ",
  "Canonical effect is up-regulates/down-regulates for existing helper compatibility; native_effect/effect_class retained. ",
  "Native SCORE preserved verbatim as native_score, NOT interpreted as an edge-existence probability. confidence and evidence_count intentionally unavailable. ",
  "Existing E11 edge_prob will substitute heuristic defaults for missing confidence/count; its weights are NOT calibrated posterior probabilities. DIRECT NO must be handled as indirect in the reviewed builder. ",
  "No accession remapping, complex/family expansion, self-loop removal, sign-conflict resolution, case-specific filtering, graph creation or propagation is performed here."
)
write_table(prepared, a$out); write_table(aud, a$audit)
jsonlite::write_json(list(
  adapter = "e03_prepare_signor_v1", native_sha256 = sha256_file(a$native),
  prepared_sha256 = sha256_file(a$out), source_version = a$source_version,
  n_native = nrow(x), n_explicit_human = sum(human), n_prepared = nrow(prepared),
  exact_native_duplicates_total = sum(exact_duplicate), sequential_exclusion_counts = summary,
  prepared_effect_classes = prepared[, .N, by = effect_class],
  prepared_directness = prepared[, .N, by = direct], native_score_diagnostics = score_diag,
  source_url = a$source_url, license = a$license, native_header = native_columns,
  notes = notes, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
), a$summary, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null")
chk <- read_gz(a$out, colClasses = "character", na.strings = NULL)
if (nrow(chk) != nrow(prepared) || !identical(names(chk), names(prepared))) stop("Prepared output read-back mismatch")
write_resource_manifest(a$out, "signor_human_protein_signed_prepared", a$source_url, a$source_version, a$license, notes = notes)
assert_frozen_resource(a$out, "signor_human_protein_signed_prepared")
cat("\nSIGNOR NORMALIZATION COMPLETE\n")
cat("Native rows:", nrow(x), "\nExplicit human rows before other filters:", sum(human), "\nPrepared signed protein-to-protein rows:", nrow(prepared), "\n")
print(summary)
cat("\nPrepared effect classes:\n"); print(prepared[, .N, by = effect_class])
cat("\nPrepared DIRECT values:\n"); print(prepared[, .N, by = direct])
cat("\nNative SCORE diagnostics (not probability calibration):\n"); print(score_diag)
cat("\nPrepared SHA256:", sha256_file(a$out), "\nOutputs:", a$out, a$audit, a$summary, "\n")
cat("No graph or readiness/enablement sentinel created. Review E11 directness and edge-weight assumptions before propagation.\n")
