# 01_helpers_qc.R  -- shared helpers and the single QC definition (tested)
suppressPackageStartupMessages(library(data.table))

assert_columns <- function(x, required, label = "table") {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " missing column(s): ", paste(missing, collapse = ", "))
  invisible(x)
}

# One candidate id form per table; "GRCh38:" is added only when joining.
make_candidate_id <- function(chrom, pos, ref, alt)
  paste0(sub("^chr", "", chrom), ":", pos, ":", ref, ":", alt)

normalize_candidate_id <- function(x) {
  x <- trimws(as.character(x)); x[x == ""] <- NA_character_
  has_build <- grepl("^GRCh38:", x, ignore.case = TRUE)
  x[!is.na(x) & !has_build] <- paste0("GRCh38:", x[!is.na(x) & !has_build])
  x
}

missing_value <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA", ".", "-")

# Core QC. Every flag is kept, not only the composite. Returns a COPY.
apply_core_qc <- function(dt, qc) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  assert_columns(dt, c("candidate_id", "FILTER", "DP", "GQ", "AD_ref", "AD_alt", "MAX_AF", "IMPACT"), "QC input")
  num <- function(x) suppressWarnings(as.numeric(x))
  dt[, `:=`(DP = num(DP), GQ = num(GQ), AD_ref = num(AD_ref), AD_alt = num(AD_alt), MAX_AF = num(MAX_AF))]
  dt[, VAF := data.table::fifelse(is.finite(AD_ref + AD_alt) & (AD_ref + AD_alt) > 0,
                                  AD_alt / (AD_ref + AD_alt), NA_real_)]   # informative reads only
  dt[, pass_filter := toupper(FILTER) %in% c("PASS", ".")]
  dt[, pass_dp     := is.finite(DP) & DP >  qc$dp_gt]
  dt[, pass_gq     := is.finite(GQ) & GQ >= qc$min_gq]
  dt[, pass_ad     := is.finite(AD_alt) & AD_alt >= qc$min_alt_reads]
  dt[, pass_af     := is.na(MAX_AF) | MAX_AF <= qc$max_af]
  dt[, pass_impact := toupper(trimws(IMPACT)) %in% c("HIGH", "MODERATE")]
  dt[, pass_core_qc := pass_filter & pass_dp & pass_gq & pass_ad & pass_af & pass_impact]
  if (anyDuplicated(dt$candidate_id)) stop("candidate_id must be unique after candidate construction")
  dt[]
}

# First failing reason, in a declared order, so exclusion percentages sum to 100.
exclusion_reason <- function(dt) {
  dt[, data.table::fcase(
    !pass_af, "gnomAD_AF", !pass_impact, "impact", !pass_filter, "caller_filter",
    !pass_dp, "DP", !pass_gq, "GQ", !pass_ad, "AD_alt", default = "pass")]
}

# Weights matched BY NAME. Missing evidence = missing coverage, never zero.
weighted_available_mean <- function(scores, weights) {
  scores <- as.matrix(scores)
  stopifnot(!is.null(colnames(scores)), all(colnames(scores) %in% names(weights)))
  w <- weights[colnames(scores)]
  avail <- is.finite(scores)
  num <- rowSums(ifelse(avail, sweep(scores, 2L, w, "*"), 0))
  den <- rowSums(sweep(avail * 1, 2L, w, "*"))
  out <- num / den
  out[den == 0] <- NA_real_
  attr(out, "n_views") <- rowSums(avail)
  out
}
