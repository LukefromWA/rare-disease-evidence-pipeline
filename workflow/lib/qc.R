# workflow/lib/qc.R -- the single QC definition (A1/A2), exclusion reasons, shortlist lanes, genomic channel
suppressPackageStartupMessages(library(data.table))

# Every flag is kept. dp_rule "gt" uses DP > min_dp; "ge" uses DP >= min_dp for sensitivity analysis. Returns a COPY.
apply_core_qc <- function(dt, qc) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  assert_columns(dt, c("candidate_id", "FILTER", "DP", "GQ", "AD_ref", "AD_alt", "MAX_AF", "IMPACT"), "QC input")
  num <- function(x) suppressWarnings(as.numeric(x)); dp_rule <- qc$dp_rule %||% "gt"; min_dp <- qc$min_dp %||% qc$dp_gt
  allowed_f <- unlist(qc$allowed_filters %||% list("PASS", ".")); allowed_i <- toupper(unlist(qc$allowed_impacts %||% list("HIGH", "MODERATE")))
  dt[, `:=`(DP = num(DP), GQ = num(GQ), AD_ref = num(AD_ref), AD_alt = num(AD_alt), MAX_AF = num(MAX_AF))]
  dt[, VAF := data.table::fifelse(is.finite(AD_ref + AD_alt) & (AD_ref + AD_alt) > 0, AD_alt / (AD_ref + AD_alt), NA_real_)]   # informative reads only
  dt[, pass_filter := toupper(FILTER) %in% toupper(allowed_f)]
  dt[, pass_dp := is.finite(DP) & (if (dp_rule == "gt") DP > min_dp else DP >= min_dp)]
  dt[, pass_gq := is.finite(GQ) & GQ >= qc$min_gq]
  dt[, pass_ad := is.finite(AD_alt) & AD_alt >= qc$min_alt_reads]
  dt[, af_status := data.table::fifelse(is.na(MAX_AF), "missing_af_kept", "af_present")]       # never invent AF; missing is flagged
  dt[, pass_af := is.na(MAX_AF) | MAX_AF <= qc$max_af]
  dt[, pass_impact := toupper(trimws(IMPACT)) %in% allowed_i]
  dt[, pass_core_qc := pass_filter & pass_dp & pass_gq & pass_ad & pass_af & pass_impact]
  if (anyDuplicated(dt$candidate_id)) stop("candidate_id must be unique after candidate construction", call. = FALSE)
  dt[]
}
# First failing reason in a declared order, so exclusion percentages sum to 100.
exclusion_reason <- function(dt) dt[, data.table::fcase(!pass_af, "gnomAD_AF", !pass_impact, "impact", !pass_filter, "caller_filter", !pass_dp, "DP", !pass_gq, "GQ", !pass_ad, "AD_alt", default = "pass")]

# A2 -- lanes, triage score; keeps every pass_core_qc row
build_shortlist <- function(em, qc) {
  em <- data.table::copy(em); em[, impact_norm := toupper(trimws(IMPACT))]
  sl <- em[pass_core_qc == TRUE]
  sl[, `:=`(is_high = impact_norm == "HIGH", is_mod = impact_norm == "MODERATE", rare = is.na(MAX_AF) | MAX_AF <= qc$max_af)]
  rng <- unlist(qc$mosaic_vaf_range %||% list(0.05, 0.35)); ma <- qc$mosaic_min_alt_reads %||% 4
  sl[, shortlist_lane := fcase(is.finite(VAF) & VAF >= rng[1] & VAF < rng[2] & AD_alt >= ma, "mosaic_review", is_high, "high_priority", default = "functional_rare")]
  sl[, triage_score := 4 * is_high + 2 * is_mod + 1 * rare + 0.25 * (GQ >= 50)]
  data.table::setorder(sl, -triage_score, gene_symbol, candidate_id); sl[, shortlist_rank := .I]
  stopifnot(nrow(sl) == sum(em$pass_core_qc))                         # stop rule: no pass_core_qc row omitted
  em[, exclusion := exclusion_reason(em)]
  list(shortlist = sl[], exclusion_summary = em[, .(n = .N, pct = round(100 * .N / nrow(em), 3)), by = exclusion][order(-n)], filter_audit = em)
}
genomic_channel <- function(sl) pmin(1, sl$triage_score / 5.25)        # declared transform (max triage = 5.25); compare with the recorded run

# Pair table within gene. Phase is UNKNOWN unless PGT + PS/PID support it; never asserted from VAF.
pair_candidates <- function(sl) {
  g <- sl[!is.na(gene_symbol), .N, by = gene_symbol][N >= 2, gene_symbol]; if (!length(g)) return(data.table(gene_symbol = character(), id1 = character(), id2 = character(), phase_status = character(), phase_evidence = character()))
  rbindlist(lapply(g, function(gn) { d <- sl[gene_symbol == gn]; cb <- utils::combn(nrow(d), 2)
    rbindlist(lapply(seq_len(ncol(cb)), function(k) { a <- d[cb[1, k]]; b <- d[cb[2, k]]
      has <- all(c("PGT", "PS") %in% names(d)) && !is.na(a$PS) && !is.na(b$PS) && a$PS == b$PS && !is.na(a$PGT) && !is.na(b$PGT) && a$PGT != "." && b$PGT != "."
      data.table(gene_symbol = gn, id1 = a$candidate_id, id2 = b$candidate_id,
                 phase_status = if (has) (if (identical(a$PGT, b$PGT)) "cis_by_read_backed_PGT_PS" else "trans_by_read_backed_PGT_PS") else "unknown",
                 phase_evidence = if (has) "PGT+PS same phase set" else "none (no parental or read-backed evidence)") })) }))
}
