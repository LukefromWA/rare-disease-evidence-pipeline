#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A1: evidence matrix from the A0 handoff VCF
a <- parse_args(c("config", "sample", "vcf", "outdir")); cfg <- read_config(a$config)
r <- build_evidence_matrix(a$vcf, cfg$qc, sample_id = cfg$vcf$sample_id, sample_index = cfg$vcf$sample_index, af_fields = unlist(cfg$qc$af_policy$fields), require_format = unlist(cfg$vcf$required_format_fields))
em <- r$matrix; od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
qs <- rbind(data.table(metric = c("n_vcf_records", "n_candidates", "n_pass_core_qc", "n_missing_af"), value = c(r$n_vcf_records, nrow(em), sum(em$pass_core_qc), sum(is.na(em$MAX_AF)))),
            em[, .(metric = paste0("n_fail_", c("filter", "dp", "gq", "ad", "af", "impact")), value = c(sum(!pass_filter), sum(!pass_dp), sum(!pass_gq), sum(!pass_ad), sum(!pass_af), sum(!pass_impact)))])
qs <- rbind(qs, data.table(metric = paste0("n_pass_core_qc_dp_", c("gt", "ge")), value = c(sum(apply_core_qc(em, modifyList(cfg$qc, list(dp_rule = "gt")))$pass_core_qc), sum(apply_core_qc(em, modifyList(cfg$qc, list(dp_rule = "ge")))$pass_core_qc))))
write_table(em, file.path(od, "evidence_matrix.tsv.gz")); write_table(qs, file.path(od, "qc_summary.tsv"))
stopifnot(nrow(em) == r$n_vcf_records, !anyDuplicated(em$candidate_id))
write_node_manifest(od, "a1", cfg, inputs = a$vcf, outputs = file.path(od, c("evidence_matrix.tsv.gz", "qc_summary.tsv")),
  extra = list(af_fields_used = as.list(r$af_fields_used), dp_rule = cfg$qc$dp_rule, n_pass_core_qc = sum(em$pass_core_qc), phase_fields_present = as.list(intersect(r$format_keys, c("PGT", "PID", "PS")))))
logmsg("A1: ", nrow(em), " candidates; ", sum(em$pass_core_qc), " pass core QC (dp_rule=", cfg$qc$dp_rule, ")")
