#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A2: shortlist, filter audit, exclusion summary, pair/phase table
a <- parse_args(c("config", "sample", "matrix", "outdir")); cfg <- read_config(a$config); em <- read_gz(a$matrix)
for (j in intersect(c("DP", "GQ", "AD_ref", "AD_alt", "MAX_AF", "VAF"), names(em))) set(em, j = j, value = as.numeric(em[[j]]))
for (j in grep("^pass_", names(em), value = TRUE)) set(em, j = j, value = as.logical(em[[j]]))
em[, pass_core_qc := as.logical(pass_core_qc)]
s <- build_shortlist(em, cfg$qc); od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
pr <- pair_candidates(s$shortlist)
write_table(s$shortlist, file.path(od, "shortlist.tsv.gz")); write_table(s$exclusion_summary, file.path(od, "exclusion_summary.tsv"))
write_table(s$filter_audit, file.path(od, "filter_audit.tsv.gz")); write_table(pr, file.path(od, "pair_candidates.tsv"))
stopifnot(!any(pr$phase_status %in% c("trans") & pr$phase_evidence == "none (no parental or read-backed evidence)"))
write_node_manifest(od, "a2", cfg, inputs = a$matrix, outputs = file.path(od, c("shortlist.tsv.gz", "exclusion_summary.tsv", "filter_audit.tsv.gz", "pair_candidates.tsv")),
  extra = list(n_shortlist = nrow(s$shortlist), lane_counts = as.list(table(s$shortlist$shortlist_lane)), n_pairs = nrow(pr), n_pairs_phase_unknown = sum(pr$phase_status == "unknown")))
logmsg("A2: shortlist ", nrow(s$shortlist), " rows | genes ", uniqueN(s$shortlist$gene_symbol), " | pairs ", nrow(pr), " (phase unknown preserved)")
