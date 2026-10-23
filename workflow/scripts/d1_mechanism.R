#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D1: mechanism characterization (rule-based). No external mechanism evidence is configured, so p_m^(1) = p_m^(0) and that is recorded.
a <- parse_args(c("config", "sample", "scenarios", "outdir")); cfg <- read_config(a$config); sc <- read_gz(a$scenarios)
sc[, `:=`(direction = direction_of(mechanism), p_m_updated = p_m, mechanism_evidence = "none_configured")]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(sc, file.path(od, "mechanism.tsv.gz"))
write_node_manifest(od, "d1", cfg, inputs = a$scenarios, outputs = file.path(od, "mechanism.tsv.gz"), extra = list(mechanism_evidence_sources = "none_configured", p_m_updated_equals_prior = TRUE))
