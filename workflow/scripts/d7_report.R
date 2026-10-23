#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D7: final Track 2 table: categorical tiers with uncertainty, limitations, and the non-clinical disclaimer on every row
a <- parse_args(c("config", "sample", "gated", "sensitivity", "outdir")); cfg <- read_config(a$config); d <- read_gz(a$gated); s <- read_gz(a$sensitivity); g <- d[eligibility == "eligible_for_research_tiering"]; od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
if (!nrow(g) || !"drug_id" %in% names(s)) { rep <- data.table(note = "no drug is currently eligible for research tiering; inspect held/excluded rows", report_disclaimer = "Computational prioritization only; not evidence of clinical efficacy or a treatment recommendation.") } else {
  rep <- g[, .(best_scenario = scenario_id[which.max(evidence_score)], gene = gene[which.max(evidence_score)], mechanism = mechanism[which.max(evidence_score)], evidence_score = max(evidence_score, na.rm = TRUE),
               n_families_available = max(n_families_available), eligibility = eligibility[1], drug_name = drug_name[1]), by = drug_id]
  rep <- merge(rep, s, by = "drug_id")
  rep[, final_tier := fcase(eligibility == "eligible_for_research_tiering" & top_k_prob >= 0.5 & n_families_available >= 2, "high-priority computational hypothesis", eligibility == "eligible_for_research_tiering", "secondary hypothesis", default = "exploratory hypothesis")]
  rep[, report_disclaimer := "Computational prioritization only; not evidence of clinical efficacy or a treatment recommendation."]; setorder(rep, -top_k_prob, -evidence_score) }
write_table(rep, file.path(od, "final_report.tsv")); write_node_manifest(od, "d7", cfg, inputs = c(a$gated, a$sensitivity), outputs = file.path(od, "final_report.tsv"), extra = list(n_rows = nrow(rep)))
