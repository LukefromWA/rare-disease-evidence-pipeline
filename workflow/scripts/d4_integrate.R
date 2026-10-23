#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# D4: coverage-aware integration across evidence families. Pathway evidence comes from D3.3 only (never a copy of mechanism evidence).
a <- parse_args(c("config", "sample", "targets", "method1", "method3", "method6", "outdir")); cfg <- read_config(a$config)
tu <- read_gz(a$targets); m1 <- read_gz(a$method1); m3 <- read_gz(a$method3); m6 <- read_gz(a$method6); fw <- unlist(cfg$track2$family_weights)
d <- merge(merge(merge(tu, m1[, .(scenario_id, drug_id, direction_gate)], by = c("scenario_id", "drug_id")), m3[, .(scenario_id, drug_id, pathway_jaccard)], by = c("scenario_id", "drug_id"), all.x = TRUE),
           m6[, .(scenario_id, drug_id, rescue_mean, p_rescue_gt0, propagation_status)], by = c("scenario_id", "drug_id"), all.x = TRUE)
d[, `:=`(e_mechanism = fcase(direction_gate == "pass", 0.90, direction_gate == "uncertain", 0.50, default = 0.10), e_pathway = normalize01(pathway_jaccard), e_propagation = p_rescue_gt0, e_signature = NA_real_, e_phenotype = NA_real_, e_clinical = NA_real_)]
E <- as.matrix(d[, paste0("e_", names(fw)), with = FALSE]); colnames(E) <- names(fw); agg <- weighted_available_mean(E, fw); d[, `:=`(evidence_score = as.numeric(agg), n_families_available = attr(agg, "n_views"))]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(d, file.path(od, "evidence.tsv.gz"))
write_node_manifest(od, "d4", cfg, inputs = c(a$targets, a$method1, a$method3, a$method6), outputs = file.path(od, "evidence.tsv.gz"), extra = list(family_weights = as.list(fw), families_unavailable = list("signature", "phenotype", "clinical")))
