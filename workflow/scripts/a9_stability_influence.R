#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A9: seeded perturbation reranks (top-m frequency, Kuncheva) and leave-one-channel-out influence. Runs only after integration.
a <- parse_args(c("config", "sample", "integrated", "outdir")); cfg <- read_config(a$config); sl <- read_gz(a$integrated)
W <- unlist(cfg$integration$channel_weights); sc <- as.matrix(sl[, .(genomic = genomic, hpo = hpo_lin_bma_pct, string = string_informed_pct, multiplex = multiplex_informed_pct)]); sc[!is.finite(sc)] <- NA
st <- run_stability(sc, sl$candidate_id, W, B = cfg$tier15$B, top_m = cfg$tier15$top_m, kappa = cfg$tier15$kappa %||% 100, seed = cfg$run$seed)
loco <- leave_one_channel_out(sc, W); loco[, candidate_id := sl$candidate_id[candidate]]
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(st$top_m_frequency, file.path(od, "top_m_frequency.tsv")); write_table(loco, file.path(od, "leave_one_channel_out.tsv.gz"))
write_table(data.table(metric = c("kuncheva_mean_vs_baseline", "B", "top_m"), value = c(st$kuncheva_mean, cfg$tier15$B, cfg$tier15$top_m)), file.path(od, "stability_summary.tsv"))
write_node_manifest(od, "a9", cfg, inputs = a$integrated, outputs = file.path(od, c("top_m_frequency.tsv", "leave_one_channel_out.tsv.gz", "stability_summary.tsv")), extra = list(seed_set_once = cfg$run$seed, B = cfg$tier15$B))
