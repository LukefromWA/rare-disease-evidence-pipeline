#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A7: coverage-aware integration. posterior_informed uses the known-gene view and soft prior; posterior_agnostic uses neither.
a <- parse_args(c("config", "sample", "shortlist", "hpo", "network", "outdir")); cfg <- read_config(a$config)
sl <- read_gz(a$shortlist); hpo <- read_gz(a$hpo); nw <- read_gz(a$network)
sl <- merge(sl, hpo[, .(gene_symbol, hpo_lin_bma, hpo_lin_bma_pct, n_gene_terms, direct_overlap, annotation_status)], by = "gene_symbol", all.x = TRUE, sort = FALSE)
sl <- merge(sl, nw, by = "gene_symbol", all.x = TRUE, sort = FALSE); data.table::setorder(sl, shortlist_rank); sl[, genomic := genomic_channel(sl)]
W <- unlist(cfg$integration$channel_weights); ab <- cfg$integration$alpha_boost; kn <- unlist(cfg$known_genes)
if (abs(sum(W) - 1) > 1e-8) stop("channel weights must sum to 1", call. = FALSE)
chan <- function(mode) data.table(candidate_id = sl$candidate_id, gene = sl$gene_symbol, genomic = sl$genomic, hpo = sl$hpo_lin_bma_pct, string = sl[[paste0("string_", mode, "_pct")]], multiplex = sl[[paste0("multiplex_", mode, "_pct")]])
ii <- integrate_tier1(chan("informed"), W, kn, ab); ia <- integrate_tier1(chan("agnostic"), W, kn, ab)
sl[, `:=`(S_bar_informed = ii$S_bar, n_views = ii$n_views, posterior_informed = ii$posterior_informed, rank_informed = ii$rank_informed,
          S_bar_agnostic = ia$S_bar, posterior_agnostic = ia$posterior_agnostic, rank_agnostic = ia$rank_agnostic)]
stopifnot(abs(sum(sl$posterior_informed) - 1) < 1e-8, abs(sum(sl$posterior_agnostic) - 1) < 1e-8)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); write_table(sl, file.path(od, "integrated_scores.tsv.gz"))
write_node_manifest(od, "a7", cfg, inputs = c(a$shortlist, a$hpo, a$network), outputs = file.path(od, "integrated_scores.tsv.gz"), extra = list(channel_weights = as.list(W), alpha_boost = ab, known_genes = as.list(kn), modes = list("informed", "agnostic")))
