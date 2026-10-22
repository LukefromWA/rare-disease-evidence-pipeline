#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A5.2: random walks with restart (STRING, multiplex STRING+Reactome); informed and agnostic views; restart grid
a <- parse_args(c("config", "sample", "shortlist", "hpo", "graphs", "outdir")); cfg <- read_config(a$config)
sl <- read_gz(a$shortlist); hpo <- read_gz(a$hpo); gr <- readRDS(a$graphs); genes <- unique(sl$gene_symbol[!is.na(sl$gene_symbol)])
seeds <- hpo[is.finite(hpo_lin_bma)][order(-hpo_lin_bma, gene_symbol)][seq_len(min(cfg$ppr$phenotype_top_n, .N)), gene_symbol]
kn <- unlist(cfg$known_genes); layers <- list(string = gr$string, pathway = gr$reactome)
views_for <- function(r) {
  v <- run_network_views(gr$string, seeds, kn, r); mp <- multiplex_rwr(layers, seeds, r, cfg$ppr$coupling); mk_ <- multiplex_rwr(layers, kn, r, cfg$ppr$coupling)
  list(string = list(informed = consensus_view(v, "informed"), agnostic = consensus_view(v, "agnostic")),
       multiplex = list(informed = if (!is.null(mk_)) 0.5 * mp$score + 0.5 * mk_$score[names(mp$score)] else mp$score, agnostic = mp$score))
}
tab_for <- function(r) { nv <- views_for(r); out <- data.table(gene_symbol = genes)
  for (view in names(nv)) for (mode in c("informed", "agnostic")) { v <- nv[[view]][[mode]]
    set(out, j = paste0(view, "_", mode, "_rwr"), value = as.numeric(v[genes])); set(out, j = paste0(view, "_", mode, "_pct"), value = percentile_in_universe(as.numeric(v[genes]), v)) }
  out }
primary <- tab_for(cfg$ppr$restart_probability)
grid <- rbindlist(lapply(unlist(cfg$ppr$restart_grid), function(r) cbind(restart = r, tab_for(r)[, .(gene_symbol, string_agnostic_pct, string_informed_pct, multiplex_agnostic_pct, multiplex_informed_pct)])))
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(primary, file.path(od, "network_gene_scores.tsv")); write_table(grid, file.path(od, "rwr_grid.tsv"))
write_node_manifest(od, "a5", cfg, inputs = c(a$shortlist, a$hpo, a$graphs), outputs = file.path(od, c("network_gene_scores.tsv", "rwr_grid.tsv")),
  extra = list(restart_primary = cfg$ppr$restart_probability, restart_grid = as.list(unlist(cfg$ppr$restart_grid)), phenotype_seeds = as.list(seeds), known_seeds_in_graph = as.list(intersect(kn, igraph::V(gr$string)$name)),
               known_seeds_missing = as.list(setdiff(kn, igraph::V(gr$string)$name)), n_genes_unscorable_string = sum(is.na(primary$string_agnostic_rwr)), coupling = cfg$ppr$coupling))
