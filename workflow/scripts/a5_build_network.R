#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A5.1: STRING and Reactome graphs (public resources; cached as one RDS)
a <- parse_args(c("config", "sample", "outdir")); cfg <- read_config(a$config)
sg <- load_string_graph(need_files(rp(cfg, cfg$resources$string_links), "STRING links"), need_files(rp(cfg, cfg$resources$string_info), "STRING info"), cfg$ppr$string_min_score)
rg <- load_reactome_graph(need_files(rp(cfg, cfg$resources$reactome_mapping), "Reactome mapping"), need_files(rp(cfg, cfg$resources$uniprot_gene_names), "UniProt gene names"), cfg$ppr$max_pathway_size)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE); saveRDS(list(string = sg$graph, reactome = rg$graph), file.path(od, "graphs.rds"), compress = "xz")
write_node_manifest(od, "a5_graphs", cfg, inputs = rp(cfg, c(cfg$resources$string_links, cfg$resources$string_info, cfg$resources$reactome_mapping, cfg$resources$uniprot_gene_names)), outputs = file.path(od, "graphs.rds"),
  file = "manifest_graphs.json", extra = list(string_min_score = cfg$ppr$string_min_score, string_nodes = igraph::vcount(sg$graph), string_edges = igraph::ecount(sg$graph), reactome_nodes = igraph::vcount(rg$graph),
  reactome_edges = igraph::ecount(rg$graph), string_file = basename(cfg$resources$string_links)))
