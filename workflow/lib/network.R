# workflow/lib/network.R -- A5: STRING graph, Reactome pathway layer, sparse and multiplex random walks with restart
suppressPackageStartupMessages({ library(data.table); library(igraph); library(Matrix) })

# edges: data.table(a, b, combined_score) with gene symbols; links files list BOTH directions
build_string_graph <- function(edges, min_score = 700) {
  e <- data.table::as.data.table(edges)[combined_score >= min_score & a != b]
  e <- e[, .(weight = max(combined_score) / 1000), by = .(a = pmin(a, b), b = pmax(a, b))]   # undirected, no doubling
  igraph::graph_from_data_frame(e[, .(from = a, to = b, weight)], directed = FALSE)
}

# restart = probability of jumping back to the seed (guide: alpha = 0.85; igraph damping = 1 - restart)
weighted_rwr <- function(graph, seed_nodes, restart_probability, iterations = 100L, tol = 1e-10) {
  stopifnot(is.numeric(restart_probability), restart_probability > 0, restart_probability < 1)
  nodes <- igraph::V(graph)$name
  seed_nodes <- intersect(unique(seed_nodes), nodes)
  if (!length(seed_nodes)) return(NULL)                  # unscorable: caller records M = 0, never zeros
  A  <- igraph::as_adjacency_matrix(graph, attr = "weight", sparse = TRUE)
  cs <- Matrix::colSums(A)
  P  <- A %*% Matrix::Diagonal(x = ifelse(cs > 0, 1 / cs, 0))
  dangling <- which(cs == 0)
  seed <- stats::setNames(numeric(length(nodes)), nodes); seed[seed_nodes] <- 1 / length(seed_nodes)
  score <- seed
  for (i in seq_len(iterations)) {
    nxt <- restart_probability * seed +
      (1 - restart_probability) * (as.numeric(P %*% score) + sum(score[dangling]) * seed)
    done <- sum(abs(nxt - score)) < tol
    score <- nxt
    if (done) break
  }
  stats::setNames(score, nodes)
}

# Views: phenotype always; known-gene view only in "informed" mode.
run_network_views <- function(graph, phenotype_seeds, known_genes, restart_probability, modes = c("informed", "agnostic")) {
  views <- list(phenotype = weighted_rwr(graph, phenotype_seeds, restart_probability))
  if ("informed" %in% modes && length(known_genes))
    views$known <- weighted_rwr(graph, known_genes, restart_probability)
  views
}
consensus_view <- function(views, mode) {
  if (mode == "informed" && !is.null(views$known)) 0.5 * views$phenotype + 0.5 * views$known
  else views$phenotype                                   # agnostic: no known-gene seed view
}
percentile_in_universe <- function(scores, universe_scores) {       # declare the universe; record it in the manifest
  vapply(scores, function(s) if (is.na(s)) NA_real_ else mean(universe_scores <= s), numeric(1))
}

.mva_cache <- new.env()                                          # big resources are loaded once per R process
load_string_graph <- function(links_path, info_path, min_score = 700) {
  key <- paste0("string_", links_path, min_score); if (!is.null(.mva_cache[[key]])) return(.mva_cache[[key]])
  links <- read_gz(links_path, select = 1:3); data.table::setnames(links, c("protein1", "protein2", "combined_score"))
  links <- links[combined_score >= min_score]
  info <- read_gz(info_path, select = 1:2); data.table::setnames(info, c("protein", "gene_symbol"))
  info <- unique(info[!duplicated(protein)]); g1 <- info$gene_symbol[match(links$protein1, info$protein)]; g2 <- info$gene_symbol[match(links$protein2, info$protein)]
  e <- data.table(a = g1, b = g2, combined_score = links$combined_score)[!is.na(a) & !is.na(b)]
  logmsg("A5: STRING edges >= ", min_score, ": ", nrow(e), " (directed rows) | genes mapped: ", uniqueN(c(e$a, e$b)))
  out <- list(graph = build_string_graph(e, 0), edges = e); .mva_cache[[key]] <- out; out   # score filter already applied above
}
load_uniprot_gene_names <- function(path) {
  x <- data.table::fread(path, colClasses = "character", quote = ""); nm <- tolower(names(x))
  up <- names(x)[grep("uniprot|entry|accession|^from$", nm)[1]]; gn <- names(x)[grep("gene", nm)[1]]
  if (is.na(up) || is.na(gn)) stop("cannot detect UniProt/gene columns in ", path, "; found: ", paste(names(x), collapse = ", "))
  out <- x[, .(uniprot_id = get(up), gene_symbol = sub("[ ;].*$", "", get(gn)))][nzchar(gene_symbol) & !is.na(gene_symbol)]   # first symbol per entry
  unique(out)
}
load_reactome_graph <- function(mapping_path, names_path, max_pathway_size = 100) {
  r <- data.table::fread(mapping_path, header = FALSE, sep = "\t", quote = "", colClasses = "character",
                         col.names = c("uniprot_id", "pathway_id", "url", "pathway_name", "evidence", "species"))
  r <- r[species == "Homo sapiens"]; gm <- load_uniprot_gene_names(names_path)
  sets <- unique(merge(r[, .(uniprot_id, pathway_id)], gm, by = "uniprot_id")[, .(pathway_id, gene_symbol)])
  sets <- sets[, if (.N >= 2L && .N <= max_pathway_size) .SD, by = pathway_id]
  pairs <- sets[, { g <- sort(unique(gene_symbol)); cb <- utils::combn(g, 2); list(a = cb[1, ], b = cb[2, ]) }, by = pathway_id]
  e <- pairs[, .(weight = min(1, .N / 5)), by = .(a, b)]                 # more shared pathways -> stronger edge (declared, capped)
  logmsg("A5: Reactome pathway layer edges: ", nrow(e), " from ", uniqueN(sets$pathway_id), " pathways (size cap ", max_pathway_size, ")")
  list(graph = igraph::graph_from_data_frame(e[, .(from = a, to = b, weight)], directed = FALSE), sets = sets)
}

# Multiplex RWR: supra-adjacency over layers, inter-layer coupling lambda, dangling mass returns to the seed.
multiplex_rwr <- function(layers, seeds, restart, coupling = 0.5, iterations = 100L, tol = 1e-10) {
  nodes <- sort(unique(unlist(lapply(layers, function(g) igraph::V(g)$name)))); n <- length(nodes); L <- length(layers)
  seeds <- intersect(unique(seeds), nodes); if (!length(seeds)) return(NULL)
  Ps <- lapply(layers, function(g) { el <- igraph::as_edgelist(g); w <- igraph::E(g)$weight
    i <- match(c(el[, 1], el[, 2]), nodes); j <- match(c(el[, 2], el[, 1]), nodes)
    A <- Matrix::sparseMatrix(i, j, x = c(w, w), dims = c(n, n)); cs <- Matrix::colSums(A); list(P = A %*% Matrix::Diagonal(x = ifelse(cs > 0, 1 / cs, 0)), dangling = which(cs == 0)) })
  s0 <- numeric(n); s0[match(seeds, nodes)] <- 1; s <- rep(s0 / sum(s0), L) / L               # seed on every layer
  x <- s
  for (it in seq_len(iterations)) {
    blocks <- split(seq_len(n * L), rep(seq_len(L), each = n)); nxt <- numeric(n * L)
    for (l in seq_len(L)) { xl <- x[blocks[[l]]]
      intra <- (1 - coupling) * (as.numeric(Ps[[l]]$P %*% xl) + sum(xl[Ps[[l]]$dangling]) * (s0 / sum(s0)))
      nxt[blocks[[l]]] <- nxt[blocks[[l]]] + intra
      if (L > 1) for (m in setdiff(seq_len(L), l)) nxt[blocks[[m]]] <- nxt[blocks[[m]]] + coupling / (L - 1) * xl }
    nxt <- restart * s + (1 - restart) * nxt; done <- sum(abs(nxt - x)) < tol; x <- nxt; if (done) break }
  per_layer <- matrix(x, nrow = n); dimnames(per_layer) <- list(nodes, names(layers) %||% paste0("L", seq_len(L)))
  list(score = stats::setNames(rowSums(per_layer), nodes), per_layer = per_layer)
}
