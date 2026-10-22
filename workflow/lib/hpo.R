# workflow/lib/hpo.R -- A4: ontology, information content, Lin / Resnik / Wu-Palmer / ancestor Jaccard (pure R)
parse_obo <- function(path) {
  need_files(path, "HPO obo"); L <- readLines(path, warn = FALSE)
  st <- which(L == "[Term]"); en <- c(st[-1] - 1L, length(L))
  rows <- lapply(seq_along(st), function(i) { b <- L[(st[i] + 1L):en[i]]
    id <- sub("^id: ", "", b[startsWith(b, "id: ")][1]); isa <- sub(" .*$", "", sub("^is_a: ", "", b[startsWith(b, "is_a: ")]))
    list(id = id, parents = isa, obsolete = any(b == "is_obsolete: true"), alt = sub("^alt_id: ", "", b[startsWith(b, "alt_id: ")]),
         name = sub("^name: ", "", b[startsWith(b, "name: ")][1])) })
  names(rows) <- vapply(rows, `[[`, "", "id"); rows <- rows[!is.na(names(rows))]
  alt <- unlist(lapply(rows, function(r) stats::setNames(rep(r$id, length(r$alt)), r$alt)))
  if (is.null(alt)) alt <- stats::setNames(character(), character())
  list(terms = rows[!vapply(rows, `[[`, TRUE, "obsolete")], alt = alt)
}
build_ontology <- function(obo) {
  terms <- obo$terms; ids <- names(terms); anc <- new.env(hash = TRUE); dep <- new.env(hash = TRUE)
  get_anc <- function(id) { if (!is.null(anc[[id]])) return(anc[[id]])
    p <- intersect(terms[[id]]$parents, ids); a <- unique(c(id, unlist(lapply(p, get_anc), use.names = FALSE))); anc[[id]] <- a; a }
  get_dep <- function(id) { if (!is.null(dep[[id]])) return(dep[[id]])
    p <- intersect(terms[[id]]$parents, ids); d <- if (!length(p)) 1L else 1L + max(vapply(p, get_dep, 0L)); dep[[id]] <- d; d }
  for (i in ids) { get_anc(i); get_dep(i) }
  list(ids = ids, ancestors = as.list(anc), depth = unlist(as.list(dep)), alt = obo$alt, names = vapply(terms, function(t) t$name %||% NA_character_, ""))
}
clean_terms <- function(x, ont) { x <- unique(trimws(x[nzchar(x)])); a <- ont$alt[x]; x[!is.na(a)] <- unname(a[!is.na(a)]); x[x %in% ont$ids] }

# gene -> direct terms; propagated sets for IC; IC(t) = -log(n_t / N)
gene_term_sets <- function(g2p, ont) {
  need_cols(g2p, c("gene_symbol", "hpo_id"), "genes_to_phenotype")
  g2p <- unique(g2p[, .(gene_symbol, hpo_id)]); g2p[, hpo_id := clean_terms_vec(hpo_id, ont)]; g2p <- g2p[!is.na(hpo_id)]
  direct <- split(g2p$hpo_id, g2p$gene_symbol)
  prop <- lapply(direct, function(t) unique(unlist(ont$ancestors[t], use.names = FALSE)))
  N <- length(prop); cnt <- table(unlist(prop, use.names = FALSE))
  ic <- stats::setNames(-log(as.numeric(cnt) / N), names(cnt)); ic[setdiff(ont$ids, names(ic))] <- max(ic) # unseen terms: maximal IC
  list(direct = direct, propagated = prop, ic = ic)
}
clean_terms_vec <- function(x, ont) { a <- ont$alt[x]; x[!is.na(a)] <- unname(a[!is.na(a)]); ifelse(x %in% ont$ids, x, NA_character_) }

term_pair_sims <- function(a, B, ont, ic) {           # one patient term vs a vector of terms
  anc_a <- ont$ancestors[[a]]
  r <- vapply(B, function(b) { com <- intersect(anc_a, ont$ancestors[[b]]); if (!length(com)) return(c(0, 0, 0))
    m <- com[which.max(ic[com])]; ra <- ic[m]; den <- ic[a] + ic[b]
    c(if (den > 0) 2 * ra / den else as.numeric(a == b), ra, 2 * ont$depth[[com[which.max(ont$depth[com])]]] / (ont$depth[[a]] + ont$depth[[b]])) }, numeric(3))
  matrix(r, nrow = 3, dimnames = list(c("lin", "resnik", "wup"), B))
}
bma <- function(S) if (!length(S) || nrow(S) == 0L || ncol(S) == 0L) NA_real_ else mean(c(mean(apply(S, 1, max)), mean(apply(S, 2, max))))

# A4 -- gene-level HPO evidence for the shortlist genes
compute_hpo_evidence <- function(genes, patient_terms, ont, gts) {
  pt <- clean_terms(patient_terms, ont); logmsg("A4: patient HPO terms provided ", length(patient_terms), ", active in ontology ", length(pt))
  if (!length(pt)) stop("no active patient HPO terms")
  genes <- unique(genes); have <- genes[genes %in% names(gts$direct)]
  U <- unique(unlist(gts$direct[have], use.names = FALSE)); sims <- lapply(pt, function(a) term_pair_sims(a, U, ont, gts$ic))
  anc_p <- unique(unlist(ont$ancestors[pt], use.names = FALSE))
  res <- rbindlist(lapply(genes, function(g) {
    d <- gts$direct[[g]]
    if (is.null(d)) return(data.table(gene_symbol = g, hpo_lin_bma = NA_real_, hpo_resnik_bma = NA_real_, hpo_wup_bma = NA_real_, hpo_anc_jaccard = NA_real_,
                                      n_gene_terms = 0L, direct_overlap = 0L, annotation_status = "no_annotation"))
    pick <- function(m) bma(do.call(rbind, lapply(sims, function(s) s[m, d, drop = FALSE])))
    ag <- gts$propagated[[g]]
    data.table(gene_symbol = g, hpo_lin_bma = pick("lin"), hpo_resnik_bma = pick("resnik"), hpo_wup_bma = pick("wup"),
               hpo_anc_jaccard = length(intersect(anc_p, ag)) / length(union(anc_p, ag)), n_gene_terms = length(d),
               direct_overlap = length(intersect(pt, ag)), annotation_status = if (length(d) >= 10) "adequate_10plus_terms" else "sparse_annotation")
  }))
  for (m in c("hpo_lin_bma", "hpo_resnik_bma", "hpo_wup_bma", "hpo_anc_jaccard")) {
    v <- res[[m]]; ok <- is.finite(v)
    set(res, j = paste0(m, "_rank"), value = ifelse(ok, frank(-v, na.last = "keep", ties.method = "min"), NA_integer_))
    set(res, j = paste0(m, "_pct"), value = ifelse(ok, vapply(v, function(s) if (is.na(s)) NA_real_ else mean(v[ok] <= s), 0), NA_real_)) }
  res[]
}

# proband-only HPO terms (family-history terms are never patient terms); validates syntax
read_proband_terms <- function(proband_path, family_path = NULL) {
  x <- trimws(readLines(need_files(proband_path, "proband HPO"), warn = FALSE)); x <- x[nzchar(x) & !startsWith(x, "#")]
  bad <- x[!grepl("^HP:[0-9]{7}$", x)]; if (length(bad)) stop("invalid HPO ids: ", paste(bad, collapse = ", "), call. = FALSE)
  if (!is.null(family_path) && file.exists(family_path)) { f <- trimws(readLines(family_path, warn = FALSE)); ov <- intersect(x, f[nzchar(f)])
    if (length(ov)) logmsg("note: ", length(ov), " proband term(s) also in the family-history file; proband file is authoritative") }
  unique(x)
}
