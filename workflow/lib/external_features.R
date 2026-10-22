# workflow/lib/external_features.R -- Layer B: splits, feature definition shared by B9 and C1, conditional softmax, stacking, conformal
suppressPackageStartupMessages(library(data.table))

# B10 -- whole duplicate groups are assigned together; runs BEFORE feature fitting
define_casewise_splits <- function(cand, seed, frac = c(train = 0.70, validation = 0.15, locked_test = 0.15)) {
  stopifnot(all(c("case_id", "duplicate_group") %in% names(cand)), abs(sum(frac) - 1) < 1e-9)
  set.seed(seed)
  groups <- sample(unique(cand$duplicate_group)); n <- length(groups)
  if (n < 3L) stop("need >= 3 duplicate groups for three splits")
  n_tr <- max(1L, floor(frac[["train"]] * n)); n_va <- max(1L, floor(frac[["validation"]] * n))
  if (n_tr + n_va >= n) n_tr <- n - n_va - 1L
  grp <- data.table(duplicate_group = groups, split = rep(c("train", "validation", "locked_test"), c(n_tr, n_va, n - n_tr - n_va)))
  man <- unique(merge(unique(cand[, .(case_id, duplicate_group)]), grp, by = "duplicate_group"))
  for (a in c("case_id", "duplicate_group")) {            # disjointness is a hard check, not a hope
    chk <- man[, uniqueN(split), by = c(a)]; stopifnot(all(chk$V1 == 1L))
  }
  man[]
}

index_resources <- function(hpo_by_gene, hpo_genes_by_term, string_edges) {
  both <- rbind(string_edges[, .(g = a, nb = b, score)], string_edges[, .(g = b, nb = a, score)])
  list(hpo_by_gene = hpo_by_gene, hpo_genes_by_term = hpo_genes_by_term, string_by_gene = split(both[, .(nb, score)], both$g))
}

# raw G, H, N (+ availability flags). Leave-one-case-out gene prior: a case never counts itself.
raw_features <- function(cand, case_terms, res, state) {
  n <- nrow(cand); G <- H <- N <- M_H <- M_N <- numeric(n)
  for (i in seq_len(n)) {
    g <- cand$gene[i]; cid <- cand$case_id[i]; terms <- unique(case_terms[[cid]])
    ann <- res$hpo_by_gene[[g]]
    M_H[i] <- as.numeric(!is.null(ann)); H[i] <- if (!is.null(ann) && length(terms)) length(intersect(terms, ann)) / length(terms) else 0
    nb <- res$string_by_gene[[g]]; M_N[i] <- as.numeric(!is.null(nb))
    seeds <- unique(unlist(res$hpo_genes_by_term[terms], use.names = FALSE))
    N[i] <- if (!is.null(nb) && length(seeds)) max(c(0, nb$score[nb$nb %in% seeds])) / 1000 else 0
    cnt <- state$gene_count[g]; cnt <- if (is.na(cnt)) 0 else cnt
    own <- as.integer(identical(unname(state$case_gene[cid]), g)); G[i] <- log1p(max(cnt - own, 0))
  }
  data.table(G = G, H = H, N = N, M_H = M_H, M_N = M_N)
}

scale_features <- function(raw, state) {
  out <- data.table::copy(raw)
  for (f in c("G", "H", "N")) set(out, j = f, value = pmin(pmax((raw[[f]] - state$min[[f]]) / (state$max[[f]] - state$min[[f]]), 0), 1))
  out[, GxH := G * H][]                                  # interactions AFTER scaling
}

fit_feature_state <- function(train, case_terms, res) {  # the ONLY place a state is estimated from data
  pos <- train[label == 1L]
  st <- list(gene_count = table(pos$gene) |> unclass() |> (\(x) stats::setNames(as.numeric(x), names(x)))(),
             case_gene = stats::setNames(pos$gene, pos$case_id), schema = c("G", "H", "N", "GxH"))
  raw <- raw_features(train, case_terms, res, st)
  st$min <- vapply(raw[, .(G, H, N)], min, numeric(1)); st$max <- vapply(raw[, .(G, H, N)], max, numeric(1))
  st$max[st$max <= st$min] <- st$min[st$max <= st$min] + 1
  st
}

# B9 and C1 call THIS function; C1 passes the frozen state and private candidates
compute_features <- function(cand, case_terms, res, state) {
  cbind(cand[, .(case_id, candidate_id, gene)], scale_features(raw_features(cand, case_terms, res, state), state))
}

predict_softmax <- function(X, beta, case_id) {
  eta <- as.numeric(X %*% beta); p <- exp(eta - ave(eta, case_id, FUN = max)); p / ave(p, case_id, FUN = sum)
}

# B11 -- conditional softmax (no intercept), ridge, analytic gradient
fit_conditional_softmax <- function(X, case_id, label, lambda = 1e-6) {
  X <- as.matrix(X); idx <- split(seq_along(case_id), case_id)
  pos <- vapply(idx, function(i) { w <- which(label[i] == 1L); stopifnot(length(w) == 1L); i[w] }, integer(1))   # exactly one positive per case
  nll <- function(b) { eta <- as.numeric(X %*% b)
    -sum(vapply(names(idx), function(k) { e <- eta[idx[[k]]]; m <- max(e); eta[pos[[k]]] - m - log(sum(exp(e - m))) }, 0)) + lambda * sum(b^2) }
  grd <- function(b) { eta <- as.numeric(X %*% b); g <- 2 * lambda * b
    for (k in names(idx)) { i <- idx[[k]]; e <- eta[i]; p <- exp(e - max(e)); p <- p / sum(p); g <- g + colSums(X[i, , drop = FALSE] * p) - X[pos[[k]], ] }
    g }
  fit <- optim(rep(0, ncol(X)), nll, grd, method = "BFGS", hessian = TRUE, control = list(maxit = 500, reltol = 1e-12))
  list(beta = stats::setNames(fit$par, colnames(X)), nll = fit$value, converged = fit$convergence == 0, vcov = tryCatch(solve(fit$hessian), error = function(e) NULL))
}

case_metrics <- function(p, case_id, label, ks = c(1, 5)) {
  d <- data.table(case_id, p, label)[, rk := frank(-p, ties.method = "min"), by = case_id][label == 1L]
  c(list(mrr = mean(1 / d$rk), nll = -mean(log(pmax(d$p, 1e-12)))), stats::setNames(lapply(ks, function(k) mean(d$rk <= k)), paste0("top", ks)))
}

stack_weights <- function(prob_causal) {                 # matrix: cases x models, causal-candidate probability
  lpd <- log(pmax(as.matrix(prob_causal), 1e-12))
  w <- as.numeric(loo::stacking_weights(lpd)); stats::setNames(w, colnames(lpd))
}
apply_stack <- function(pred_long, w) {                  # pred_long: case_id, candidate_id, model_id, probability
  d <- data.table::as.data.table(pred_long); d[, wp := probability * unname(w[model_id])]
  ans <- d[, .(probability = sum(wp)), by = .(case_id, candidate_id)]
  stopifnot(all(abs(ans[, sum(probability), by = case_id]$V1 - 1) < 1e-8)); ans[]
}
conformal_qhat <- function(nonconf, alpha = 0.1) {
  n <- length(nonconf); k <- ceiling((n + 1) * (1 - alpha)); if (k > n) return(Inf)   # too few calibration cases
  sort(nonconf)[k]
}

# ---- external cohort helpers (B1-B8) ----
load_external_resources <- function(hpo_g2p_path, string_edges) {      # string_edges: data.table(a, b, combined_score)
  g2p <- data.table::fread(need_files(hpo_g2p_path, "gene_to_hpo"), sep = "\t", quote = "")
  need_cols(g2p, c("gene_symbol", "hpo_id"), "genes_to_phenotype"); g2p <- unique(g2p[, .(gene_symbol, hpo_id)])
  e <- unique(string_edges[, .(a, b, score = combined_score)])
  list(res = index_resources(split(g2p$hpo_id, g2p$gene_symbol), split(g2p$gene_symbol, g2p$hpo_id), e), universe = intersect(unique(g2p$gene_symbol), unique(c(e$a, e$b))))
}
read_external_cases <- function(path, min_observed_hpo = 3) {
  x <- read_gz(path, colClasses = "character")
  need_cols(x, c("case_id", "duplicate_group_key", "causal_gene_symbol_normalized", "observed_hpo_ids_normalized", "n_observed_hpo_terms", "has_hpo_and_string_mapping"), "deduplicated cases")
  n0 <- nrow(x); x[, n_observed_hpo_terms := as.integer(n_observed_hpo_terms)]
  x <- x[n_observed_hpo_terms >= min_observed_hpo & has_hpo_and_string_mapping %in% c("1", "TRUE", "true") & nzchar(causal_gene_symbol_normalized)]
  logmsg("B1-B7: retained ", nrow(x), " / ", n0, " deduplicated cases")
  if (anyDuplicated(x$case_id)) stop("duplicate case_id in external cases"); assert_no_private(x)
  x[, `:=`(duplicate_group = duplicate_group_key, true_gene = causal_gene_symbol_normalized)]; x[]
}
build_candidate_sets <- function(cases, universe, n_neg, seed) {        # background = HPO+STRING-mapped genes minus the positive
  set.seed(seed)
  rbindlist(lapply(seq_len(nrow(cases)), function(i) { tg <- cases$true_gene[i]; neg <- sample(setdiff(universe, tg), min(n_neg, length(universe) - 1L)); genes <- c(tg, neg)
    data.table(case_id = cases$case_id[i], duplicate_group = cases$duplicate_group[i], gene = genes, candidate_id = paste0(cases$case_id[i], "_", seq_along(genes)), label = c(1L, rep(0L, length(neg)))) }))
}
assert_no_private <- function(x) if (any(grepl("data/private|results/private|WGS_EX|MVA_WGS/data", unlist(x), ignore.case = TRUE))) stop("private reference in a Layer B artifact", call. = FALSE)
