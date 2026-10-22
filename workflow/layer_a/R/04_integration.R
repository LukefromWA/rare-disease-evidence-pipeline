# 04_integration.R -- A7 integration (informed/agnostic), A7b dependence, A7c Shapley, A8 mosaicism
suppressPackageStartupMessages(library(data.table))

# channels: data.table(candidate_id, gene, genomic, hpo, string, multiplex)   (each in [0,1] or NA)
integrate_tier1 <- function(channels, weights, known_genes, alpha_boost = 2) {
    logmsg("04: tier1 integration started")
  assert_columns(channels, c("candidate_id", "gene", names(weights)), "channels")
  ch <- data.table::copy(data.table::as.data.table(channels))
  cols <- names(weights)
  sc <- as.matrix(ch[, ..cols]); sc[!is.finite(sc)] <- NA
  S  <- weighted_available_mean(sc, weights)
  ch[, `:=`(S_bar = as.numeric(S), n_views = attr(S, "n_views"))]
  ch[, scorable := is.finite(S_bar) & n_views > 0]
  prior_informed <- data.table::fifelse(ch$gene %in% known_genes, alpha_boost, 1)
  post <- function(prior) { l <- prior * data.table::fifelse(ch$scorable, ch$S_bar, 0); if (sum(l) == 0) l else l / sum(l) }
  ch[, `:=`(posterior_informed = post(prior_informed), posterior_agnostic = post(rep(1, .N)))]
  ch[, `:=`(rank_informed = frank(-posterior_informed, ties.method = "min"),
            rank_agnostic = frank(-posterior_agnostic, ties.method = "min"))]
  ch[]
}

# Gaussian-copula rho on complete cases (A7b)
estimate_dependence <- function(x, y) {
    logmsg("04: dependence estimation started")
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 5L) return(list(rho = NA_real_, spearman = NA_real_, n = n))
  zx <- qnorm(rank(x[ok]) / (n + 1)); zy <- qnorm(rank(y[ok]) / (n + 1))
  list(rho = cor(zx, zy), spearman = cor(x[ok], y[ok], method = "spearman"), n = n)
}

# Exact Shapley over channels for ONE candidate (A7c). v(S) = coverage-aware mean over S; v(empty) = 0.
shapley_attribution <- function(score_row, weights) {
    logmsg("04: Shapley attribution started")
  ch <- names(weights); K <- length(ch)
  v <- function(S) {
    if (!length(S)) return(0)
    z <- weighted_available_mean(matrix(score_row[S], 1, dimnames = list(NULL, S)), weights)[1]
    if (is.na(z)) 0 else z                               # NA subset (nothing observed) contributes 0
  }
  phi <- stats::setNames(numeric(K), ch)
  for (k in ch) {
    others <- setdiff(ch, k)
    for (m in 0:length(others)) {
      subs <- if (m == 0) list(character()) else combn(others, m, simplify = FALSE)
      wgt <- factorial(m) * factorial(K - m - 1) / factorial(K)
      for (S in subs) phi[k] <- phi[k] + wgt * (v(c(S, k)) - v(S))
    }
  }
  attr(phi, "full_minus_empty") <- v(ch) - v(character())      # efficiency check: sum(phi) equals this
  phi
}

# A8 -- beta-binomial posterior on informative depth, and germline-het vs mosaic log Bayes factor
run_tier9 <- function(tbl, qc, alpha0 = 1, beta0 = 1) {
  t <- data.table::copy(data.table::as.data.table(tbl)); n <- t$AD_ref + t$AD_alt
  t[, included_tier9 := DP > qc$dp_gt & GQ >= qc$min_gq & n > 0]
  t[, VAF_observed := ifelse(n > 0, AD_alt / n, NA_real_)]
  t[, `:=`(post_alpha = alpha0 + AD_alt, post_beta = beta0 + AD_ref)]
  t[, VAF_post_mean := post_alpha / (post_alpha + post_beta)]
  t[, `:=`(VAF_95CrI_lower = qbeta(0.025, post_alpha, post_beta), VAF_95CrI_upper = qbeta(0.975, post_alpha, post_beta))]
  # log BF > 0 favors a free (mosaic) allele fraction; < 0 favors a germline het (p = 0.5)
  t[, log_bf_mosaic_vs_het := lchoose(n, AD_alt) + lbeta(AD_alt + alpha0, n - AD_alt + beta0) - lbeta(alpha0, beta0) -
                              dbinom(AD_alt, n, 0.5, log = TRUE)]
  t[, mosaic_review := included_tier9 & AD_alt >= qc$mosaic_min_alt_reads & VAF_observed >= qc$mosaic_vaf_range[1] &
                       VAF_observed < qc$mosaic_vaf_range[2]]
  t[]
}
