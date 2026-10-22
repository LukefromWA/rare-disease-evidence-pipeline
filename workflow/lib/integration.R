# workflow/lib/integration.R -- A7 coverage-aware integration (informed/agnostic), A7b dependence, A7c Shapley, A8 mosaicism, A10 aggregation
suppressPackageStartupMessages(library(data.table))

# Weights matched BY NAME. Missing evidence = missing coverage, never zero.
weighted_available_mean <- function(scores, weights) {
  scores <- as.matrix(scores)
  stopifnot(!is.null(colnames(scores)), all(colnames(scores) %in% names(weights)))
  w <- weights[colnames(scores)]
  avail <- is.finite(scores)
  num <- rowSums(ifelse(avail, sweep(scores, 2L, w, "*"), 0))
  den <- rowSums(sweep(avail * 1, 2L, w, "*"))
  out <- num / den
  out[den == 0] <- NA_real_
  attr(out, "n_views") <- rowSums(avail)
  out
}


# channels: data.table(candidate_id, gene, genomic, hpo, string, multiplex)   (each in [0,1] or NA)
integrate_tier1 <- function(channels, weights, known_genes, alpha_boost = 2) {
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
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 5L) return(list(rho = NA_real_, spearman = NA_real_, n = n))
  zx <- qnorm(rank(x[ok]) / (n + 1)); zy <- qnorm(rank(y[ok]) / (n + 1))
  list(rho = cor(zx, zy), spearman = cor(x[ok], y[ok], method = "spearman"), n = n)
}

# Exact Shapley over channels for ONE candidate (A7c). v(S) = coverage-aware mean over S; v(empty) = 0.
shapley_attribution <- function(score_row, weights) {
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

# A10 -- gene level (best variant per gene), coverage and stability classes, review tier. Same column names as the recorded A10.
aggregate_ranked_candidates <- function(sl, freq) {
  d <- merge(sl, freq[, .(candidate_id, p_top_m = freq)], by = "candidate_id", all.x = TRUE, sort = FALSE)
  d[, `:=`(has_genomic = is.finite(genomic), has_hpo = is.finite(hpo_lin_bma_pct), has_string = is.finite(string_informed_pct), has_multiplex = is.finite(multiplex_informed_pct))]
  d[, coverage_class := fcase(has_genomic & has_hpo & has_string & has_multiplex & genomic >= .7 & hpo_lin_bma_pct >= .7 & pmin(string_informed_pct, multiplex_informed_pct) >= .7, "three_view_concordant",
                              has_genomic & has_hpo & has_string & has_multiplex, "three_view_mixed",
                              has_genomic & has_hpo & (has_string | has_multiplex), "two_view_plus_partial_network", default = "other_or_unscorable")]
  d[, rank_stability_class := fcase(p_top_m >= .90, "topm_robust", p_top_m >= .5, "topm_moderate", default = "scenario_sensitive")]
  d[, review_tier := fcase(coverage_class == "three_view_concordant" & rank_stability_class == "topm_robust", "A_robust_lead",
                           coverage_class %in% c("three_view_concordant", "three_view_mixed"), "B_supported", default = "C_extended_review")]
  data.table::setorder(d, -posterior_informed, gene_symbol)
  g <- d[!is.na(gene_symbol), .SD[1L], by = gene_symbol]; g[, n_shortlist_variants := d[!is.na(gene_symbol), .N, by = gene_symbol][match(g$gene_symbol, gene_symbol), N]]
  g[, gene_priority_rank := frank(-posterior_informed, ties.method = "first")]
  out <- g[, .(gene_priority_rank, gene_symbol, best_candidate_id = candidate_id, n_shortlist_variants, best_integrated_score = S_bar_informed,
               best_HPO_semantic_score = hpo_lin_bma_pct, best_PPI_network_score = string_informed_pct, best_multiplex_score = multiplex_informed_pct,
               best_GQ = GQ, best_MAX_AF = MAX_AF, best_VAF = VAF, best_DP = DP, best_AD_alt = AD_alt, best_IMPACT = IMPACT, best_consequence = consequence,
               posterior_informed, posterior_agnostic, rank_informed, rank_agnostic, n_views, coverage_class, rank_stability_class, review_tier,
               best_low_vaf_review = mosaic_review)]
  data.table::setorder(out, gene_priority_rank); out[]
}
