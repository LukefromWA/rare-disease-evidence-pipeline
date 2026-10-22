# workflow/lib/stability.R -- A9: perturb, rerank, top-m frequency, Kuncheva, leave-one-channel-out
suppressPackageStartupMessages(library(data.table))

rdirichlet_base <- function(n, alpha) {
  d <- matrix(rgamma(n * length(alpha), shape = pmax(alpha, 1e-6)), nrow = n, byrow = TRUE)
  d / rowSums(d)
}

kuncheva_index <- function(a, b, n) {                   # two top-k id sets from n candidates
  k <- length(a); stopifnot(length(b) == k, k > 0, k < n)
  r <- length(intersect(a, b)); (r * n - k^2) / (k * (n - k))
}

# score matrix (rows = candidates, cols = channels); weights perturbed on the simplex
rerank <- function(sc, w) { s <- weighted_available_mean(sc, w); s[is.na(s)] <- -Inf; rank(-s, ties.method = "min") }

run_stability <- function(sc, ids, weights, B = 1000L, top_m = 5L, kappa = 100, seed) {
  set.seed(seed)                                         # seed once; never inside the loop
  W <- rdirichlet_base(B, kappa * weights / sum(weights)); colnames(W) <- names(weights)
  base_top <- ids[order(rerank(sc, weights))][seq_len(top_m)]
  tops <- vector("list", B)
  for (b in seq_len(B)) tops[[b]] <- ids[order(rerank(sc, W[b, ]))][seq_len(top_m)]
  freq <- table(factor(unlist(tops), levels = ids)) / B
  ki   <- vapply(tops, kuncheva_index, numeric(1), b = base_top, n = length(ids))
  list(top_m_frequency = data.table(candidate_id = ids, freq = as.numeric(freq[ids])),
       kuncheva_mean = mean(ki), base_top = base_top)
}

leave_one_channel_out <- function(sc, weights) {
  full <- weighted_available_mean(sc, weights)
  rbindlist(lapply(colnames(sc), function(ch) {
    keep <- setdiff(colnames(sc), ch)
    s <- weighted_available_mean(sc[, keep, drop = FALSE], weights[keep])
    data.table(channel_removed = ch, candidate = seq_len(nrow(sc)), score_delta = as.numeric(s) - as.numeric(full),
               rank_shift = rank(-ifelse(is.na(s), -Inf, s), ties.method = "min") - rank(-ifelse(is.na(full), -Inf, full), ties.method = "min"))
  }))
}
