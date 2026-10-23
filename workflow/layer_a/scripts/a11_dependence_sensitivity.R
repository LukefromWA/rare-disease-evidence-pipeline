#!/usr/bin/env Rscript

# Dirichlet/Beta Monte Carlo sensitivity analysis for N-of-1 candidate ranking.
# Quantifies rank robustness under score, weight, network-view, and dependence
# assumptions. It is not a calibrated pathogenicity or causal probability model.

options(stringsAsFactors = FALSE)

.source_file <- normalizePath(sub("^--file=", "", grep("^--file=",commandArgs(FALSE),value=TRUE)[1]))
.runtime_file <- normalizePath(file.path(dirname(.source_file),"..","..","lib","a11_runtime.R"))
source(.runtime_file)
.cli <- a11_cli(c("input","outdir"),list(seed="20260912",draws="10000",`focus-gene`=""))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Required preinstalled package: data.table")
}
suppressPackageStartupMessages(library(data.table))

INPUT <- .cli[["input"]]
OUT_DIR <- .cli[["outdir"]]
FOCUS_GENE <- .cli[["focus-gene"]]
N_DRAWS <- a11_integer(.cli[["draws"]],"draws",1L)
SEED <- a11_integer(.cli[["seed"]],"seed",0L)
BASE_WEIGHTS <- c(genomic = 0.45, hpo = 0.35, network = 0.20)
NETWORK_MODELS <- c("string", "multiplex")
GAMMA_VALUES <- c(0, 0.5, 1)
SCENARIOS <- list(
  conservative = list(concentration = 20, precision = c(genomic = 75, hpo = 40, network = 40)),
  moderate = list(concentration = 80, precision = c(genomic = 250, hpo = 100, network = 100)),
  optimistic = list(concentration = 200, precision = c(genomic = 500, hpo = 250, network = 250))
)
CHANNELS <- c("genomic", "hpo", "network")

.ctx <- a11_begin(OUT_DIR,c(evidence=INPUT),.source_file,.runtime_file)
set.seed(SEED)
.rng_initial <- .Random.seed

clamp01 <- function(x, eps = 1e-6) pmin(pmax(x, eps), 1 - eps)

rdirichlet_one <- function(alpha) {
  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  names(x) <- names(alpha)
  x / sum(x)
}

draw_beta_score <- function(x, precision) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (any(ok)) {
    mu <- clamp01(x[ok])
    out[ok] <- rbeta(sum(ok), mu * precision, (1 - mu) * precision)
  }
  out
}

rank_deterministic <- function(score, candidate_id) {
  out <- rep(NA_real_, length(score))
  ok <- is.finite(score)
  if (!any(ok)) return(out)
  ord <- order(-score[ok], candidate_id[ok])
  out[which(ok)[ord]] <- seq_along(ord)
  out
}

weighted_score <- function(G, H, N, weights) {
  X <- cbind(genomic = G, hpo = H, network = N)
  observed <- is.finite(X)
  W <- matrix(weights, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
  numerator <- rowSums(X * W, na.rm = TRUE)
  denominator <- rowSums(W * observed)
  out <- numerator / denominator
  out[denominator == 0] <- NA_real_
  out
}

gaussian_copula_cor <- function(X) {
  X <- as.matrix(X)
  colnames(X) <- CHANNELS
  keep <- complete.cases(X)
  n_complete <- sum(keep)

  if (n_complete < 20L) {
    R <- diag(length(CHANNELS))
    dimnames(R) <- list(CHANNELS, CHANNELS)
    return(list(R = R, n_complete = n_complete))
  }

  U <- apply(X[keep, , drop = FALSE], 2, function(x) {
    rank(x, ties.method = "average") / (length(x) + 1)
  })

  U <- as.matrix(U)
  colnames(U) <- CHANNELS
  R <- cor(qnorm(clamp01(U)), use = "pairwise.complete.obs")
  R <- as.matrix(R)
  R[!is.finite(R)] <- 0
  diag(R) <- 1
  dimnames(R) <- list(CHANNELS, CHANNELS)
  list(R = R, n_complete = n_complete)
}

shrink_correlation <- function(R, shrinkage = 0.50) {
  R <- as.matrix(R)
  dimnames(R) <- list(CHANNELS, CHANNELS)
  R2 <- (1 - shrinkage) * R + shrinkage * diag(length(CHANNELS))
  eig <- eigen(R2, symmetric = TRUE)
  eig$values[eig$values < 1e-6] <- 1e-6
  R2 <- cov2cor(eig$vectors %*% diag(eig$values) %*% t(eig$vectors))
  dimnames(R2) <- list(CHANNELS, CHANNELS)
  R2
}

adjust_weights <- function(w, rho_hn, gamma) {
  if (!identical(names(w), names(BASE_WEIGHTS))) {
    stop("Dirichlet weight names are invalid.")
  }
  w["network"] <- w["network"] / (1 + gamma * abs(rho_hn))
  w / sum(w)
}

summarize_draws <- function(score_matrix, candidate_id, gene_symbol) {
  rank_matrix <- apply(score_matrix, 2, rank_deterministic, candidate_id = candidate_id)
  rank_matrix <- as.matrix(rank_matrix)

  data.table(
    candidate_id = candidate_id,
    gene_symbol = gene_symbol,
    posterior_score_mean = rowMeans(score_matrix, na.rm = TRUE),
    posterior_score_q025 = apply(score_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
    posterior_score_q975 = apply(score_matrix, 1, quantile, probs = 0.975, na.rm = TRUE),
    posterior_rank_mean = rowMeans(rank_matrix, na.rm = TRUE),
    posterior_rank_median = apply(rank_matrix, 1, median, na.rm = TRUE),
    posterior_rank_q025 = apply(rank_matrix, 1, quantile, probs = 0.025, na.rm = TRUE),
    posterior_rank_q975 = apply(rank_matrix, 1, quantile, probs = 0.975, na.rm = TRUE),
    p_rank1 = rowMeans(rank_matrix == 1, na.rm = TRUE),
    p_top3 = rowMeans(rank_matrix <= 3, na.rm = TRUE),
    p_top5 = rowMeans(rank_matrix <= 5, na.rm = TRUE),
    p_top10 = rowMeans(rank_matrix <= 10, na.rm = TRUE)
  )
}

if (!file.exists(INPUT)) stop("Input file not found: ", INPUT)
dat <- fread(INPUT)
required <- c("candidate_id", "gene_symbol", "genomic", "hpo", "string", "multiplex")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

for (nm in c("genomic", "hpo", "string", "multiplex")) {
  raw <- dat[[nm]]
  numeric_value <- suppressWarnings(as.numeric(raw))
  if (any(!is.na(raw) & is.na(numeric_value)) ||
      any(is.nan(numeric_value) | is.infinite(numeric_value))) {
    stop(nm, " contains invalid numeric or nonfinite values.")
  }
  dat[, (nm) := numeric_value]
  bad <- is.finite(dat[[nm]]) & (dat[[nm]] < 0 | dat[[nm]] > 1)
  if (any(bad)) stop(nm, " contains values outside [0, 1].")
}

if (anyNA(dat$candidate_id) || any(!nzchar(trimws(as.character(dat$candidate_id)))) ||
    anyNA(dat$gene_symbol) || any(!nzchar(trimws(as.character(dat$gene_symbol)))) || anyDuplicated(dat$candidate_id)) stop("Invalid candidate/gene keys")
.coverage <- dat[,.(candidate_id,gene_symbol,entered_analysis=rowSums(cbind(is.finite(genomic),is.finite(hpo),is.finite(string)|is.finite(multiplex)))>=2)]
fwrite(.coverage,file.path(OUT_DIR,"candidate_coverage.csv"))
dat <- dat[rowSums(cbind(is.finite(genomic), is.finite(hpo), is.finite(string) | is.finite(multiplex))) >= 2]
if (nrow(dat)<2L) stop("At least two eligible candidates required by this legacy implementation")
if (!nrow(dat)) stop("No candidates have at least two usable evidence views.")
if (anyDuplicated(dat$candidate_id)) stop("candidate_id must be unique.")

copula_string <- gaussian_copula_cor(cbind(genomic = dat$genomic, hpo = dat$hpo, network = dat$string))
copula_multiplex <- gaussian_copula_cor(cbind(genomic = dat$genomic, hpo = dat$hpo, network = dat$multiplex))
R_string <- shrink_correlation(copula_string$R)
R_multiplex <- shrink_correlation(copula_multiplex$R)
rho_string <- R_string["hpo", "network"]
rho_multiplex <- R_multiplex["hpo", "network"]

fwrite(as.data.table(R_string, keep.rownames = "channel"), file.path(OUT_DIR, "string_shrunk_copula_R.csv"))
fwrite(as.data.table(R_multiplex, keep.rownames = "channel"), file.path(OUT_DIR, "multiplex_shrunk_copula_R.csv"))
fwrite(data.table(
  network_model = c("string", "multiplex"),
  hpo_network_gaussian_copula_rho = c(rho_string, rho_multiplex),
  absolute_rho = abs(c(rho_string, rho_multiplex)),
  n_complete = c(copula_string$n_complete, copula_multiplex$n_complete)
), file.path(OUT_DIR, "copula_dependence_summary.csv"))

message("Candidates entering analysis: ", nrow(dat))
message(sprintf("HPO-STRING rho: %.4f", rho_string))
message(sprintf("HPO-multiplex rho: %.4f", rho_multiplex))

all_results <- list()
weight_summary_rows <- list()
model_no <- 0L
total_models <- length(SCENARIOS) * length(NETWORK_MODELS) * length(GAMMA_VALUES)

for (scenario_name in names(SCENARIOS)) {
  scenario <- SCENARIOS[[scenario_name]]
  alpha_w <- BASE_WEIGHTS * scenario$concentration

  for (network_model in NETWORK_MODELS) {
    N_base <- dat[[network_model]]
    rho_hn <- if (network_model == "string") rho_string else rho_multiplex

    for (gamma in GAMMA_VALUES) {
      model_no <- model_no + 1L
      key <- paste(scenario_name, network_model, sprintf("gamma%.2f", gamma), sep = "__")
      message(sprintf("[%d/%d] %s", model_no, total_models, key))

      score_matrix <- matrix(NA_real_, nrow(dat), N_DRAWS)
      weight_matrix <- matrix(NA_real_, N_DRAWS, length(CHANNELS), dimnames = list(NULL, CHANNELS))

      for (d in seq_len(N_DRAWS)) {
        w_raw <- rdirichlet_one(alpha_w)
        names(w_raw) <- names(BASE_WEIGHTS)
        w <- adjust_weights(w_raw, rho_hn, gamma)

        if (length(w) != 3L || !identical(names(w), names(BASE_WEIGHTS))) {
          stop("Invalid adjusted weight vector in model: ", key)
        }

        weight_matrix[d, ] <- w
        score_matrix[, d] <- weighted_score(
          draw_beta_score(dat$genomic, scenario$precision["genomic"]),
          draw_beta_score(dat$hpo, scenario$precision["hpo"]),
          draw_beta_score(N_base, scenario$precision["network"]),
          w
        )
      }

      result <- summarize_draws(score_matrix, dat$candidate_id, dat$gene_symbol)
      result[, `:=`(scenario = scenario_name, network_model = network_model, gamma = gamma)]
      all_results[[key]] <- result
      fwrite(result, file.path(OUT_DIR, paste0("posterior_", key, ".csv")))

      W <- as.data.table(weight_matrix)
      weight_summary_rows[[key]] <- data.table(
        model = key,
        genomic_mean = mean(W$genomic), genomic_q025 = quantile(W$genomic, 0.025), genomic_q975 = quantile(W$genomic, 0.975),
        hpo_mean = mean(W$hpo), hpo_q025 = quantile(W$hpo, 0.025), hpo_q975 = quantile(W$hpo, 0.975),
        network_mean = mean(W$network), network_q025 = quantile(W$network, 0.025), network_q975 = quantile(W$network, 0.975)
      )
    }
  }
}

combined <- rbindlist(all_results, idcol = "model", fill = TRUE)
consensus <- combined[
  ,
  .(
    posterior_score_mean = mean(posterior_score_mean, na.rm = TRUE),
    posterior_rank_mean = mean(posterior_rank_mean, na.rm = TRUE),
    p_rank1 = mean(p_rank1, na.rm = TRUE),
    p_top3 = mean(p_top3, na.rm = TRUE),
    p_top5 = mean(p_top5, na.rm = TRUE),
    p_top10 = mean(p_top10, na.rm = TRUE),
    score_lower_envelope = min(posterior_score_q025, na.rm = TRUE),
    score_upper_envelope = max(posterior_score_q975, na.rm = TRUE),
    rank_lower_envelope = min(posterior_rank_q025, na.rm = TRUE),
    rank_upper_envelope = max(posterior_rank_q975, na.rm = TRUE),
    n_model_specifications = .N
  ),
  by = .(candidate_id, gene_symbol)
]

setorder(consensus, -p_rank1, posterior_rank_mean, candidate_id)
consensus[, posterior_rank := seq_len(.N)]
top5 <- consensus[seq_len(min(5L, .N))]
focus <- if (nzchar(FOCUS_GENE)) consensus[toupper(gene_symbol)==toupper(FOCUS_GENE)] else consensus[0L]

fwrite(consensus, file.path(OUT_DIR, "bayesian_dependence_sensitivity_ensemble_consensus.csv"))
fwrite(top5, file.path(OUT_DIR, "bayesian_dependence_sensitivity_ensemble_top5.csv"))
fwrite(focus, file.path(OUT_DIR, "focus_gene_bayesian_dependence_sensitivity_ensemble.csv"))

model_comparison <- combined[
  ,
  .(
    mean_p_rank1 = mean(p_rank1, na.rm = TRUE),
    mean_p_top5 = mean(p_top5, na.rm = TRUE),
    mean_rank = mean(posterior_rank_mean, na.rm = TRUE)
  ),
  by = .(scenario, network_model, gamma)
]
fwrite(model_comparison, file.path(OUT_DIR, "bayesian_dependence_model_comparison.csv"))
fwrite(rbindlist(weight_summary_rows), file.path(OUT_DIR, "dependence_adjusted_weight_summary.csv"))

fwrite(data.table(
  seed = SEED,
  input = INPUT,
  n_candidates = nrow(dat),
  n_draws_per_model = N_DRAWS,
  n_model_specifications = length(all_results),
  base_genomic_weight = BASE_WEIGHTS["genomic"],
  base_hpo_weight = BASE_WEIGHTS["hpo"],
  base_network_weight = BASE_WEIGHTS["network"]
), file.path(OUT_DIR, "run_manifest.csv"))

cat("\n============================================================\n")
cat("DEPENDENCE-AWARE MONTE CARLO RANKING COMPLETE\n")
cat("============================================================\n\n")
print(top5[, .(posterior_rank, candidate_id, gene_symbol, posterior_score_mean, posterior_rank_mean, p_rank1, p_top3, p_top5, p_top10)])
if (nrow(focus)) {
  cat("\nfocus_gene\n")
  print(focus[, .(candidate_id, gene_symbol, posterior_score_mean, posterior_rank_mean, p_rank1, p_top3, p_top5, p_top10, score_lower_envelope, score_upper_envelope, rank_lower_envelope, rank_upper_envelope)])
}
cat("\nOutputs: ", normalizePath(OUT_DIR), "\n", sep = "")

saveRDS(.rng_initial,file.path(OUT_DIR,"rng_state_initial.rds"))
saveRDS(.Random.seed,file.path(OUT_DIR,"rng_state_final.rds"))
a11_finish(.ctx,OUT_DIR,"a11_dependence_sensitivity",list(seed=SEED,draws=N_DRAWS,
 base_weights=as.list(BASE_WEIGHTS),network_models=NETWORK_MODELS,gamma_values=GAMMA_VALUES,
 scenarios=lapply(SCENARIOS,function(s) list(concentration=s$concentration,precision=as.list(s$precision))),
 focus_gene=FOCUS_GENE,minimum_complete_correlation_rows=20L,correlation_shrinkage=.50,
 beta_clamp_epsilon=1e-6,minimum_evidence_views=2L,rank_tie_break="candidate_id",
 consensus="equal arithmetic average across sensitivity specifications"))
