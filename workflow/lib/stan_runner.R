# R/stan_runner.R -- CmdStanR implementation for M2 grouped OOF and compact freezing.
suppressPackageStartupMessages({ library(data.table); library(cmdstanr) })

STAN_FILES <- c(m2a = "m2a.stan", m2b = "m2b.stan", m3a = "m3a.stan", m3b = "m3b.stan")
MODEL_PARS <- list(m2a = "beta", m2b = c("beta", "rho_raw", "beta_cop"),
                   m3a = c("beta", "z_beta", "tau_beta"),
                   m3b = c("beta", "z_beta", "tau_beta", "mu_rho_raw", "tau_rho", "z_rho", "beta_cop", "tau_cop", "z_cop"))
PROFILES <- list(quick = list(chains = 4, warmup = 1000, sampling = 1000, adapt_delta = 0.99, thin = 1),
                 final = list(chains = 4, warmup = 1000, sampling = 1000, adapt_delta = 0.95, thin = 2))
GATES <- list(quick = list(rhat_max = 1.05, ess_min = 100, ebfmi_min = 0.3), final = list(rhat_max = 1.01, ess_min = 400, ebfmi_min = 0.3))

stan_dir <- function() { d <- Sys.getenv("MVA_STAN_DIR", unset = ""); if (nzchar(d)) d else file.path(getwd(), "stan") }

compile_models <- function(models = names(STAN_FILES), dir = stan_dir(), cache_dir = NULL, verbose = TRUE) {
  out <- list()
  for (id in models) {
    f <- file.path(dir, STAN_FILES[[id]])
    if (!file.exists(f)) stop("missing Stan file: ", f)
    if (verbose) message("compiling ", id, " ...")
    out[[id]] <- cmdstanr::cmdstan_model(f, compile = TRUE)
  }
  out
}

clear_cache <- function(cache_dir) { if (!is.null(cache_dir) && dir.exists(cache_dir)) unlink(cache_dir, recursive = TRUE); invisible(TRUE) }

fit_prep_state <- function(train, feature_cols, hpo_col, net_col) {
  probs <- seq(0, 1, length.out = 513)
  grid <- function(x) { x <- x[is.finite(x)]; list(q = quantile(x, probs, names = FALSE, type = 7), n = length(x)) }
  list(mu = vapply(train[, ..feature_cols], mean, 0), sd = pmax(vapply(train[, ..feature_cols], sd, 0), 1e-8), feature_cols = feature_cols,
       probs = probs, hpo = grid(train[[hpo_col]]), net = grid(train[[net_col]]), hpo_col = hpo_col, net_col = net_col)
}

apply_prep_state <- function(dt, st) {
  d <- copy(as.data.table(dt))
  for (f in st$feature_cols) set(d, j = f, value = (d[[f]] - st$mu[[f]]) / st$sd[[f]])
  u <- function(x, g) { v <- approx(g$q, st$probs, xout = x, ties = mean, rule = 2)$y * g$n / (g$n + 1); pmin(pmax(v, 1e-6), 1 - 1e-6) }
  ok <- is.finite(d[[st$hpo_col]]) & is.finite(d[[st$net_col]])
  d[, copula_available := as.integer(ok)]
  d[, u_hpo := fifelse(ok, u(d[[st$hpo_col]], st$hpo), 0.5)][, u_network := fifelse(ok, u(d[[st$net_col]], st$net), 0.5)]
  d[]
}

default_priors <- function(K) list(beta_prior_scale = rep(1, K), tau_beta_prior_scale = rep(0.5, K), mu_rho_prior_scale = 1, tau_rho_prior_scale = 0.5,
                                   beta_copula_prior_scale = 1, tau_copula_prior_scale = 0.5)

SCALAR_DATA <- c("N", "C", "K", "G", "N_new", "C_new", "mu_rho_prior_scale", "tau_rho_prior_scale", "beta_copula_prior_scale", "tau_copula_prior_scale")
as_stan_arrays <- function(dat) { nm <- setdiff(names(dat), SCALAR_DATA)
  dat[nm] <- lapply(dat[nm], function(x) if (is.matrix(x)) x else as.array(x)); dat }

.group_index <- function(g, levels, allow_new) { i <- match(g, levels); if (anyNA(i)) { if (!allow_new) stop("unknown group in training data"); i[is.na(i)] <- length(levels) + 1L }; i }

pack_stan_data <- function(train, valid = NULL, feature_cols, group_levels, prior = default_priors(length(feature_cols))) {
  setorder(train, case_id, candidate_id); train[, .row := .I]
  cs <- train[, .(s = min(.row), e = max(.row), npos = sum(label == 1L), y = which(label == 1L)[1L]), by = case_id]
  if (any(cs$npos != 1L)) stop("each training case needs exactly one positive candidate")
  dat <- c(list(N = nrow(train), C = nrow(cs), K = length(feature_cols), G = length(group_levels), X = as.matrix(train[, ..feature_cols]),
                group_id = .group_index(train$group, group_levels, FALSE), u_hpo = train$u_hpo, u_network = train$u_network,
                copula_available = train$copula_available, case_start = cs$s, case_end = cs$e, y_pos = cs$y), prior)
  train[, .row := NULL]
  if (is.null(valid)) return(as_stan_arrays(c(dat, list(N_new = 0L, C_new = 0L, X_new = matrix(0, 0, length(feature_cols)), group_new = integer(), u_hpo_new = numeric(),
                                         u_network_new = numeric(), copula_available_new = integer(), case_start_new = integer(), case_end_new = integer()))))
  setorder(valid, case_id, candidate_id); valid[, .row := .I]; vs <- valid[, .(s = min(.row), e = max(.row)), by = case_id]; valid[, .row := NULL]
  as_stan_arrays(c(dat, list(N_new = nrow(valid), C_new = nrow(vs), X_new = as.matrix(valid[, ..feature_cols]), group_new = .group_index(valid$group, group_levels, TRUE),
              u_hpo_new = valid$u_hpo, u_network_new = valid$u_network, copula_available_new = valid$copula_available, case_start_new = vs$s, case_end_new = vs$e)))
}

fit_model <- function(model, id, data, profile = "quick", seed = 1L) {
  pr <- PROFILES[[profile]]
  model$sample(data = data, iter_warmup = pr$warmup, iter_sampling = pr$sampling, chains = pr$chains,
               parallel_chains = min(pr$chains, 2L), thin = pr$thin, adapt_delta = pr$adapt_delta,
               seed = seed, refresh = 100, show_messages = TRUE, save_warmup = FALSE)
}

sampler_diag <- function(fit) {
  d <- fit$diagnostic_summary()
  list(divergences = sum(d$num_divergent), max_treedepth_hits = sum(d$num_max_treedepth),
       min_ebfmi = min(unlist(d$ebfmi)))
}

diagnostic_gate <- function(fit, pars, profile = "final") {
  g <- GATES[[profile]]
  draws <- fit$draws(variables = pars)
  rhat <- posterior::rhat(draws)
  ess <- posterior::ess_bulk(draws)
  d <- sampler_diag(fit)
  chk <- c(rhat = all(rhat < g$rhat_max, na.rm = TRUE),
           ess = all(ess > g$ess_min, na.rm = TRUE),
           div = d$divergences == 0,
           ebfmi = d$min_ebfmi > g$ebfmi_min,
           depth = d$max_treedepth_hits == 0)
  list(pass = all(chk), checks = chk, max_rhat = max(rhat, na.rm = TRUE),
       min_ess = min(ess, na.rm = TRUE), diag = d)
}

param_summary <- function(fit, pars, id) {
  dm <- fit$draws(variables = pars, format = "draws_matrix")
  draws <- fit$draws(variables = pars)
  q <- apply(dm, 2, quantile, probs = c(0.025, 0.5, 0.975), names = FALSE)
  rhat <- posterior::rhat(draws)
  ess <- posterior::ess_bulk(draws)
  data.table(model_id = id, parameter = colnames(dm), mean = signif(colMeans(dm), 4),
             sd = signif(apply(dm, 2, sd), 4), q025 = signif(q[1, ], 4),
             q50 = signif(q[2, ], 4), q975 = signif(q[3, ], 4),
             n_eff = round(ess), Rhat = round(rhat, 3))
}

posterior_mean_probs <- function(fit) {
  d <- fit$draws(variables = "p_new", format = "draws_matrix")
  colMeans(d)
}


# ---- adaptive refit: escalate warmup/sampling and adapt_delta when diagnostics fail ----
fit_model_adaptive <- function(model, id, data, profile = "quick", seed = 1L, max_attempts = 3L) {
  attempts <- list(
    list(warmup = 1000L, sampling = 1000L, adapt_delta = 0.99),
    list(warmup = 1500L, sampling = 1500L, adapt_delta = 0.995),
    list(warmup = 2000L, sampling = 2000L, adapt_delta = 0.999)
  )
  pr <- PROFILES[[profile]]
  best <- NULL
  for (i in seq_len(max_attempts)) {
    a <- attempts[[min(i, length(attempts))]]
    message(sprintf("adaptive fit: %s attempt %d/%d (warmup=%d, sampling=%d, adapt_delta=%.3f)",
                    id, i, max_attempts, a$warmup, a$sampling, a$adapt_delta))
    fit <- model$sample(
      data = data,
      iter_warmup = a$warmup,
      iter_sampling = a$sampling,
      chains = pr$chains,
      parallel_chains = min(pr$chains, 4L),
      thin = pr$thin,
      adapt_delta = a$adapt_delta,
      seed = seed + i,
      refresh = 100,
      show_messages = TRUE,
      save_warmup = FALSE
    )
    gate <- diagnostic_gate(fit, MODEL_PARS[[id]], profile)
    if (gate$pass) {
      return(list(fit = fit, gate = gate, attempt = i, passed = TRUE))
    }
    if (is.null(best) || gate$max_rhat < best$gate$max_rhat) {
      best <- list(fit = fit, gate = gate, attempt = i, passed = FALSE)
    }
  }
  best
}

run_oof <- function(dev, feature_cols, group_levels, compiled, models = names(compiled), profile = "quick", seed = 1L, out_dir = NULL,
                    hpo_col = "hpo_score", net_col = "network_score", prior = default_priors(length(feature_cols))) {
  oof <- list(); dg <- list(); ps <- list()
  for (fold in sort(unique(dev$inner_fold))) {
    tr <- copy(dev[inner_fold != fold]); va <- copy(dev[inner_fold == fold])
    st <- fit_prep_state(tr, feature_cols, hpo_col, net_col)
    tr <- apply_prep_state(tr, st); va <- apply_prep_state(va, st)
    dat <- pack_stan_data(tr, va, feature_cols, group_levels, prior)
    setorder(va, case_id, candidate_id)
    for (id in models) {
      adaptive <- fit_model_adaptive(compiled[[id]], id, dat, profile, seed + fold)
      fit <- adaptive$fit
      gate <- adaptive$gate

      chain_dir <- file.path(out_dir, "chains")
      dir.create(chain_dir, recursive = TRUE, showWarnings = FALSE)
      chain_file <- file.path(chain_dir, sprintf("%s_fold%02d.rds", id, fold))
      param_draws <- fit$draws(variables = MODEL_PARS[[id]], format = "draws_matrix")
      p_new_draws <- if (any(grepl("^p_new", colnames(fit$draws(format = "draws_df"))))) {
        fit$draws(variables = "p_new", format = "draws_matrix")
      } else {
        NULL
      }
      saveRDS(list(
        model_id = id,
        inner_fold = fold,
        attempt = adaptive$attempt,
        passed = adaptive$passed,
        parameters = MODEL_PARS[[id]],
        param_draws = param_draws,
        p_new_draws = p_new_draws
      ), chain_file, compress = "gzip")
      p <- posterior_mean_probs(fit)
      oof[[length(oof) + 1L]] <- data.table(case_id = va$case_id, candidate_id = va$candidate_id, label = va$label, model_id = id,
                                             probability = signif(p, 5), inner_fold = fold)
      dg[[length(dg) + 1L]] <- data.table(model_id = id, inner_fold = fold, pass = gate$pass, max_rhat = gate$max_rhat,
                                           min_ess = gate$min_ess, divergences = gate$diag$divergences, min_ebfmi = gate$diag$min_ebfmi)
      ps[[length(ps) + 1L]] <- cbind(param_summary(fit, MODEL_PARS[[id]], id), inner_fold = fold)
      rm(fit); invisible(gc(verbose = FALSE))
    }
  }
  res <- list(oof = rbindlist(oof), diagnostics = rbindlist(dg), param_summaries = rbindlist(ps))
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(res$oof, file.path(out_dir, "oof_predictions.tsv.gz"))
    fwrite(res$diagnostics, file.path(out_dir, "oof_diagnostics.tsv"))
    fwrite(res$param_summaries, file.path(out_dir, "oof_param_summaries.tsv.gz"))
  }
  res
}

eligible_models <- function(diagnostics) diagnostics[, .(ok = all(pass)), by = model_id][ok == TRUE, model_id]

freeze_model <- function(model, id, dev, feature_cols, group_levels, profile = "final", seed = 1L, n_keep = 200L, digits = 4L,
                         hpo_col = "hpo_score", net_col = "network_score", prior = default_priors(length(feature_cols)), out_dir = NULL) {
  st <- fit_prep_state(dev, feature_cols, hpo_col, net_col)
  tr <- apply_prep_state(copy(dev), st)
  fit <- fit_model(model, id, pack_stan_data(tr, NULL, feature_cols, group_levels, prior), profile, seed)
  gate <- diagnostic_gate(fit, MODEL_PARS[[id]], profile)
  M <- fit$draws(variables = MODEL_PARS[[id]], format = "draws_matrix")
  idx <- unique(round(seq(1, nrow(M), length.out = min(n_keep, nrow(M)))))
  draws <- signif(M[idx, , drop = FALSE], digits)
  fr <- list(model_id = id, draws = draws, prep_state = st, group_levels = group_levels, feature_cols = feature_cols, prior = prior,
             gate = gate[c("pass", "checks", "max_rhat", "min_ess")],
             stan_sha = digest::digest(file = file.path(stan_dir(), STAN_FILES[[id]]), algo = "sha256"),
             seed = seed, cmdstanr = as.character(packageVersion("cmdstanr")))
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(fr, file.path(out_dir, paste0("frozen_", id, ".rds")), compress = "xz")
  }
  rm(fit, M); invisible(gc(verbose = FALSE)); fr
}

score_frozen <- function(fr, model, new_dt) {
  va <- apply_prep_state(copy(new_dt), fr$prep_state)
  setorder(va, case_id, candidate_id)
  one <- data.table(case_id = "d", candidate_id = "d", label = 1L, group = fr$group_levels[1], u_hpo = 0.5, u_network = 0.5, copula_available = 0L)
  one[, (fr$feature_cols) := 0]
  dat <- pack_stan_data(one, va, fr$feature_cols, fr$group_levels, fr$prior)
  out <- model$generate_quantities(data = dat, draws = fr$draws)
  data.table(case_id = va$case_id, candidate_id = va$candidate_id, probability = colMeans(out$draws(vars = "p_new", format = "draws_matrix")))
}

bma_sensitivity <- function(dev, feature_cols, group_levels, compiled, models = names(compiled), profile = "final", seed = 1L, reps = 5L, out_dir = NULL,
                            hpo_col = "hpo_score", net_col = "network_score", prior = default_priors(length(feature_cols))) {
  stop("BMA sensitivity is not enabled in the CmdStanR runner; use eligible_models(), stack_from_oof(), and freeze_model() instead.")
}

dir_size_kb <- function(path) round(sum(file.info(list.files(path, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE) / 1024, 1)
