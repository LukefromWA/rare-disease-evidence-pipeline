# workflow/lib/gates.R -- Layer D: scenarios, direction + safety gates, evidence integration, sensitivity
suppressPackageStartupMessages(library(data.table))

expand_scenarios <- function(hyp, mech) {                # join ON the hypothesis; never a cross join
  assert_columns(hyp, c("gene", "variant", "q_j"), "hypotheses"); assert_columns(mech, c("variant", "mechanism", "p_m"), "mechanisms")
  stopifnot(abs(sum(unique(hyp[, .(variant, q_j)])$q_j) - 1) < 1e-8)
  s <- merge(hyp, mech, by = "variant"); s[, `:=`(scenario_id = paste(gene, variant, mechanism, sep = "|"), w_raw = q_j * p_m)]
  dropped <- s[w_raw <= 0, .(scenario_id, reason = "zero_weight")]; s <- s[w_raw > 0]
  s[, scenario_weight := w_raw / sum(w_raw)]
  list(scenarios = s[], audit = dropped)
}

direction_of <- function(mechanism) data.table::fcase(
  mechanism %in% c("loss_of_function", "haploinsufficiency", "biallelic_loss"), "restore_or_compensate",
  mechanism %in% c("gain_of_function", "activating"), "inhibit_or_counteract",
  mechanism == "dominant_negative", "reduce_mutant_effect_or_bypass", default = "uncertain_direction")

# uses the NORMALIZED action vocabulary from the resource bundle: activating / inhibitory / uncertain
direction_gate <- function(direction, action_normalized) {
    d <- tolower(trimws(as.character(direction)))
    a <- tolower(trimws(as.character(action_normalized)))

    supported_direction <- d %in% c(
        "restore_or_compensate", "inhibit_or_counteract"
    )
    supported_action <- a %in% c("activating", "inhibitory")

    compatible <- (
        d %in% "restore_or_compensate" & a %in% "activating"
    ) | (
        d %in% "inhibit_or_counteract" & a %in% "inhibitory"
    )

    # Dominant-negative rescue needs intervention-specific evidence.
    # Generic inhibition is not treated as mutant-selective suppression.
    data.table::fcase(
        compatible, "pass",
        supported_direction & supported_action, "fail",
        default = "uncertain"
    )
}

safety_gate <- function(ped_status, contra, route_feasible) {
    ped <- tolower(trimws(as.character(ped_status)))
    con <- tolower(trimws(as.character(contra)))
    route <- toupper(trimws(as.character(route_feasible)))

    fail <- ped %in% "contraindicated" |
        con %in% "absolute" |
        route %in% c("FALSE", "0")

    complete <- ped %in% c("supportive", "limited") &
        con %in% c("none", "relative") &
        route %in% c("TRUE", "1")

    data.table::fcase(
        fail, "fail",
        complete, "pass",
        default = "not_configured"
    )
}
final_eligibility <- function(direction_status, safety_status, pediatric_mode = TRUE) data.table::fcase(
  direction_status == "fail" | safety_status == "fail", "excluded",
  direction_status == "pass" & safety_status == "pass", "eligible_for_research_tiering",
  safety_status == "not_configured" & pediatric_mode, "uncertain_hold", default = "uncertain")

# D4: coverage-aware mean over evidence families; keep n_families_available
integrate_track2_evidence <- function(E, w) { s <- weighted_available_mean(E[, names(w), drop = FALSE], w)
  data.table(score = as.numeric(s), n_families_available = attr(s, "n_views")) }

# D6: ALL drugs are scored per draw, then ranked ACROSS drugs within that draw
scenario_sensitivity <- function(M, w0, B = 500L, kappa = 100, top_k = 2L, seed) {
  stopifnot(ncol(M) == length(w0), abs(sum(w0) - 1) < 1e-8); set.seed(seed)    # M: drugs x scenarios (0 = no evidence)
  W <- rdirichlet_base(B, kappa * w0); S <- W %*% t(M)                         # B x drugs
  R <- t(apply(-S, 1, rank, ties.method = "min")); colnames(R) <- rownames(M)
  data.table(drug_id = rownames(M), top_k_prob = colMeans(R <= top_k), median_rank = apply(R, 2, median),
             best_rank = apply(R, 2, min), worst_rank = apply(R, 2, max))
}

# ---- C4 / D helpers ----
consequence_class <- function(consequence, impact) {
  x <- tolower(as.character(consequence)); imp <- toupper(as.character(impact))
  fcase(grepl("stop_gained|frameshift|splice_donor|splice_acceptor|start_lost|stop_lost", x), "lof_like", grepl("missense|inframe", x), "missense",
        is.na(consequence) & imp == "HIGH", "lof_like", is.na(consequence) & imp == "MODERATE", "missense", default = "other")
}

mechanism_prior_table <- function(mp, classes_present) {
  if (identical(mp$mode, "declared")) { d <- mp$declared; if (is.null(d)) stop("mechanism_prior.mode is declared but no table given")
    bad <- vapply(d, function(r) abs(sum(unlist(r)) - 1) > 1e-8, TRUE); if (any(bad)) stop("declared mechanism_prior rows must sum to 1: ", paste(names(d)[bad], collapse = ", "))
    return(rbindlist(lapply(names(d), function(k) data.table(cclass = k, mechanism = names(d[[k]]), p_m = as.numeric(unlist(d[[k]]))))))
  }
  warning("mechanism_prior.mode = uniform: neutral default, NOT a declared biological prior (D-3)", call. = FALSE)
  cl <- unlist(mp$classes); rbindlist(lapply(classes_present, function(k) data.table(cclass = k, mechanism = cl, p_m = 1 / length(cl))))
}

propagate_state <- function(A, x0, alpha, iterations = 40L) { x <- x0; for (i in seq_len(iterations)) x <- tanh(alpha * as.numeric(Matrix::t(A) %*% x) + x0); x }

rescue_mc <- function(g, target, gene, mech_dir, drug_sign, mutation_sign, alpha, draws) {
  nodes <- igraph::V(g)$name
  if (!(target %in% nodes) || !(gene %in% nodes)) return(list(status = "not_in_graph"))
  fwd <- igraph::subcomponent(g, target, mode = "out")$name; bwd <- igraph::subcomponent(g, gene, mode = "in")$name; keep <- intersect(fwd, bwd)
  if (!(gene %in% fwd)) return(list(status = "no_path"))
  sg <- igraph::induced_subgraph(g, keep); nm <- igraph::V(sg)$name; n <- length(nm); el <- igraph::as_edgelist(sg)
  sgn <- igraph::E(sg)$sign; pr <- igraph::E(sg)$edge_probability; dirsign <- if (mech_dir == "restore_or_compensate") 1 else -1
  x0m <- numeric(n); x0m[match(gene, nm)] <- mutation_sign; x0d <- x0m; x0d[match(target, nm)] <- x0d[match(target, nm)] + drug_sign
  gi <- match(gene, nm); res <- numeric(draws)
  for (b in seq_len(draws)) { inc <- if (length(pr)) stats::rbinom(length(pr), 1, pr) else integer()
    A <- if (length(pr)) Matrix::sparseMatrix(match(el[, 1], nm), match(el[, 2], nm), x = sgn * inc, dims = c(n, n)) else Matrix::Matrix(0, n, n, sparse = TRUE)
    res[b] <- dirsign * (propagate_state(A, x0d, alpha)[gi] - propagate_state(A, x0m, alpha)[gi]) }
  list(status = "ok", mean = mean(res), lo = unname(stats::quantile(res, 0.025)), hi = unname(stats::quantile(res, 0.975)), p_gt0 = mean(res > 0), n_nodes = n)
}

normalize01 <- function(x) { r <- range(x, na.rm = TRUE); if (!all(is.finite(r))) return(rep(NA_real_, length(x))); if (diff(r) == 0) return(rep(0.5, length(x))); (x - r[1]) / diff(r) }

# D6 -- scenario AND family-weight perturbation; all drugs scored per draw, ranked across drugs within the draw


# BEGIN full_universe_sensitivity_v1
track2_sensitivity <- function(
    rows, fam_w, B, kappa, top_k, seed, scenarios
) {
    need_cols(scenarios, c(
        "scenario_id", "scenario_weight"
    ), "full sensitivity scenario universe")
    need_cols(rows, c(
        "drug_id", "scenario_id"
    ), "sensitivity rows")

    sc <- data.table::copy(scenarios[, .(
        scenario_id, scenario_weight
    )])
    if (!nrow(sc) || anyNA(sc$scenario_id) ||
        anyDuplicated(sc$scenario_id))
        stop("Invalid full scenario identifiers")
    if (!is.numeric(sc$scenario_weight) ||
        any(!is.finite(sc$scenario_weight)) ||
        any(sc$scenario_weight < 0) ||
        abs(sum(sc$scenario_weight) - 1) > 1e-8)
        stop("Full scenario weights must be nonnegative and sum to one")

    data.table::setorder(sc, scenario_id)

    if (!nrow(rows))
        stop("No rows supplied to sensitivity helper")
    if (anyNA(rows$drug_id) || anyNA(rows$scenario_id))
        stop("Missing drug/scenario identifiers")
    if (anyDuplicated(rows[, .(drug_id, scenario_id)]))
        stop("Duplicate drug/scenario rows: resolve upstream; no overwrite")
    if (!all(rows$scenario_id %in% sc$scenario_id))
        stop("Evidence contains scenarios outside the full universe")

    fam_w <- unlist(fam_w)
    if (!is.numeric(fam_w) ||
        is.null(names(fam_w)) ||
        anyNA(names(fam_w)) ||
        anyDuplicated(names(fam_w)) ||
        any(!nzchar(names(fam_w))) ||
        any(!is.finite(fam_w)) ||
        any(fam_w < 0) ||
        sum(fam_w) <= 0)
        stop("Invalid family weights")

    fam_w <- fam_w / sum(fam_w)
    fams <- names(fam_w)
    evidence_columns <- paste0("e_", fams)
    need_cols(rows, evidence_columns, "family evidence")

    if (length(B) != 1L || !is.finite(B) ||
        B < 1 || B != as.integer(B) ||
        length(top_k) != 1L || !is.finite(top_k) ||
        top_k < 1 || top_k != as.integer(top_k) ||
        length(kappa) != 1L || !is.finite(kappa) || kappa <= 0)
        stop("Invalid sensitivity settings")

    B <- as.integer(B)
    top_k <- as.integer(top_k)

    drugs <- sort(unique(rows$drug_id))
    di <- match(rows$drug_id, drugs)
    si <- match(rows$scenario_id, sc$scenario_id)

    if ("scenario_weight" %in% names(rows)) {
        if (!is.numeric(rows$scenario_weight) ||
            any(!is.finite(rows$scenario_weight)) ||
            any(abs(
                rows$scenario_weight - sc$scenario_weight[si]
            ) > 1e-10))
            stop("Evidence-row weights differ from canonical scenario weights")
    }

    E <- as.matrix(rows[, evidence_columns, with = FALSE])
    storage.mode(E) <- "double"
    colnames(E) <- fams
    if (any(!is.na(E) & !is.finite(E)))
        stop("Nonfinite family evidence")
    if (any(E < 0 | E > 1, na.rm = TRUE))
        stop("Family evidence must be in [0,1] or missing")

    has_evidence <- rowSums(
        !is.na(E[, fam_w > 0, drop = FALSE])
    ) > 0

    sample_weights <- function(weights) {
        active <- weights > 0
        result <- matrix(
            0, nrow = B, ncol = length(weights)
        )
        result[, active] <- rdirichlet_base(
            B, kappa * weights[active]
        )
        colnames(result) <- names(weights)

        if (any(!is.finite(result)) ||
            any(result < 0) ||
            any(abs(rowSums(result) - 1) > 1e-8))
            stop("Invalid Dirichlet draws; no replacement values invented")
        result
    }

    set.seed(seed)
    Ws <- sample_weights(sc$scenario_weight)
    Wf <- sample_weights(fam_w)

    baseline_rows <- as.numeric(
        weighted_available_mean(E, fam_w)
    )
    if (any(!is.finite(baseline_rows[has_evidence])))
        stop("Invalid baseline evidence aggregation")
    baseline_rows[!has_evidence] <- 0

    M0 <- matrix(0, length(drugs), nrow(sc))
    M0[cbind(di, si)] <- baseline_rows
    baseline_score <- as.numeric(M0 %*% sc$scenario_weight)

    coverage <- numeric(length(drugs))
    for (j in seq_along(drugs)) {
        covered_scenarios <- unique(si[di == j & has_evidence])
        coverage[j] <- sum(sc$scenario_weight[covered_scenarios])
    }

    S <- matrix(
        NA_real_, B, length(drugs),
        dimnames = list(NULL, drugs)
    )
    R <- S

    for (draw in seq_len(B)) {
        draw_family_weights <- stats::setNames(
            as.numeric(Wf[draw, ]),
            fams
        )
        values <- as.numeric(
            weighted_available_mean(E, draw_family_weights)
        )
        if (any(!is.finite(values[has_evidence])))
            stop("Invalid sampled evidence aggregation")
        values[!has_evidence] <- 0

        M <- matrix(0, length(drugs), nrow(sc))
        M[cbind(di, si)] <- values
        scores <- as.numeric(M %*% Ws[draw, ])

        S[draw, ] <- scores
        R[draw, ] <- rank(-scores, ties.method = "min")
    }

    data.table::data.table(
        drug_id = drugs,
        median_score = apply(S, 2, median),
        score_lo = apply(S, 2, quantile, 0.025),
        score_hi = apply(S, 2, quantile, 0.975),
        top_k_prob = colMeans(R <= top_k),
        median_rank = apply(R, 2, median),
        best_rank = apply(R, 2, min),
        worst_rank = apply(R, 2, max),
        baseline_score = baseline_score,
        scenario_evidence_weight = coverage,
        uncovered_scenario_weight = pmax(
            0, sum(sc$scenario_weight) - coverage
        ),
        n_full_scenarios = nrow(sc)
    )
}
# END full_universe_sensitivity_v1

# explicit statuses: a method without data returns a NAMED status, never invented evidence
write_method_status <- function(od, method, status, reason, cfg, inputs = character(), node = method) {
  stopifnot(status %in% c("ok", "unavailable_data", "not_applicable", "disabled_pending_enablement", "disabled_pending_rescue_validation", "not_run_cut_list"))
  dir.create(od, recursive = TRUE, showWarnings = FALSE); f <- file.path(od, paste0(method, "_status.tsv"))
  write_table(data.table(method = method, status = status, reason = reason), f); write_node_manifest(od, node, cfg, inputs = inputs, outputs = f, file = paste0("manifest_", method, ".json"), extra = list(status = status)); invisible(f)
}
bundle_paths <- function(cfg) { root <- rp(cfg, cfg$track2$root); d <- file.path(root, "derived")
  list(root = root, ready = file.path(root, ".ready"), enabled = file.path(root, ".signed_propagation_enabled"), drugs = file.path(d, "drug_target_action_map.tsv.gz"),
       safety = file.path(d, "pediatric_safety_map.tsv.gz"), pathways = file.path(d, "reactome_gene_sets.tsv.gz"), graph = file.path(root, "graph", "signed_directed_graph.rds")) }
require_bundle <- function(cfg) {
    b <- bundle_paths(cfg)
    expected <- cfg$track2$resource_scope
    if (is.null(expected)) expected <- "core"

    if (!expected %in% c("core", "full"))
        stop("track2.resource_scope must be core or full")

    if (!file.exists(b$ready))
        stop("Resource readiness sentinel missing: ", b$ready)

    manifest_path <- file.path(
        dirname(b$ready), "manifests", "resource_manifest.json"
    )
    if (!file.exists(manifest_path))
        stop("Resource validation manifest missing")

    manifest <- jsonlite::read_json(manifest_path)

    if (is.null(manifest$scope) ||
        !identical(manifest$scope, expected))
        stop("Resource scope mismatch or old manifest; rerun E22")

    if (!length(manifest$artifacts))
        stop("Resource manifest has no artifacts")

    for (entry in manifest$artifacts) {
        if (!file.exists(entry$path) ||
            !identical(entry$sha256, sha256_file(entry$path)))
            stop("Validated resource changed: ", entry$path)
    }

    b
}
