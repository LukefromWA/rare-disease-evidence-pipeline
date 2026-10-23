#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

a <- parse_args(
    c("config", "sample", "c4", "outdir"),
    c("frozen_union_dir", "a10_views", "dependence_views")
)
cfg <- read_config(a$config)
require_bundle(cfg)

frozen <- !is.null(a$frozen_union_dir) && nzchar(a$frozen_union_dir)

if (frozen) {
    ud <- a$frozen_union_dir
    mp <- file.path(ud, "manifest.json")
    mf <- jsonlite::read_json(mp)

    if (!identical(mf$version, "A10_C4_union_Stan_v1_2"))
        stop("Unsupported union version")

    names_needed <- c(
        "current_Stan_final_hypothesis_set.tsv",
        "current_Stan_scenarios.tsv",
        "shortlisted_variant_evidence.tsv",
        "candidate_registry.tsv",
        "protected_unscored_candidates.tsv"
    )
    paths <- file.path(ud, names_needed)

    for (i in seq_along(paths)) {
        path <- paths[i]
        name <- names_needed[i]
        found <- vapply(mf$outputs, function(x) {
            identical(x$path, name)
        }, logical(1))

        if (sum(found) != 1L ||
            !file.exists(path) ||
            !identical(
                mf$outputs[[which(found)]]$sha256,
                sha256_file(path)
            ))
            stop("Frozen union checksum mismatch: ", name)
    }

    if (!identical(
        normalizePath(a$c4),
        normalizePath(paths[1L])
    )) stop("D0 input does not match configured frozen union")

    c4 <- data.table::fread(paths[1L])
    sc <- data.table::fread(paths[2L])

    need_cols(c4, c(
        "hypothesis_id", "gene", "variant", "q_j", "mechanism", "p_m"
    ), "frozen current-Stan hypotheses")
    need_cols(sc, c(
        "scenario_id", "hypothesis_id", "gene", "variant",
        "mechanism", "scenario_weight"
    ), "frozen current-Stan scenarios")

    h <- unique(c4[, .(hypothesis_id, gene, variant, q_j)])
    stopifnot(
        !anyDuplicated(h$hypothesis_id),
        all(is.finite(h$q_j)),
        all(h$q_j >= 0),
        abs(sum(h$q_j) - 1) < 1e-8,
        all(is.finite(c4$p_m)),
        all(c4$p_m >= 0 & c4$p_m <= 1),
        all(abs(c4[, .(mass = sum(p_m)),
                   by = hypothesis_id]$mass - 1) < 1e-8),
        !anyDuplicated(c4[, .(hypothesis_id, mechanism)]),
        !anyDuplicated(sc$scenario_id),
        all(is.finite(sc$scenario_weight)),
        all(sc$scenario_weight >= 0),
        abs(sum(sc$scenario_weight) - 1) < 1e-8,
        setequal(sc$hypothesis_id, h$hypothesis_id)
    )

    c4_key <- paste(c4$hypothesis_id, c4$mechanism, sep = "::")
    sc_key <- paste(sc$hypothesis_id, sc$mechanism, sep = "::")
    index <- match(sc_key, c4_key)

    if (anyNA(index)) stop("Scenario does not match a frozen hypothesis")
    if (!all(sc$gene == c4$gene[index]) ||
        !all(sc$variant == c4$variant[index]))
        stop("Frozen scenario gene/variant mismatch")

    expected <- c4$q_j[index] * c4$p_m[index]
    if (sum(expected) <= 0) stop("Zero frozen scenario mass")
    expected <- expected / sum(expected)

    if (any(abs(expected - sc$scenario_weight) > 1e-8))
        stop("Frozen scenario weights disagree with q_j * p_m")

    if (nrow(h) != mf$counts$primary_current_Stan ||
        nrow(sc) != mf$counts$primary_scenarios)
        stop("Frozen handoff counts disagree with manifest")

    # Preserve existing scenario IDs, weights and metadata unchanged.
    audit <- data.table::data.table(
        status = "consumed_verified_frozen_union",
        weighted_hypotheses = nrow(h),
        scenarios = nrow(sc),
        protected_unscored = mf$counts$protected_unscored,
        supporting_variant_rows = mf$counts$shortlisted_evidence_rows
    )
    inputs <- c(mp, paths)
} else {
    c4 <- read_gz(a$c4)
    need_cols(c4, c(
        "gene", "variant", "q_j", "mechanism", "p_m"
    ), "C4 hypothesis set")
    ex <- expand_scenarios(
        unique(c4[, .(gene, variant, q_j)]),
        unique(c4[, .(variant, mechanism, p_m)])
    )
    sc <- ex$scenarios
    audit <- ex$audit
    inputs <- a$c4
}

stopifnot(
    abs(sum(sc$scenario_weight) - 1) < 1e-8,
    !anyDuplicated(sc$scenario_id)
)


# BEGIN track2_model_views_v1
if (!frozen)
    stop("This multi-view attachment requires the configured frozen union")

need_cols(sc, c(
    "hypothesis_id", "gene", "variant",
    "a10_variant_original", "q_j",
    "stan_probability", "stan_rank_mean_probability"
), "frozen model-view source")

h <- unique(sc[, .(
    hypothesis_id, gene, variant, a10_variant_original,
    q_j, stan_probability, stan_rank_mean_probability
)])
if (anyDuplicated(h$hypothesis_id))
    stop("Inconsistent hypothesis metadata across scenarios")
if (any(abs(h$q_j - h$stan_probability) > 1e-8))
    stop("q_j is not the declared Stan-view weight")

make_view <- function(view, type, score, rank, status,
                      source_field, source_hash,
                      p_rank1 = rep(NA_real_, nrow(h)),
                      p_top3 = rep(NA_real_, nrow(h)),
                      p_top5 = rep(NA_real_, nrow(h)),
                      p_top10 = rep(NA_real_, nrow(h))) {
    data.table::data.table(
        hypothesis_id = h$hypothesis_id,
        gene = h$gene,
        variant = h$variant,
        model_view = view,
        view_status = status,
        score_type = type,
        score = score,
        rank = rank,
        p_rank1 = p_rank1,
        p_top3 = p_top3,
        p_top5 = p_top5,
        p_top10 = p_top10,
        source_field = source_field,
        source_sha256 = source_hash
    )
}

view_tables <- list(make_view(
    "stan",
    "normalized_current_Stan_hypothesis_weight",
    h$stan_probability,
    h$stan_rank_mean_probability,
    rep("scored", nrow(h)),
    "stan_probability",
    sha256_file(a$c4)
))

if (is.null(a$a10_views) || is.null(a$dependence_views))
    stop("Configured multi-view run requires A10 and dependence source files")

av <- read_gz(a$a10_views)
dv <- data.table::fread(a$dependence_views)

need_cols(av, c(
    "best_candidate_id", "gene_symbol",
    "posterior_informed", "posterior_agnostic",
    "rank_informed", "rank_agnostic"
), "A10 model views")
need_cols(dv, c(
    "candidate_id", "gene_symbol", "posterior_score_mean",
    "posterior_rank", "p_rank1", "p_top3", "p_top5", "p_top10"
), "dependence ensemble")

if (anyDuplicated(av$best_candidate_id) ||
    anyDuplicated(dv$candidate_id))
    stop("Duplicate source candidate IDs")

ai <- match(h$a10_variant_original, av$best_candidate_id)
if (anyNA(ai) || !all(h$gene == av$gene_symbol[ai]))
    stop("A10 identities do not match the frozen representatives")
if (!setequal(h$a10_variant_original, av$best_candidate_id))
    stop("A10 candidate universe differs from the frozen representatives")

di <- match(h$a10_variant_original, dv$candidate_id)
matched <- !is.na(di)
if (!all(h$gene[matched] == dv$gene_symbol[di[matched]]))
    stop("Dependence ensemble gene mismatch")
if (!all(dv$candidate_id %in% h$a10_variant_original))
    stop("Dependence source contains unmatched candidate IDs")

for (kind in c("informed", "agnostic")) {
    field <- paste0("posterior_", kind)
    rank_field <- paste0("rank_", kind)
    values <- as.numeric(av[[field]][ai])

    if (any(!is.finite(values)) || any(values < 0))
        stop("Invalid A10 source values")

    view_tables[[length(view_tables) + 1L]] <- make_view(
        paste0("a10_", kind),
        "original_A10_representative_weight_not_renormalized",
        values,
        as.numeric(av[[rank_field]][ai]),
        rep("scored", nrow(h)),
        field,
        sha256_file(a$a10_views)
    )
}

view_tables[[length(view_tables) + 1L]] <- make_view(
    "dependence_ensemble",
    "posterior_score_mean_not_a_normalized_hypothesis_probability",
    as.numeric(dv$posterior_score_mean[di]),
    as.numeric(dv$posterior_rank[di]),
    ifelse(matched, "scored", "not_scored_in_source"),
    "posterior_score_mean",
    sha256_file(a$dependence_views),
    p_rank1 = as.numeric(dv$p_rank1[di]),
    p_top3 = as.numeric(dv$p_top3[di]),
    p_top5 = as.numeric(dv$p_top5[di]),
    p_top10 = as.numeric(dv$p_top10[di])
)

model_views <- data.table::rbindlist(view_tables)
stopifnot(!anyDuplicated(
    model_views[, .(hypothesis_id, model_view)]
))

view_coverage <- model_views[, .(
    hypotheses = .N,
    scored = sum(view_status == "scored"),
    not_scored = sum(view_status != "scored"),
    source_value_sum_diagnostic_only = sum(score, na.rm = TRUE)
), by = .(model_view, score_type)]

inputs <- c(inputs, a$a10_views, a$dependence_views)
# END track2_model_views_v1


od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(sc, file.path(od, "scenarios.tsv.gz"))
write_table(audit, file.path(od, "scenario_audit.tsv"))

write_table(model_views, file.path(od, "candidate_model_views.tsv.gz"))
write_table(view_coverage, file.path(od, "model_view_coverage.tsv"))

write_node_manifest(
    od, "d0", cfg,
    inputs = inputs,
    outputs = file.path(od, c("scenarios.tsv.gz", "scenario_audit.tsv", "candidate_model_views.tsv.gz", "model_view_coverage.tsv")),
    extra = list(
        frozen_union = frozen,
        q_j_scope = "current_Stan_view_only_not_cross_model_consensus",
        model_views_mixed_into_new_probability = FALSE,
        n_scenarios = nrow(sc),
        weights_sum = sum(sc$scenario_weight),
        supporting_variant_evidence_preserved_as_input = frozen,
        protected_unscored_candidates_assigned_weights = FALSE
    )
)
logmsg(
    "D0: frozen_union=", frozen,
    "; scenarios=", nrow(sc),
    "; scenario-weight sum=", sum(sc$scenario_weight)
)
