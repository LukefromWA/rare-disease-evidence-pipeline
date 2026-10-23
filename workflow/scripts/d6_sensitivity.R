#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

# D6: full-universe weight sensitivity, not treatment-benefit calibration.
a <- parse_args(c(
    "config", "sample", "gated", "scenarios", "outdir"
))
cfg <- read_config(a$config)
d <- read_gz(a$gated)
sc <- read_gz(a$scenarios)

need_cols(sc, c(
    "scenario_id", "scenario_weight"
), "canonical D0 scenarios")
need_cols(d, c(
    "scenario_id", "drug_id", "eligibility"
), "gated evidence")

if (!nrow(sc) || anyNA(sc$scenario_id) ||
    anyDuplicated(sc$scenario_id) ||
    !is.numeric(sc$scenario_weight) ||
    any(!is.finite(sc$scenario_weight)) ||
    any(sc$scenario_weight < 0) ||
    abs(sum(sc$scenario_weight) - 1) > 1e-8)
    stop("Invalid canonical full scenario universe")

if (anyNA(d$scenario_id) || anyNA(d$drug_id) ||
    !all(d$scenario_id %in% sc$scenario_id))
    stop("Gated evidence contains missing or unknown identifiers")

if ("scenario_weight" %in% names(d)) {
    canonical_weights <- sc$scenario_weight[
        match(d$scenario_id, sc$scenario_id)
    ]
    if (!is.numeric(d$scenario_weight) ||
        any(!is.finite(d$scenario_weight)) ||
        any(abs(d$scenario_weight - canonical_weights) > 1e-10))
        stop("Gated weights differ from canonical D0 weights")
}

g <- d[eligibility == "eligible_for_research_tiering"]
counts <- d[, .(
    n_drugs_in_gated_table = uniqueN(drug_id),
    n_eligible_drugs = uniqueN(
        drug_id[which(eligibility == "eligible_for_research_tiering")]
    )
), by = scenario_id]

coverage <- merge(
    sc[, .(scenario_id, scenario_weight)],
    counts,
    by = "scenario_id", all.x = TRUE, sort = FALSE
)
for (column in c("n_drugs_in_gated_table", "n_eligible_drugs")) {
    set(
        coverage,
        i = which(is.na(coverage[[column]])),
        j = column,
        value = 0L
    )
}
coverage[, coverage_status := fcase(
    n_drugs_in_gated_table == 0,
        "no_rows_in_direct_gated_table",
    n_eligible_drugs == 0,
        "rows_present_but_none_eligible",
    default = "eligible_rows_present"
)]

if (nrow(g)) {
    result <- track2_sensitivity(
        g,
        unlist(cfg$track2$family_weights),
        cfg$track2$mc_draws,
        cfg$tier15$kappa %||% 100,
        cfg$track2$top_k,
        cfg$run$seed,
        scenarios = sc
    )
    state <- "full_universe_weight_sensitivity_completed"
} else {
    result <- data.table(
        drug_id = character(),
        median_score = numeric(),
        score_lo = numeric(),
        score_hi = numeric(),
        top_k_prob = numeric(),
        median_rank = numeric(),
        best_rank = numeric(),
        worst_rank = numeric(),
        baseline_score = numeric(),
        scenario_evidence_weight = numeric(),
        uncovered_scenario_weight = numeric(),
        n_full_scenarios = integer()
    )
    state <- "no_eligible_drugs"
}

od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(result, file.path(od, "sensitivity.tsv"))
write_table(coverage, file.path(od, "scenario_coverage.tsv"))
write_table(data.table(
    status = state,
    n_eligible_drugs = uniqueN(g$drug_id),
    n_full_scenarios = nrow(sc),
    full_scenario_weight = sum(sc$scenario_weight),
    weight_with_eligible_rows = sum(
        coverage[n_eligible_drugs > 0, scenario_weight]
    ),
    weight_without_eligible_rows = sum(
        coverage[n_eligible_drugs == 0, scenario_weight]
    )
), file.path(od, "status.tsv"))

write_node_manifest(
    od, "d6", cfg,
    inputs = c(a$gated, a$scenarios),
    outputs = file.path(od, c(
        "sensitivity.tsv", "scenario_coverage.tsv", "status.tsv"
    )),
    extra = list(
        status = state,
        top_k = cfg$track2$top_k,
        draws = cfg$track2$mc_draws,
        weight_scope =
            "declared_D0_weights_not_cross_model_consensus",
        observed_scenario_subset_renormalized = FALSE,
        absent_pair_convention =
            "zero_evidence_score_contribution_not_no_biological_effect",
        interpretation =
            "declared_weight_sensitivity_not_benefit_probability",
        rank_ties = "min",
        zero_weight_scenarios_preserved = TRUE
    )
)
