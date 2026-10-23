#!/usr/bin/env Rscript
# Install in written_final_scripts/workflow/scripts beside c4_hypotheses.R.
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("config", "sample", "a10", "a2", "legacy_c4", "features", "predictions", "outdir"))
for (pkg in c("data.table", "digest", "jsonlite")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
suppressPackageStartupMessages(library(data.table))
sha <- function(p) digest::digest(file = p, algo = "sha256")
need <- function(x, n, label) if (!all(n %in% names(x))) stop(label, " missing: ", paste(setdiff(n, names(x)), collapse = ", "))
blank <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", ".", "NA", "NaN")
key <- function(g, v) paste(g, v, sep = "\r")
truth <- function(x) toupper(trimws(as.character(x))) %in% c("TRUE", "T", "1")
read_tab <- function(p) fread(cmd = if (grepl("[.]gz$", p)) paste("gzip -cd", shQuote(p)) else paste("cat", shQuote(p)), colClasses = "character", na.strings = NULL)
if (file.exists(a$outdir)) stop("Refusing to reuse output directory: ", a$outdir)
input_paths <- unlist(a[c("config", "a10", "a2", "legacy_c4", "features", "predictions")])
if (any(!file.exists(input_paths))) stop("Missing input file")
input_hashes <- vapply(input_paths, sha, character(1))
h <- paste0(a$predictions, ".sha256")
if (!file.exists(h) || !identical(trimws(readLines(h, n = 1L)), sha(a$predictions))) stop("Prediction checksum missing or mismatched")
cfg <- read_config(a$config)
x <- read_tab(a$a10)
a2 <- read_tab(a$a2)
old <- read_tab(a$legacy_c4)
f <- read_tab(a$features)
pred <- readRDS(a$predictions)
s <- as.data.table(copy(pred$summary))
need(x, c("gene_symbol", "best_candidate_id", "best_integrated_score", "best_IMPACT"), "A10")
need(a2, c("candidate_id", "gene_symbol", "consequence", "IMPACT", "pass_core_qc"), "A2")
need(old, c("gene", "variant", "q_j", "mechanism", "p_m"), "Historical C4")
need(f, c("case_id", "candidate_id", "gene", "G", "H", "N"), "C1")
need(s, c("case_id", "candidate_id", "probability", "rank_mean_probability"), "Stan predictions")
if (is.null(pred$predictor_id) || !nzchar(pred$predictor_id)) stop("Prediction artifact lacks predictor_id")
if (anyDuplicated(x$gene_symbol) || anyDuplicated(x$best_candidate_id) || any(blank(x$gene_symbol)) || any(blank(x$best_candidate_id))) stop("A10 keys must be unique and nonblank")
if (anyDuplicated(a2$candidate_id) || anyDuplicated(f$candidate_id) || anyDuplicated(s$candidate_id)) stop("Duplicate allele/representative IDs; inspect before joining")
if (any(blank(old$gene)) || any(blank(old$variant))) stop("Historical C4 contains blank keys")
if (!nrow(x) || !nrow(s)) stop("Empty current candidate set")
if (uniqueN(f$case_id) != 1L || uniqueN(s$case_id) != 1L || !identical(unique(f$case_id), as.character(unique(s$case_id)))) stop("C1/Stan case mismatch")
if (nrow(s) != nrow(x) || nrow(f) != nrow(x) || !setequal(s$candidate_id, x$best_candidate_id) || !setequal(f$candidate_id, x$best_candidate_id)) stop("Current A10/C1/Stan candidate sets differ")
fi <- match(x$best_candidate_id, f$candidate_id)
if (!identical(x$gene_symbol, f$gene[fi])) stop("Current A10/C1 gene mismatch")
if (!is.numeric(s$probability) || any(!is.finite(s$probability)) || any(s$probability < 0 | s$probability > 1) || abs(sum(s$probability) - 1) > 1e-8) stop("Invalid current Stan probability vector")

# Membership is a full union. Historical C4 weights remain historical metadata.
a10_keys <- data.table(gene = x$gene_symbol, variant = x$best_candidate_id)
old_keys <- unique(old[, .(gene, variant)])
reg <- unique(rbindlist(list(a10_keys, old_keys)))
reg[, hypothesis_id := paste0("hyp_", vapply(key(gene, variant), digest::digest, character(1), algo = "sha256", serialize = FALSE))]
reg[, in_A10 := key(gene, variant) %in% key(a10_keys$gene, a10_keys$variant)]
reg[, in_historical_C4 := key(gene, variant) %in% key(old_keys$gene, old_keys$variant)]
xi <- match(key(reg$gene, reg$variant), key(a10_keys$gene, a10_keys$variant))
set(reg, j = "a10_best_integrated_score", value = x$best_integrated_score[xi])
set(reg, j = "a10_best_IMPACT", value = x$best_IMPACT[xi])
si <- match(reg$variant, as.character(s$candidate_id))
si[!reg$in_A10] <- NA_integer_
for (nm in setdiff(names(s), c("case_id", "candidate_id"))) set(reg, j = paste0("stan_", nm), value = s[[nm]][si])
reg[, stan_predictor_id := fifelse(in_A10, pred$predictor_id, NA_character_)]
reg[, model_status := fifelse(in_A10, "current_Stan_scored", "historical_only_not_scored_by_current_Stan")]
reg[, historical_C4_status := fifelse(in_historical_C4, "historical_source_not_current_model_evidence", "absent")]
reg[, candidate_retained := TRUE]
if (!all(key(a10_keys$gene, a10_keys$variant) %in% key(reg$gene, reg$variant)) || !all(key(old_keys$gene, old_keys$variant) %in% key(reg$gene, reg$variant))) stop("Union lost input hypotheses")

# Link only A2 evidence, not all A1 rows. Multiple variants do not duplicate q_j.
links <- copy(a2)
links[, gene := gene_symbol]
links <- links[gene %in% reg$gene]
links[, current_A2_evidence := TRUE]
links[, representative_in_union := key(gene, candidate_id) %in% key(reg$gene, reg$variant)]
reg[, n_current_A2_variants := vapply(gene, function(g) uniqueN(links[gene == g, candidate_id]), integer(1))]
ai <- match(key(reg$gene, reg$variant), key(a2$gene_symbol, a2$candidate_id))
reg[, representative_in_current_A2 := !is.na(ai)]
reg[, current_representative_qc_pass := truth(a2$pass_core_qc[ai])]
if (any(reg$in_A10 & (!reg$representative_in_current_A2 | !reg$current_representative_qc_pass))) stop("Current A10 representative absent from A2 or not core-QC passing")

# Primary C4 weights come directly from the existing trained Stan predictions.
primary <- copy(reg[in_A10])
primary[, q_j := stan_probability]
primary[, consequence := a2$consequence[match(key(gene, variant), key(a2$gene_symbol, a2$candidate_id))]]
primary[, impact := a2$IMPACT[match(key(gene, variant), key(a2$gene_symbol, a2$candidate_id))]]
primary[, cclass := consequence_class(consequence, impact)]
mp <- as.data.table(mechanism_prior_table(cfg$mechanism_prior, unique(primary$cclass)))
need(mp, c("cclass", "mechanism", "p_m"), "Configured mechanism prior")
mp[, p_m := as.numeric(p_m)]
if (anyDuplicated(mp[, .(cclass, mechanism)]) || any(!is.finite(mp$p_m)) || any(mp$p_m < 0 | mp$p_m > 1)) stop("Invalid configured mechanism prior")
if (any(mp[, abs(sum(p_m) - 1) > 1e-8, by = cclass]$V1) || !all(primary$cclass %in% mp$cclass)) stop("Mechanism priors do not cover classes with normalized rows")
fh <- merge(primary, mp, by = "cclass", all.x = TRUE, allow.cartesian = TRUE, sort = FALSE)
if (anyNA(fh$mechanism) || anyNA(fh$p_m)) stop("Missing mechanism prior")
sc <- copy(fh)
sc[, scenario_id := paste0(hypothesis_id, "__", mechanism)]
sc[, scenario_weight_raw := q_j * p_m]
if (sum(sc$scenario_weight_raw) <= 0) stop("Zero scenario mass")
sc[, scenario_weight := scenario_weight_raw / sum(scenario_weight_raw)]
sc[, weighting_scheme := "current_Stan_CV_mixture_x_configured_mechanism_prior"]
if (anyDuplicated(sc$scenario_id) || abs(sum(sc$scenario_weight) - 1) > 1e-8) stop("Invalid D0 scenario table")
# Zero-weight scenarios remain present; no top-k or probability threshold is applied.
explore <- copy(reg[!in_A10])
explore[, exploration_status := "retained_unweighted_requires_current_scoring_and_mechanism_review"]
explore[, q_j := NA_real_]
explore[, scenario_weight := NA_real_]
old_keep <- copy(old)
old_keep[, historical_only_do_not_use_as_current_Stan := TRUE]
if (!identical(vapply(input_paths, sha, character(1)), input_hashes)) stop("Input changed during execution")
if (!dir.create(a$outdir, recursive = TRUE)) stop("Cannot create new output directory")
paths <- c(registry = "candidate_registry.tsv", variants = "shortlisted_variant_evidence.tsv", primary_weights = "current_Stan_hypothesis_weights.tsv", current_C4 = "current_Stan_final_hypothesis_set.tsv", current_D0 = "current_Stan_scenarios.tsv", exploration = "protected_unscored_candidates.tsv", historical_C4 = "historical_C4_preserved.tsv")
tables <- list(reg, links, primary, fh, sc, explore, old_keep)
for (i in seq_along(paths)) fwrite(tables[[i]], file.path(a$outdir, paths[i]), sep = "\t", na = "")
contract <- c(
  "Candidate registry = A10 UNION historical C4, keyed by gene and representative variant.",
  "All A10 candidates retained; no probability or top-k pruning.",
  "Primary weights are unchanged current Stan CV probabilities, not a new blended posterior.",
  "A10 scores and historical C4 probabilities are separate metadata, not additional independent likelihoods.",
  "Historical-only candidates remain protected but unweighted; this artifact does not assess case/build compatibility of historical-only rows.",
  "No zero imputation or invented Bayesian weights for unscored candidates.",
  "Only matching-run A2 variant evidence is linked; full A1 archive is not an active Track 2 input.",
  "No trans phase, variant pathogenicity, biallelic mechanism or clinical efficacy is inferred.",
  "Primary mechanisms use the existing configured consequence-based prior for the representative, not a new pair-aware mechanism model.",
  "Track 2 must consume the registry AND primary scenario table AND protected-unscored table, retaining coverage status.",
  "Existing D0/Track 2 scripts are not edited. Do not feed the registry directly into the old scalar-weight D0 interface.",
  "New downstream code must not discard protected-unscored candidates via inner joins or zero imputation.",
  "Unweighted exploration must be reported separately from probability-weighted aggregation."
)
writeLines(contract, file.path(a$outdir, "HANDOFF_CONTRACT.txt"))
manifest <- list(version = "A10_C4_union_Stan_v1", timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), sample_label = a$sample,
 inputs = lapply(seq_along(input_paths), function(i) list(role = names(input_paths)[i], path = normalizePath(input_paths[i]), sha256 = unname(input_hashes[i]))),
 outputs = lapply(paths, function(p) list(path = p, sha256 = sha(file.path(a$outdir, p)))),
 counts = list(A10 = nrow(x), historical_C4_hypotheses = nrow(old_keys), union = nrow(reg), primary_current_Stan = nrow(primary), protected_unscored = nrow(explore), shortlisted_evidence_rows = nrow(links), primary_scenarios = nrow(sc)),
 predictor_id = pred$predictor_id, predictor_interpretation = pred$interpretation, provisional_divergences = pred$provisional_divergences,
 mechanism_prior_mode = cfg$mechanism_prior$mode, model_refitted = FALSE, source_files_modified = FALSE, downstream_workflow_rewired = FALSE,
 contract = contract, session_info = capture.output(sessionInfo()))
jsonlite::write_json(manifest, file.path(a$outdir, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("\nUNION HANDOFF COMPLETE\n")
cat("A10 hypotheses:", nrow(x), "\nHistorical C4 hypotheses:", nrow(old_keys), "\nUnion hypotheses:", nrow(reg), "\nCurrent Stan weighted hypotheses:", nrow(primary), "\nProtected unscored hypotheses:", nrow(explore), "\nA2 supporting variant rows:", nrow(links), "\nPrimary scenarios:", nrow(sc), "\nScenario weight sum:", sum(sc$scenario_weight), "\n")
cat("No refit, no score blending, no source edits. Downstream workflow not yet rewired.\n")
cat("Outputs:", a$outdir, "\n")
