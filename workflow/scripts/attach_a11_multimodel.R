#!/usr/bin/env Rscript
# Attach existing A11 results to a validated union; never refit or blend models.
a <- commandArgs(trailingOnly = TRUE)
if (!length(a) || length(a) %% 2L) stop("Use --union DIR --consensus CSV --a11_input CSV --outdir NEW_DIR")
keys <- a[seq(1L,length(a),2L)]
allowed <- c("--union","--consensus","--a11_input","--outdir")
if (anyDuplicated(keys) || !setequal(keys,allowed)) {
  stop("Supply exactly --union, --consensus, --a11_input and --outdir, once each.")
}
o <- as.list(a[seq(2L,length(a),2L)])
names(o) <- sub("^--","",keys)
req <- function(k) { if (is.null(o[[k]]) || !nzchar(o[[k]])) stop("Missing --",k); o[[k]] }
for (p in c("data.table","digest","jsonlite")) if (!requireNamespace(p,quietly=TRUE)) stop("Required package: ",p)
suppressPackageStartupMessages(library(data.table))
u <- req("union"); dest <- req("outdir")
if (file.exists(dest)) stop("Refusing to reuse output directory")
sha <- function(p) digest::digest(file=p,algo="sha256")
need <- function(x,n,label) if (!all(n %in% names(x))) stop(label," missing: ",paste(setdiff(n,names(x)),collapse=", "))
read_tab <- function(p) {
  if (!file.exists(p)) stop("Missing file: ",p)
  h <- fread(p,nrows=0L)
  ch <- intersect(c("gene","variant","hypothesis_id","scenario_id","candidate_id","gene_symbol"),names(h))
  if (length(ch)) fread(p,colClasses=list(character=ch)) else fread(p)
}
canon <- function(x) {
  z <- trimws(as.character(x)); z <- sub("^GRCh38:","",z,ignore.case=TRUE)
  if (anyNA(z) || any(!nzchar(z)) || any(grepl("^(GRCh|hg)[0-9]",z,ignore.case=TRUE))) stop("Missing ID or unexpected assembly prefix")
  z <- sub("^chr","",z,ignore.case=TRUE)
  sub("^M:", "MT:", z, ignore.case=TRUE)
}
key <- function(g,v) paste(g,v,sep="\r")
flag <- function(x) toupper(as.character(x)) %in% c("TRUE","T","1")
manifest_path <- file.path(u,"manifest.json")
base <- jsonlite::read_json(manifest_path,simplifyVector=FALSE)
if (!identical(base$current_build,"GRCh38") || !identical(base$legacy_build,"GRCh38")) stop("Require GRCh38-declared union v1.1 or later")
files <- c("candidate_registry.tsv","current_Stan_scenarios.tsv","shortlisted_variant_evidence.tsv","protected_unscored_candidates.tsv")
for (n in files) {
  entries <- base$outputs
  hit <- which(vapply(entries,function(e) identical(as.character(e$path),n),logical(1)))
  if (length(hit)!=1L || !identical(as.character(entries[[hit]]$sha256),sha(file.path(u,n)))) stop("Union manifest checksum mismatch: ",n)
}
inputs <- c(union_manifest=manifest_path,registry=file.path(u,files[1]),stan_scenarios=file.path(u,files[2]),
            variants=file.path(u,files[3]),protected=file.path(u,files[4]),consensus=req("consensus"),a11_input=req("a11_input"))
mp <- file.path(dirname(req("consensus")),"run_manifest.csv")
if (file.exists(mp)) inputs <- c(inputs,a11_run_manifest=mp)
input_sha <- vapply(inputs,sha,character(1))
reg <- read_tab(inputs[["registry"]]); stan <- read_tab(inputs[["stan_scenarios"]]); cons <- read_tab(inputs[["consensus"]]); inp <- read_tab(inputs[["a11_input"]])
need(reg,c("hypothesis_id","gene","variant","in_A10","stan_probability"),"Union registry")
need(stan,c("hypothesis_id","scenario_id","gene","variant","q_j","mechanism","p_m","scenario_weight"),"Stan scenarios")
need(cons,c("candidate_id","gene_symbol","posterior_score_mean","p_rank1","p_top3","p_top5","p_top10"),"A11 consensus")
need(inp,c("candidate_id","gene_symbol","genomic","hpo","string","multiplex"),"A11 input")
assert_ids <- function(x,cols,label) {
  for (nm in cols) {
    v <- as.character(x[[nm]])
    if (anyNA(v) || any(!nzchar(trimws(v)))) {
      stop(label, " has missing or blank ", nm)
    }
  }
}
assert_ids(reg,c("hypothesis_id","gene","variant"),"Registry")
assert_ids(stan,c("hypothesis_id","scenario_id","gene","variant","mechanism"),"Stan scenarios")
assert_ids(cons,c("candidate_id","gene_symbol"),"Consensus")
assert_ids(inp,c("candidate_id","gene_symbol"),"A11 input")
if (!nrow(reg) || !nrow(stan) || !nrow(cons) || !nrow(inp)) {
  stop("Registry, scenarios, consensus and A11 input must be nonempty.")
}
if (anyDuplicated(reg$hypothesis_id)) stop("Duplicate registry hypotheses")
if (anyDuplicated(stan$scenario_id)) stop("Duplicate Stan scenario IDs")
rk <- key(reg$gene,canon(reg$variant)); ck <- key(cons$gene_symbol,canon(cons$candidate_id)); ik <- key(inp$gene_symbol,canon(inp$candidate_id))
if (anyDuplicated(rk) || anyDuplicated(ck) || anyDuplicated(ik)) stop("Duplicate gene/allele keys")
membership <- toupper(trimws(as.character(reg$in_A10)))
if (anyNA(membership) ||
    any(!membership %in% c("TRUE","T","1","FALSE","F","0"))) {
  stop("Invalid in_A10 membership flag.")
}
reg[, in_A10 := membership %in% c("TRUE","T","1")]
current <- reg$in_A10
if (!setequal(ik,rk[current])) stop("A11 input does not match complete current A10 registry")
if (!all(ck %in% rk[current])) stop("Consensus has candidates outside current A10")
for (nm in c("p_rank1","p_top3","p_top5","p_top10")) {
  v <- cons[[nm]]
  if (!is.numeric(v) || any(!is.finite(v)) || any(v<0 | v>1)) stop("Invalid A11 probability field: ",nm)
}
if (any(cons$p_rank1 > cons$p_top3 + 1e-12 |
        cons$p_top3 > cons$p_top5 + 1e-12 |
        cons$p_top5 > cons$p_top10 + 1e-12)) {
  stop("A11 rank frequencies must satisfy rank1 <= top3 <= top5 <= top10.")
}
if (abs(sum(cons$p_rank1)-1)>1e-6) stop("A11 p_rank1 does not sum to 1; do not silently renormalize")
if (!is.numeric(reg$stan_probability) ||
    any(!is.finite(reg$stan_probability[current])) ||
    any(reg$stan_probability[current] < 0 | reg$stan_probability[current] > 1) ||
    abs(sum(reg$stan_probability[current])-1)>1e-8) {
  stop("Invalid current Stan weights")
}
if (!all(unique(stan$hypothesis_id) %in% reg$hypothesis_id[current]) || !setequal(unique(stan$hypothesis_id),reg$hypothesis_id[current])) stop("Stan scenario/registry mismatch")
if (anyDuplicated(stan[,.(hypothesis_id,mechanism)])) stop("Duplicate mechanism scenarios")
if (any(!is.finite(stan$p_m)) || any(stan$p_m<0 | stan$p_m>1) || any(stan[,abs(sum(p_m)-1)>1e-8,by=hypothesis_id]$V1)) stop("Invalid mechanism priors")
for (nm in c("q_j","p_m","scenario_weight")) {
  v <- stan[[nm]]
  if (!is.numeric(v) || any(!is.finite(v)) || any(v < 0 | v > 1)) {
    stop("Invalid Stan scenario field: ",nm)
  }
}
si <- match(stan$hypothesis_id,reg$hypothesis_id)
if (any(abs(stan$q_j-reg$stan_probability[si])>1e-12)) {
  stop("Stan q_j differs from saved registry")
}
if (any(key(stan$gene,canon(stan$variant)) != rk[si])) {
  stop("Stan scenario gene/allele differs from registry hypothesis.")
}
expected_weight <- stan$q_j * stan$p_m
if (abs(sum(expected_weight)-1)>1e-8 ||
    abs(sum(stan$scenario_weight)-1)>1e-8 ||
    any(abs(stan$scenario_weight-expected_weight)>1e-8)) {
  stop("Saved Stan scenario weights do not match q_j * p_m.")
}
ci <- match(rk,ck)
for (nm in names(cons)) set(reg,j=paste0("a11_",nm),value=cons[[nm]][ci])
reg[, a11_scored := !is.na(ci)]
reg[, a11_coverage_status := fifelse(a11_scored,"present_in_existing_A11_consensus",fifelse(in_A10==TRUE,"entered_A11_but_absent_from_consensus","outside_current_A11_input"))]

# Every registry hypothesis has a row for each source; unavailable is NA, not 0.
stan_long <- reg[,.(hypothesis_id,gene,variant,weight_source="external_Stan_CV",weight_value=stan_probability,
                   source_scored=in_A10==TRUE,weight_semantics="within_candidate_set_Stan_predictive_mixture",source_universe_size=sum(current))]
a11_long <- reg[,.(hypothesis_id,gene,variant,weight_source="A11_dependence_rank_first",weight_value=a11_p_rank1,
                  source_scored=a11_scored,weight_semantics="rank_first_probability_sensitivity_NOT_causal_probability",source_universe_size=nrow(cons))]
weights <- rbindlist(list(stan_long,a11_long))
if (any(!weights$source_scored & !is.na(weights$weight_value))) stop("Unscored source received a weight")

# Preserve the Stan branch. Construct a separately labelled rank-first sensitivity branch.
set(stan,j="analysis_branch",value=rep("external_Stan_CV",nrow(stan)))
mi <- match(stan$hypothesis_id,reg$hypothesis_id)
a11 <- copy(stan[reg$a11_scored[mi]==TRUE])
aj <- match(a11$hypothesis_id,reg$hypothesis_id)
set(a11,j="stan_q_j_reference",value=a11$q_j)
set(a11,j="q_j",value=reg$a11_p_rank1[aj])
set(a11,j="scenario_weight_raw",value=a11$q_j*a11$p_m)
total <- sum(a11$scenario_weight_raw)
if (abs(total-1)>1e-6) stop("A11 sensitivity scenario mass differs from rank-probability mass")
set(a11,j="scenario_weight",value=a11$scenario_weight_raw/total)
set(a11,j="analysis_branch",value=rep("A11_rank_first_sensitivity_NOT_causal_posterior",nrow(a11)))
set(a11,j="weighting_scheme",value=rep("existing_A11_p_rank1_x_same_configured_mechanism_prior",nrow(a11)))
set(a11,j="scenario_id",value=paste0("A11rank__",a11$scenario_id))
if (anyDuplicated(a11$scenario_id) || abs(sum(a11$scenario_weight)-1)>1e-8) stop("Invalid A11 scenario table")
# The other stan_* fields remain explicitly Stan references, not A11 estimates.

coverage <- data.table(candidate_id=inp$candidate_id,gene_symbol=inp$gene_symbol,in_A11_consensus=ik %in% ck)
for (nm in c("genomic","hpo","string","multiplex")) {
  v <- suppressWarnings(as.numeric(inp[[nm]]))
  set(coverage,j=paste0(nm,"_finite"),value=is.finite(v))
}
coverage[, all_four_channels_finite := genomic_finite & hpo_finite & string_finite & multiplex_finite]
missing <- copy(coverage[in_A11_consensus==FALSE])
if (sum(reg$a11_scored)!=nrow(cons)) stop("Consensus linkage lost rows")
if (!identical(vapply(inputs,sha,character(1)),input_sha)) stop("Input changed during execution")
if (!dir.create(dest,recursive=TRUE)) stop("Could not create new output directory")
paths <- c(registry="multimodel_candidate_registry.tsv",weights="source_weights_long.tsv",stan="external_Stan_scenarios.tsv",a11="A11_rank_first_sensitivity_scenarios.tsv",coverage="A11_input_coverage.tsv",missing="A10_absent_from_A11_consensus.tsv")
tables <- list(reg,weights,stan,a11,coverage,missing)
for (i in seq_along(paths)) fwrite(tables[[i]],file.path(dest,paths[i]),sep="\t",na="")
contract <- c(
 "Union membership unchanged. A11 results attached by gene and assembly-reconciled allele key.",
 "Current A10 rows absent from A11 remain retained, with A11 probabilities unavailable (NA), not zero.",
 "Posterior score means are metadata, not normalized causal probabilities.",
 "A11 p_rank1 defines a separately labelled rank-first sensitivity branch, not a causal posterior or externally validated model.",
 "Rank-first zero for a scored candidate is distinct from an unavailable result for an unscored candidate.",
 "The two branches are not blended, multiplied, or claimed to be the same statistical estimand.",
 "Each branch is conditional on its own candidate universe; ranks/probabilities are not directly calibrated across universes.",
 "The same existing configured mechanism priors are reused as sensitivity assumptions, not learned mechanism probabilities.",
 "Finite-channel audit is diagnostic only; it does not prove why A11 excluded a row.",
 "No model fit or clinical reclassification is performed. Historical-only coverage and unresolved rows remain in the source union bundle.",
 "The existing Track 2 scripts and workflow are not edited; they must preserve branch IDs and union coverage before deployment.",
 "Source presence/sum checks are integrity checks, not statistical validation."
)
writeLines(contract,file.path(dest,"HANDOFF_CONTRACT.txt"))
manifest <- list(version="A11_multimodel_handoff_v1_2_validation",time_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
 inputs=lapply(seq_along(inputs),function(i) list(role=names(inputs)[i],path=normalizePath(inputs[i]),sha256=unname(input_sha[i]))),
 outputs=lapply(paths,function(n) list(path=n,sha256=sha(file.path(dest,n)))),
 counts=list(union=nrow(reg),current_A10=sum(current),A11_scored=nrow(cons),current_A10_missing_A11=sum(current & !reg$a11_scored),Stan_scenarios=nrow(stan),A11_rank_sensitivity_scenarios=nrow(a11)),
 original_Stan_predictor_id=base$predictor_id,mechanism_prior_mode=base$mechanism_prior_mode,
 model_refitted=FALSE,model_sources_blended=FALSE,downstream_workflow_rewired=FALSE,contract=contract,session_info=capture.output(sessionInfo()))
jsonlite::write_json(manifest,file.path(dest,"manifest.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat("\nMULTIMODEL HANDOFF COMPLETE\n")
cat("Union hypotheses preserved:",nrow(reg),"\nCurrent A10:",sum(current),"\nA11 scored:",nrow(cons),"\nCurrent A10 unavailable in A11:",nrow(missing),"\nStan scenarios:",nrow(stan),"\nA11 rank-first sensitivity scenarios:",nrow(a11),"\n")
cat("\nAmong candidates absent from A11, nonfinite input counts:\n")
for (nm in c("genomic","hpo","string","multiplex")) cat(nm,":",sum(!missing[[paste0(nm,"_finite")]]),"\n")
cat("\nThis diagnoses coverage; exclusion reasons require inspecting A11 code.\n")
cat("Separate weight sums: Stan =",sum(stan_long$weight_value,na.rm=TRUE),"; A11 rank-first =",sum(a11_long$weight_value,na.rm=TRUE),"\n")
cat("Outputs:",dest,"\nNo refit, no source blending, downstream not yet rewired.\n")
