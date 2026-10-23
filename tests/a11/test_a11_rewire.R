args <- commandArgs(trailingOnly=TRUE)
repo <- args[1]; legacy_a <- args[2]; legacy_f <- args[3]; td <- args[4]
suppressPackageStartupMessages(library(data.table))
run <- function(script, arguments=character(), success=TRUE) {
  log <- tempfile(tmpdir=td)
  status <- system2(Sys.which("Rscript"), c("--vanilla",shQuote(script),shQuote(arguments)), stdout=log,stderr=log)
  if ((status==0)!=success) stop(paste(readLines(log,warn=FALSE),collapse="\n"))
}
port_a <- file.path(repo,"workflow/layer_a/scripts/a11_dependence_sensitivity.R")
port_f <- file.path(repo,"workflow/layer_f/scripts/f_a11_coverage_report.R")
port_c <- file.path(repo,"workflow/layer_c/scripts/attach_a11_multimodel.R")
dir.create(file.path(td,"input folder"))
evidence <- file.path(td,"input folder","synthetic.csv")
n <- 34L
d <- data.table(candidate_id=paste0("1:",seq_len(n),":A:C"),gene_symbol=paste0("SYNTH_G",seq_len(n)),
 genomic=seq(.15,.95,length.out=n),hpo=((seq_len(n)*7)%%31+1)/33,
 string=((seq_len(n)*11)%%29+1)/31,multiplex=((seq_len(n)*13)%%31+1)/33)
d[33,`:=`(hpo=NA_real_,string=NA_real_)]
d[34,`:=`(hpo=NA_real_,string=NA_real_,multiplex=NA_real_)]
for (nm in c("genomic","hpo","string","multiplex")) d[,(paste0("has_",nm)):=is.finite(get(nm))]
fwrite(d,evidence)
quote_r <- function(x) encodeString(x,quote='"')
old_a <- paste(readLines(legacy_a,warn=FALSE),collapse="\n")
old_a <- sub('INPUT <- "results/private/tier01_candidate_evidence.csv"',paste0("INPUT <- ",quote_r(evidence)),old_a,fixed=TRUE)
old_a <- sub('OUT_DIR <- "bayesian_dependence"',paste0("OUT_DIR <- ",quote_r(file.path(td,"legacy_run"))),old_a,fixed=TRUE)
old_a <- sub("N_DRAWS <- 10000L","N_DRAWS <- 64L",old_a,fixed=TRUE)
old_script <- file.path(td,"legacy_calculation.R"); writeLines(old_a,old_script)
run(old_script)
newrun <- file.path(td,"new_run")
run(port_a,c("--input",evidence,"--outdir",newrun,"--draws","64"))
for (name in list.files(file.path(td,"legacy_run"),pattern="[.]csv$")) {
  if (name == "run_manifest.csv" || grepl("^[^/]+_bayesian_dependence_sensitivity_ensemble[.]csv$", name)) next
  before <- fread(file.path(td,"legacy_run",name)); after <- fread(file.path(newrun,name))
  if (!isTRUE(all.equal(before,after,tolerance=1e-12))) stop("Legacy calculation differs: ",name)
}
repeat_run <- file.path(td,"repeat_run")
run(port_a,c("--input",evidence,"--outdir",repeat_run,"--draws","64"))
cons_name <- "bayesian_dependence_sensitivity_ensemble_consensus.csv"
if (!isTRUE(all.equal(fread(file.path(newrun,cons_name)),fread(file.path(repeat_run,cons_name)),tolerance=0))) stop("Repeatability failed")
ledger <- fread(file.path(newrun,"candidate_coverage.csv"))
if (nrow(ledger)!=34L || sum(ledger$entered_analysis)!=33L) stop("Coverage ledger failed")
old_f <- paste(readLines(legacy_f,warn=FALSE),collapse="\n")
for (item in list(c('INPUT_FILE <- "results/private/tier01_candidate_evidence.csv"',paste0("INPUT_FILE <- ",quote_r(evidence))),
 c('CONSENSUS_FILE <- "bayesian_dependence/bayesian_dependence_sensitivity_ensemble_consensus.csv"',paste0("CONSENSUS_FILE <- ",quote_r(file.path(newrun,cons_name)))),
 c('COPULA_FILE <- "bayesian_dependence/copula_dependence_summary.csv"',paste0("COPULA_FILE <- ",quote_r(file.path(newrun,"copula_dependence_summary.csv")))),
 c('OUT_DIR <- "results/private/bayesian_dependence_report"',paste0("OUT_DIR <- ",quote_r(file.path(td,"legacy_report")))))) old_f <- sub(item[1],item[2],old_f,fixed=TRUE)
old_report_script <- file.path(td,"legacy_report.R"); writeLines(old_f,old_report_script); run(old_report_script)
report <- file.path(td,"new_report")
report_args <- c("--input",evidence,"--consensus",file.path(newrun,cons_name),"--copula",file.path(newrun,"copula_dependence_summary.csv"),"--a11-manifest",file.path(newrun,"run_provenance.json"),"--outdir",report)
run(port_f,report_args)
for (name in c("bayesian_dependence_ranked_candidates_with_coverage.csv","bayesian_dependence_top25_evidence.tsv","bayesian_dependence_top10_evidence.tsv","copula_dependence_summary.tsv")) if (!isTRUE(all.equal(fread(file.path(td,"legacy_report",name)),fread(file.path(report,name)),tolerance=1e-12))) stop("Report differs: ",name)
changed <- copy(d); changed[1,genomic:=.2]; changed_path <- file.path(td,"changed.csv"); fwrite(changed,changed_path)
bad <- report_args; bad[match("--input",bad)+1L] <- changed_path; bad[match("--outdir",bad)+1L] <- file.path(td,"bad_report")
run(port_f,bad,success=FALSE)
run(port_a,c("--input",evidence,"--outdir",newrun,"--draws","64"),success=FALSE)
union <- file.path(td,"synthetic_union"); dir.create(union)
reg <- data.table(hypothesis_id=paste0("SYNTH_H",seq_len(35)),gene=c(d$gene_symbol,"SYNTH_HIST"),variant=c(d$candidate_id,"2:1:G:T"),in_A10=c(rep(TRUE,34),FALSE),stan_probability=c(rep(1/34,34),NA_real_))
stan <- reg[in_A10==TRUE,.(hypothesis_id,gene,variant,q_j=stan_probability)]
stan[,`:=`(mechanism="synthetic",p_m=1,scenario_weight=q_j,scenario_id=paste0(hypothesis_id,"__synthetic"))]
files <- c("candidate_registry.tsv","current_Stan_scenarios.tsv","shortlisted_variant_evidence.tsv","protected_unscored_candidates.tsv")
tabs <- list(reg,stan,d[,.(candidate_id,gene_symbol)],reg[in_A10==FALSE])
for (i in seq_along(files)) fwrite(tabs[[i]],file.path(union,files[i]),sep="\t")
jsonlite::write_json(list(current_build="GRCh38",legacy_build="GRCh38",predictor_id="SYNTHETIC_ONLY",mechanism_prior_mode="synthetic",outputs=lapply(files,function(n) list(path=n,sha256=digest::digest(file=file.path(union,n),algo="sha256")))),file.path(union,"manifest.json"),auto_unbox=TRUE)
attachment <- file.path(td,"new_attachment")
run(port_c,c("--union",union,"--consensus",file.path(newrun,cons_name),"--a11_input",evidence,"--a11_manifest",file.path(newrun,"run_provenance.json"),"--outdir",attachment))
r <- fread(file.path(attachment,"multimodel_candidate_registry.tsv"))
if (nrow(r)!=35L || sum(r$a11_scored)!=33L || !is.na(r[gene=="SYNTH_G34",a11_p_rank1]) || !is.na(r[gene=="SYNTH_HIST",a11_p_rank1])) stop("Attachment coverage failed")
w <- fread(file.path(attachment,"source_weights_long.tsv"))
if (any(!w$source_scored & !is.na(w$weight_value)) || any(abs(w[,sum(weight_value,na.rm=TRUE),by=weight_source]$V1-1)>1e-6)) stop("Source separation failed")
cat("PASS: legacy calculation regression; fixed-seed repeatability; eligibility ledger; legacy report regression; input-mismatch rejection; overwrite rejection; union coverage; source separation\n")
