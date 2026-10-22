#!/usr/bin/env Rscript
# Separate Stan CV-ensemble path; never calls sample() or freeze_model().
a <- commandArgs(trailingOnly=TRUE)
if (!length(a) || length(a)%%2) stop("Use --mode build|score|calibrate and --key value pairs")
opt <- as.list(a[seq(2,length(a),2)]); names(opt) <- sub("^--","",a[seq(1,length(a),2)])
get <- function(k,default=NULL) {v<-opt[[k]]; if(is.null(v)) default else v}
req <- function(k) {v<-get(k); if(is.null(v)||!nzchar(v)) stop("Required --",k); v}
yes <- function(k) identical(tolower(get(k,"no")),"yes")
mode <- req("mode")
for (pkg in c("data.table","posterior","cmdstanr","digest")) if(!requireNamespace(pkg,quietly=TRUE)) stop("Required package: ",pkg)
suppressPackageStartupMessages(library(data.table))
runner <- get("runner","written_final_scripts/workflow/lib/stan_runner.R")
helpers <- get("helpers","scripts/Final Scripts/b14_stan_cv_helpers.R")
if(!file.exists(runner)||!file.exists(helpers)) stop("Runner or companion helper not found")
source(runner); source(helpers)
stanpath <- file.path(get("stan-dir","written_final_scripts/workflow/stan"),"m3b.stan")
if(!file.exists(stanpath)) stop("Missing m3b.stan")
sha <- function(p) digest::digest(file=p,algo="sha256")
write_new <- function(x,path) {
  if(file.exists(path)||file.exists(paste0(path,".sha256"))) stop("Refusing to overwrite: ",path)
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  saveRDS(x,path,compress="gzip")
  if(!file.exists(path)||file.info(path)$size<=0) stop("Write failed")
  writeLines(sha(path),paste0(path,".sha256"))
  cat("Saved:",path,"\n")
}
read_checked <- function(path) {
  h<-paste0(path,".sha256")
  if(!file.exists(h)||!identical(trimws(readLines(h,n=1)),sha(path))) stop("Missing/mismatched checksum: ",path)
  readRDS(path)
}
prepare_raw <- function(d) {
  need<-c("case_id","candidate_id","G","H","N")
  if(!all(need%in%names(d))) stop("Missing raw input columns: ",paste(setdiff(need,names(d)),collapse=", "))
  d<-copy(d); d[,`:=`(case_id=as.character(case_id),candidate_id=as.character(candidate_id))]
  if(anyDuplicated(d[,.(case_id,candidate_id)])) stop("Duplicate candidates in a case")
  if("group"%in%names(d)&&any(is.na(d$group)|d$group!="all")) stop("Only original fitted group all is supported")
  d[,`:=`(group="all",hpo_score=H,network_score=N)]
  d
}
read_input <- function(path) {
  z<-fread(path,colClasses=list(character=c("case_id","candidate_id")))
  prepare_raw(z)
}
model <- NULL
compiled <- function() {
  if(is.null(model)) model <<- cmdstanr::cmdstan_model(stanpath)
  model
}

if(mode=="build") {
  if(!yes("confirm-raw-features")) stop("Supply --confirm-raw-features yes only for the ORIGINAL B11 --features input")
  run<-req("run"); features<-req("features"); protocol<-req("candidate-protocol")
  if(!nzchar(protocol)) stop("Candidate protocol ID required")
  d<-read_input(features)
  if(!"split"%in%names(d)||!"label"%in%names(d)) stop("Original split and label columns required")
  d<-copy(d[split=="train"])
  if(!nrow(d)) stop("No train rows")
  gs<-unique(d[,.(duplicate_group=case_id)])
  folds<-rep_len(seq_len(3),nrow(gs))
  d[,inner_fold:=folds[match(case_id,gs$duplicate_group)]]
  dg<-fread(file.path(run,"stan_oof","oof_diagnostics.tsv"))
  components<-list(); checks<-list()
  verify_n<-as.integer(get("verify-cases","3"))
  if(!is.finite(verify_n)||verify_n<1) stop("verify-cases must be positive")
  for(fold in 1:3) {
    path<-file.path(run,"stan_oof","chains",sprintf("m3b_fold%02d.rds",fold))
    cat("Packaging and checking",basename(path),"\n"); flush.console()
    x<-readRDS(path)
    if(!identical(x$model_id,"m3b")||x$inner_fold!=fold) stop("Wrong chain artifact")
    ar<-posterior::as_draws_array(x$param_draws)
    if(posterior::nchains(ar)!=4) stop("Expected original four chains")
    s<-as.data.frame(posterior::summarise_draws(ar,"rhat","ess_bulk","ess_tail"))
    if(!all(is.finite(s$rhat)&s$rhat<1.01&is.finite(s$ess_bulk)&s$ess_bulk>400&is.finite(s$ess_tail)&s$ess_tail>400)) stop("Parameter gate failed")
    row<-dg[model_id=="m3b" & inner_fold==fold]
    if(nrow(row)!=1||!is.finite(row$min_ebfmi)||row$min_ebfmi<=.3) stop("Invalid sampler summary")
    if(row$divergences>0&&!yes("allow-divergences")) stop("Retained fit has divergences; explicit --allow-divergences yes required")
    if(!yes("allow-unverified-treedepth")) stop("Treedepth counts are absent; explicit --allow-unverified-treedepth yes required")
    tr<-copy(d[inner_fold!=fold]); va<-copy(d[inner_fold==fold])
    st<-fit_prep_state(tr,c("G","H","N"),"hpo_score","network_score")
    setorder(va,case_id,candidate_id)
    if(ncol(x$p_new_draws)!=nrow(va)) stop("Held-out candidate count differs from original")
    fr<-list(model_id="m3b",draws=ar,prep_state=st,group_levels="all",feature_cols=c("G","H","N"),
      prior=default_priors(3),development_case_ids=unique(d$case_id),candidate_protocol_id=protocol,
      stan_sha=sha(stanpath),provisional_divergences=row$divergences>0,
      gate=list(parameter_diagnostics=s,sampler_summary=as.data.frame(row),treedepth_verified=FALSE),
      retained_attempt=x$attempt,source_chain_sha=sha(path),draw_policy="all_retained_no_rounding")
    fr$predictor_id<-digest::digest(fr,algo="sha256")
    cases<-head(unique(va$case_id),verify_n)
    ii<-which(va$case_id%in%cases)
    scored<-score_stan_full(fr,compiled(),va[ii],stanpath)
    cols<-match(sprintf("p_new[%d]",ii),colnames(x$p_new_draws))
    if(anyNA(cols)) stop("Cannot align original p_new columns")
    expected<-x$p_new_draws[,cols,drop=FALSE]
    actual<-scored$p_new_draws
    if(!identical(dim(expected),dim(actual))) stop("Verification draw dimensions differ")
    err<-max(abs(expected-actual))
    if(!is.finite(err)||err>1e-6) stop("Prediction reproduction failed; max absolute error=",err,
      ". Check original input, row ordering, preprocessing and source version.")
    checks[[fold]]<-data.frame(fold=fold,verified_cases=length(cases),verified_candidates=length(ii),max_abs_draw_error=err)
    components[[fold]]<-fr
    rm(x,ar,tr,va,expected,actual,scored); invisible(gc())
  }
  bundle<-list(kind="equal_weight_three_fold_CV_ensemble_NOT_full_data_posterior",model_id="m3b",
    fold_weights=rep(1/3,3),components=components,development_case_ids=unique(d$case_id),
    candidate_protocol_id=protocol,training_features_sha=sha(features),runner_sha=sha(runner),
    stan_sha=sha(stanpath),verification=do.call(rbind,checks),
    provisional_divergences=TRUE,treedepth_verified=FALSE,
    note="Reconstructed original B11 split and preprocessing; reproduction verified on specified cases, not all rows")
  bundle$predictor_id<-digest::digest(bundle,algo="sha256")
  print(bundle$verification)
  write_new(bundle,req("out"))
  cat("Built CV ensemble without sampling. Calibration is NOT included.\n")
} else if(mode%in%c("score","calibrate")) {
  if(!yes("confirm-raw-features")) stop("Supply --confirm-raw-features yes only for unstandardized G,H,N using B11 definitions")
  bundle<-read_checked(req("bundle"))
  if(!identical(bundle$stan_sha,sha(stanpath))) stop("Stan source hash mismatch")
  if(!identical(req("candidate-protocol"),bundle$candidate_protocol_id)) stop("Candidate protocol ID mismatch")
  d<-read_input(req("features"))
  if(!nrow(d)) stop("No prediction rows")
  if(any(d$case_id%in%bundle$development_case_ids)) stop("Input cases overlap B11 training/development cases")
  maxcases<-as.integer(get("batch-cases","25"))
  if(!is.finite(maxcases)||maxcases<1) stop("Invalid batch-cases")
  ids<-unique(d$case_id); batches<-split(ids,ceiling(seq_along(ids)/maxcases))
  results<-list(); draw_batches<-list()
  wquant<-function(v,w,p) {o<-order(v); v[o][which(cumsum(w[o])>=p)[1]]}
  for(b in seq_along(batches)) {
    cat("Scoring batch",b,"of",length(batches),"\n"); flush.console()
    z<-copy(d[case_id%in%batches[[b]]]); setorder(z,case_id,candidate_id)
    parts<-lapply(bundle$components,function(fr) score_stan_full(fr,compiled(),z,stanpath))
    if(!all(vapply(parts,function(x) identical(x$summary[,.(case_id,candidate_id)],parts[[1]]$summary[,.(case_id,candidate_id)]),logical(1)))) stop("Component row mismatch")
    draws<-do.call(rbind,lapply(parts,function(x) as.matrix(x$p_new_draws)))
    dw<-unlist(lapply(seq_along(parts),function(j) rep(bundle$fold_weights[j]/nrow(parts[[j]]$p_new_draws),nrow(parts[[j]]$p_new_draws))))
    summary<-copy(parts[[1]]$summary)
    summary[,probability:=drop(crossprod(dw,draws))]
    summary[,probability_lower_95:=vapply(seq_len(ncol(draws)),function(j) wquant(draws[,j],dw,.025),numeric(1))]
    summary[,probability_upper_95:=vapply(seq_len(ncol(draws)),function(j) wquant(draws[,j],dw,.975),numeric(1))]
    summary[,mean_rank:=Reduce(`+`,lapply(seq_along(parts),function(j) bundle$fold_weights[j]*parts[[j]]$summary$mean_rank))]
    summary[,posterior_rank_first:=Reduce(`+`,lapply(seq_along(parts),function(j) bundle$fold_weights[j]*parts[[j]]$summary$posterior_rank_first))]
    summary[,rank_mean_probability:=frank(-probability,ties.method="min"),by=case_id]
    results[[b]]<-summary
    if(mode=="score") draw_batches[[b]]<-list(candidate_map=summary[,.(case_id,candidate_id)],component_draws=lapply(parts,`[[`,"p_new_draws"),component_weights=bundle$fold_weights)
    rm(parts,draws,dw); invisible(gc())
  }
  scored<-list(summary=rbindlist(results),draw_batches=draw_batches,predictor_id=bundle$predictor_id,
    candidate_protocol_id=bundle$candidate_protocol_id,provisional_divergences=TRUE,
    interpretation="CV-ensemble predictive mixture; not a full-data Bayesian posterior")
  if(mode=="calibrate") {
    if(!yes("confirm-untouched")) stop("Require --confirm-untouched yes after confirming labels were not used for selection, fitting, or tuning")
    if(!yes("confirm-comparable")) stop("Require --confirm-comparable yes only after auditing candidate construction and exchangeability")
    if("split"%in%names(d)&&any(d$split=="locked_test")) stop("Refusing to use locked_test as calibration")
    alpha<-as.numeric(strsplit(get("alpha","0.1,0.05"),",",fixed=TRUE)[[1]])
    calibration<-fit_split_conformal(scored,bundle,alpha)
    calibration$provisional_divergences<-TRUE
    calibration$features_sha<-sha(req("features"))
    print(calibration$thresholds)
    write_new(calibration,req("out"))
  } else {
    if(!is.null(get("calibration"))) scored$summary<-apply_split_conformal(scored,read_checked(get("calibration")))
    write_new(scored,req("out"))
    tab<-sub("[.]rds$","",req("out"))
    tab<-paste0(tab,"_summary.tsv")
    if(file.exists(tab)) stop("Summary output exists; RDS already saved; do not overwrite")
    fwrite(scored$summary,tab,sep="\t")
    cat("Summary:",tab,"\n")
    print(head(scored$summary[order(case_id,-probability)],20))
  }
} else stop("Unknown mode")
