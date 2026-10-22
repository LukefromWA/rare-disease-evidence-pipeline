# Source after stan_runner.R. No model fitting calls.
score_stan_full <- function(fr,model,new_dt,stan_file,max_draw_cells=2e7,seed=123L) {
  need<-c("draws","prep_state","group_levels","feature_cols","prior","predictor_id","candidate_protocol_id")
  if(!all(need%in%names(fr))) stop("Incomplete predictor component")
  if(!identical(fr$stan_sha,digest::digest(file=stan_file,algo="sha256"))) stop("Stan source hash mismatch")
  dt<-data.table::copy(data.table::as.data.table(new_dt))
  cols<-unique(c("case_id","candidate_id","group",fr$feature_cols,fr$prep_state$hpo_col,fr$prep_state$net_col))
  if(!all(cols%in%names(dt))) stop("Missing scoring input columns")
  if(!nrow(dt)||anyNA(dt[,c("case_id","candidate_id","group"),with=FALSE])) stop("Empty input or missing identifiers")
  if(anyDuplicated(dt[,.(case_id,candidate_id)])) stop("Duplicate candidate identifiers")
  if(any(!dt$group%in%fr$group_levels)) stop("Unseen groups not supported by this adapter")
  if(any(!vapply(dt[,fr$feature_cols,with=FALSE],function(v) is.numeric(v)&&all(is.finite(v)),logical(1)))) stop("Invalid design features")
  if(nrow(dt)*posterior::ndraws(fr$draws)>max_draw_cells) stop("Batch too large; decrease --batch-cases")
  va<-apply_prep_state(dt,fr$prep_state); data.table::setorderv(va,c("case_id","candidate_id"))
  one<-data.table::data.table(case_id="dummy",candidate_id="dummy",label=1L,group=fr$group_levels[1],u_hpo=.5,u_network=.5,copula_available=0L)
  one[,(fr$feature_cols):=0]
  dat<-pack_stan_data(one,va,fr$feature_cols,fr$group_levels,fr$prior)
  gq<-model$generate_quantities(data=dat,fitted_params=fr$draws,seed=seed)
  pd<-gq$draws(variables="p_new",format="draws_matrix")
  nums<-as.integer(gsub("[^0-9]","",colnames(pd)))
  if(anyNA(nums)||!identical(sort(nums),seq_len(nrow(va)))) stop("Unexpected p_new columns")
  pd<-pd[,order(nums),drop=FALSE]
  ans<-data.table::copy(va[,.(case_id,candidate_id)])
  if("label"%in%names(va)) ans[,label:=va$label]
  ans[,probability:=colMeans(pd)]
  q<-t(apply(pd,2,quantile,probs=c(.025,.975),names=FALSE))
  ans[,`:=`(probability_lower_95=q[,1],probability_upper_95=q[,2],mean_rank=NA_real_,posterior_rank_first=NA_real_)]
  for(ii in split(seq_len(nrow(va)),va$case_id)) {
    a<-pd[,ii,drop=FALSE]
    if(any(!is.finite(a))||any(a<0|a>1)||max(abs(rowSums(a)-1))>1e-6) stop("Invalid probability draws")
    ranks<-do.call(rbind,lapply(seq_len(nrow(a)),function(s) rank(-a[s,],ties.method="average")))
    first<-a==apply(a,1,max)
    ans[ii,`:=`(mean_rank=colMeans(ranks),posterior_rank_first=colMeans(first/rowSums(first)))]
  }
  list(summary=ans,p_new_draws=pd,predictor_id=fr$predictor_id,candidate_protocol_id=fr$candidate_protocol_id,
       provisional_divergences=fr$provisional_divergences)
}
fit_split_conformal <- function(scored,fr,alpha=c(.10,.05)) {
  if(!identical(scored$predictor_id,fr$predictor_id)||!identical(scored$candidate_protocol_id,fr$candidate_protocol_id)) stop("Predictor/protocol mismatch")
  if(!length(fr$development_case_ids)) stop("Development-case manifest required")
  z<-data.table::as.data.table(scored$summary)
  if(!"label"%in%names(z)||anyNA(z$label)||any(!z$label%in%c(0,1))) stop("Calibration needs binary labels")
  if(any(as.character(z$case_id)%in%fr$development_case_ids)) stop("Calibration/development overlap")
  ck<-z[,.(positive=sum(label==1),total=sum(probability)),by=case_id]
  if(any(ck$positive!=1)||any(abs(ck$total-1)>1e-6)) stop("Invalid calibration cases")
  if(any(!is.finite(alpha)|alpha<=0|alpha>=1)) stop("Invalid alpha")
  s<-1-z[label==1,probability]; n<-length(s)
  if(!n) stop("No calibration cases")
  k<-ceiling((n+1)*(1-alpha)); q<-vapply(k,function(j) if(j>n) Inf else sort(s)[j],numeric(1))
  list(predictor_id=fr$predictor_id,candidate_protocol_id=fr$candidate_protocol_id,
    score="1_minus_CV_ensemble_mean_probability",n_calibration=n,
    calibration_case_ids=unique(as.character(z$case_id)),scores=s,
    thresholds=data.frame(alpha=alpha,k=k,qhat=q),
    assumption="Untouched labels and exchangeable deployment-comparable candidate construction required")
}
apply_split_conformal <- function(scored,calibration) {
  if(!identical(scored$predictor_id,calibration$predictor_id)||!identical(scored$candidate_protocol_id,calibration$candidate_protocol_id)) stop("Frozen predictor/protocol mismatch")
  if(any(as.character(scored$summary$case_id)%in%calibration$calibration_case_ids)) stop("Prediction/calibration overlap")
  z<-data.table::copy(data.table::as.data.table(scored$summary))
  for(i in seq_len(nrow(calibration$thresholds))) {
    nm<-paste0("in_conformal_set_",format(100*(1-calibration$thresholds$alpha[i]),trim=TRUE))
    data.table::set(z,j=nm,value=1-z$probability<=calibration$thresholds$qhat[i])
  }
  z
}
