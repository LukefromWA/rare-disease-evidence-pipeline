# A11-specific execution/provenance helpers; no scientific calculations.
a11_sha <- function(p) digest::digest(file=p, algo="sha256")
a11_cli <- function(required, defaults=list()) {
  x <- commandArgs(trailingOnly=TRUE)
  if (length(x) %% 2L) stop("Use --key value pairs")
  k <- sub("^--", "", x[seq.int(1L, length(x), by=2L)])
  if (!length(x)) k <- character()
  if (length(x) && any(!startsWith(x[seq.int(1L,length(x),by=2L)], "--"))) stop("Invalid argument name")
  if (anyDuplicated(k) || any(!k %in% c(required,names(defaults)))) stop("Duplicate or unknown argument")
  z <- defaults
  if (length(x)) for (i in seq_along(k)) z[[k[i]]] <- x[2L*i]
  for (n in required) if (is.null(z[[n]]) || !nzchar(z[[n]])) stop("Required --",n)
  z
}
a11_integer <- function(x, name, minimum=0L) {
  n <- suppressWarnings(as.numeric(x))
  if (length(n)!=1L || !is.finite(n) || n!=floor(n) || n<minimum || n>.Machine$integer.max) stop("Invalid ",name)
  as.integer(n)
}
for (.p in c("data.table","digest","jsonlite")) if (!requireNamespace(.p,quietly=TRUE)) stop("Install required package before execution: ",.p)
a11_begin <- function(outdir, inputs, script, runtime) {
  if (file.exists(outdir)) stop("Refusing existing output path: ",outdir)
  if (any(!file.exists(inputs))) stop("Missing declared input")
  if (!dir.create(outdir,recursive=TRUE)) stop("Cannot create output directory")
  list(inputs=inputs, hashes=vapply(inputs,a11_sha,character(1)), script=script,
       script_sha=a11_sha(script), runtime=runtime, runtime_sha=a11_sha(runtime))
}
a11_finish <- function(ctx, outdir, kind, parameters) {
  if (!identical(vapply(ctx$inputs,a11_sha,character(1)),ctx$hashes) ||
      !identical(a11_sha(ctx$script),ctx$script_sha) || !identical(a11_sha(ctx$runtime),ctx$runtime_sha)) stop("Input or source changed during execution")
  files <- list.files(outdir,full.names=TRUE)
  files <- files[!dir.exists(files) & basename(files)!="run_provenance.json"]
  outputs <- setNames(lapply(files,function(p) list(sha256=a11_sha(p))),basename(files))
  inputs <- setNames(lapply(seq_along(ctx$inputs),function(i) list(sha256=unname(ctx$hashes[i]))),names(ctx$inputs))
  versions <- setNames(lapply(c("data.table","digest","jsonlite"),function(p) as.character(utils::packageVersion(p))),c("data.table","digest","jsonlite"))
  identity <- list(kind=kind,inputs=inputs,parameters=parameters,script_sha=ctx$script_sha,
                   runtime_sha=ctx$runtime_sha,R=R.version.string,RNGkind=RNGkind(),packages=versions)
  manifest <- c(identity,list(run_id=digest::digest(identity,algo="sha256"),outputs=outputs,
    time_utc=format(Sys.time(),tz="UTC",usetz=TRUE),locale=Sys.getlocale(),
    session_info=capture.output(sessionInfo())))
  jsonlite::write_json(manifest,file.path(outdir,"run_provenance.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA,null="null")
}
a11_verify <- function(manifest, evidence, consensus, copula=NULL) {
  if (!file.exists(manifest)) stop("Missing A11 producer provenance")
  m <- jsonlite::read_json(manifest,simplifyVector=FALSE)
  if (!identical(m$kind,"a11_dependence_sensitivity")) stop("Wrong A11 producer provenance kind")
  if (!identical(m$inputs$evidence$sha256,a11_sha(evidence))) stop("A11 evidence/provenance mismatch")
  expected <- c(bayesian_dependence_sensitivity_ensemble_consensus.csv=consensus)
  if (!is.null(copula)) expected <- c(expected,copula_dependence_summary.csv=copula)
  for (n in names(expected)) if (!identical(m$outputs[[n]]$sha256,a11_sha(expected[[n]]))) stop("A11 output/provenance mismatch: ",n)
  m
}
