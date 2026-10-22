# workflow/lib/classes.R -- optional S4 contracts and R6 workers
suppressPackageStartupMessages({ library(methods); library(data.table); library(R6) })
setOldClass(c("data.table", "data.frame"))

setClass("CandidateTable", slots = c(data = "data.table", stage = "character", required = "character"),
  validity = function(object) {
    d <- object@data; err <- character(); miss <- setdiff(object@required, names(d))
    if (length(miss)) err <- c(err, paste("missing columns:", paste(miss, collapse = ", ")))
    if ("candidate_id" %in% names(d) && anyDuplicated(d$candidate_id)) err <- c(err, "candidate_id not unique")
    if (length(err)) err else TRUE })
CandidateTable <- function(dt, stage, required = "candidate_id")
  new("CandidateTable", data = as.data.table(dt), stage = stage, required = required)

setClass("HypothesisSet", slots = c(data = "data.table"),
  validity = function(object) {
    d <- object@data; if (length(setdiff(c("gene", "variant", "q_j", "p_m"), names(d)))) return("missing columns")
    if (abs(sum(unique(d[, .(variant, q_j)])$q_j) - 1) > 1e-8) return("sum(q_j) != 1")
    if (any(d$p_m < 0 | d$p_m > 1)) return("p_m outside [0,1]"); TRUE })
setGeneric("scenarios", function(x, ...) standardGeneric("scenarios"))
setMethod("scenarios", "HypothesisSet", function(x, ...) {
  d <- copy(x@data)[, w := q_j * p_m][w > 0]; d[, scenario_weight := w / sum(w)][] })

Node <- R6Class("Node", public = list(
  id = NULL, inputs = list(), required_cols = "candidate_id",
  initialize = function(id, inputs = list()) { self$id <- id; self$inputs <- inputs },
  run = function(output_path) {
    missing <- unlist(self$inputs)[!file.exists(unlist(self$inputs))]; if (length(missing)) stop("Missing input(s): ", paste(missing, collapse = ", "))
    res <- self$compute(); self$validate(res); data.table::fwrite(res, output_path, sep = "\t"); invisible(res) },
  validate = function(res) { stopifnot(all(self$required_cols %in% names(res)), !anyDuplicated(res$candidate_id)) },
  compute = function() stop("abstract: implement compute() in the subclass")))

FrozenEnsemble <- R6Class("FrozenEnsemble", public = list(
  initialize = function(coef, weights) { private$coef <- coef; private$weights <- weights; lockEnvironment(private, bindings = TRUE) },
  predict = function(X, case_id) predict_softmax(X, private$coef, case_id),
  fit = function(...) stop("refit is forbidden on a frozen ensemble")),
  private = list(coef = NULL, weights = NULL))
