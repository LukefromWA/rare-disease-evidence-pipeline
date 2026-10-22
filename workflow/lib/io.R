# workflow/lib/io.R -- config, paths, safe table IO, argument parsing, loud assertions (pure helpers; no scientific logic)
suppressPackageStartupMessages({ library(data.table); library(yaml); library(jsonlite) })

`%||%` <- function(a, b) if (is.null(a)) b else a
logmsg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)), file = stderr())

# ---- arguments: every node script takes explicit --name value pairs ----
parse_args <- function(required, optional = character(), argv = commandArgs(trailingOnly = TRUE)) {
  if (length(argv) %% 2L != 0L) stop("arguments must be --name value pairs: ", paste(argv, collapse = " "), call. = FALSE)
  keys <- sub("^--", "", argv[seq(1, length(argv), 2)]); vals <- argv[seq(2, length(argv), 2)]; a <- stats::setNames(as.list(vals), keys)
  miss <- setdiff(required, names(a))
  if (length(miss)) stop("missing argument(s): ", paste0("--", miss, collapse = " "), "\nusage: --", paste(required, collapse = " <v> --"), " <v>", call. = FALSE)
  unk <- setdiff(names(a), c(required, optional)); if (length(unk)) stop("unknown argument(s): ", paste(unk, collapse = ", "), call. = FALSE)
  a
}

# ---- config: config.yaml overlaid by config.private.yaml (private absolute paths only); root = MVA_ROOT or cwd ----
merge_lists <- function(base, over) { for (k in names(over)) base[[k]] <- if (is.list(base[[k]]) && is.list(over[[k]]) && !is.null(names(over[[k]]))) merge_lists(base[[k]], over[[k]]) else over[[k]]; base }
read_config <- function(config, private = NULL) {
  need_files(config, "config"); cfg <- yaml::read_yaml(config)
  priv <- private %||% if (basename(config) == "config.yaml") file.path(dirname(config), "config.private.yaml") else NULL      # toy configs are self-contained
  if (!is.null(priv) && file.exists(priv)) cfg <- merge_lists(cfg, yaml::read_yaml(priv))
  cfg$root <- normalizePath(Sys.getenv("MVA_ROOT", unset = "."), mustWork = TRUE); cfg$config_path <- normalizePath(config)
  cfg$config_sha256 <- digest::digest(file = config, algo = "sha256"); cfg
}
rp <- function(cfg, path) { if (is.null(path)) return(character()); p <- as.character(unlist(path)); ifelse(is.na(p), NA_character_, ifelse(grepl("^(/|[A-Za-z]:)", p), p, file.path(cfg$root, p))) }
need_files <- function(paths, what = "input") {
  paths <- paths[!is.na(paths)]; miss <- paths[!file.exists(paths)]
  if (length(miss)) stop("Missing ", what, ": ", paste(miss, collapse = "; "), call. = FALSE); invisible(paths)
}
need_cols <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) stop(label, " is missing column(s): ", paste(miss, collapse = ", "), "\n  found: ", paste(names(x), collapse = ", "), call. = FALSE); invisible(x)
}
assert_columns <- function(x, required, label = "table") need_cols(x, required, label)

# ---- safe table IO (gz via pipe so R.utils is not required) ----
read_gz <- function(path, ...) {
  need_files(path)
  if (grepl("\\.gz$", path)) data.table::fread(cmd = paste("gzip -cd", shQuote(path)), ...) else data.table::fread(path, ...)
}
read_table_auto <- read_gz
write_table <- function(x, path) { dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); data.table::fwrite(x, path, sep = if (grepl("\\.csv", path)) "," else "\t"); invisible(path) }
node_dir <- function(cfg, sample, node) { d <- file.path(rp(cfg, cfg$run$output_root), sample, node); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
ext_dir  <- function(cfg, node) { d <- file.path(rp(cfg, cfg$run$external_root %||% "results/external"), node); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
