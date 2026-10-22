# workflow/lib/manifest.R -- checksums, tool versions, command capture, one manifest per node, boundary assertions
suppressPackageStartupMessages({ library(digest); library(jsonlite) })

sha256_file <- function(path) { if (!file.exists(path)) stop("Missing file for checksum: ", path); digest::digest(file = path, algo = "sha256") }
tool_version <- function(cmd, arg = "--version") tryCatch(trimws(paste(system2(cmd, arg, stdout = TRUE, stderr = TRUE)[1], collapse = " ")), error = function(e) NA_character_, warning = function(w) NA_character_)
tool_versions <- function() list(R = R.version.string, bcftools = tool_version("bcftools"), vep = tool_version("vep"), tabix = tool_version("tabix"),
                                  python = tool_version("python3"), data.table = as.character(utils::packageVersion("data.table")))

# One manifest.json per node directory. Records inputs and outputs with sha256, config hash, seed, versions, command.
write_node_manifest <- function(dir, node, cfg, inputs = character(), outputs = character(), extra = list(), file = "manifest.json") {
  fin <- function(p) { p <- p[!is.na(p) & file.exists(p)]; if (!length(p)) return(list()); lapply(stats::setNames(p, basename(p)), function(f) list(path = normalizePath(f), sha256 = sha256_file(f), size_bytes = file.info(f)$size)) }
  m <- c(list(node = node, run_id = cfg$run$run_id, mode = cfg$run$mode, seed = cfg$run$seed, config_sha256 = cfg$config_sha256, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
              command = paste(commandArgs(), collapse = " "), tool_versions = tool_versions(), inputs = fin(inputs), outputs = fin(outputs)), extra)
  assert_no_private_in_public(m, node)
  jsonlite::write_json(m, file.path(dir, file), pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"); invisible(file.path(dir, file))
}
# Layer B and Track 2 resource manifests must never mention private paths.
assert_no_private_in_public <- function(m, node) {
  if (!grepl("^(b|e)", node)) return(invisible(TRUE))
  if (any(grepl("data/private|results/private|WGS_EX|MVA_WGS/data", unlist(m[setdiff(names(m), c("command"))]), ignore.case = TRUE))) stop("private reference in a public/external manifest: ", node, call. = FALSE)
  invisible(TRUE)
}
read_manifest <- function(path) jsonlite::read_json(need_files(path), simplifyVector = FALSE)
assert_frozen <- function(path, expected_sha) { if (!identical(sha256_file(path), expected_sha)) stop("checksum mismatch (artifact changed after freeze): ", path, call. = FALSE); invisible(TRUE) }
write_bundle_checksums <- function(paths, out, root = ".") {          # `sha256sum -c` compatible
  paths <- sort(unique(paths[file.exists(paths)])); rel <- sub(paste0("^", normalizePath(root), "/"), "", normalizePath(paths))
  writeLines(sprintf("%s  %s", vapply(paths, sha256_file, ""), rel), out); invisible(out)
}

# ---- Track 2 resource manifests (case-independent; one per raw/derived artifact) ----
resource_cfg <- function(registry_path) list(run = list(run_id = "track2_resources", mode = "resources", seed = NA), config_sha256 = sha256_file(registry_path))
write_resource_manifest <- function(path, resource_id, source_url = NA, source_version = NA, license = NA, schema_version = "1", notes = NA_character_) {
  if (!file.exists(path)) stop("Cannot manifest missing resource: ", path)
  e <- list(resource_id = resource_id, local_path = normalizePath(path), source_url = source_url, source_version = source_version, license = license, schema_version = schema_version,
            sha256 = sha256_file(path), size_bytes = file.info(path)$size, timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), notes = notes)
  jsonlite::write_json(e, paste0(path, ".resource_manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null"); invisible(e)
}
assert_frozen_resource <- function(path, id = basename(path)) {
  m <- paste0(path, ".resource_manifest.json"); if (!file.exists(path)) stop("Missing resource: ", id, call. = FALSE); if (!file.exists(m)) stop("Missing resource manifest: ", id, call. = FALSE)
  if (!identical(jsonlite::read_json(m)$sha256, sha256_file(path))) stop("Checksum mismatch: ", id, call. = FALSE); invisible(TRUE)
}
