#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# E00: registry gate. Strict policy: url/version/license must be recorded; frozen raw snapshots must exist when live downloads are off.
a <- parse_args(c("registry", "root", "out"), c("scope")); reg <- yaml::read_yaml(a$registry)$track2_resources

scope <- a$scope %||% "full"
if (!scope %in% c("core", "full")) stop("Invalid resource scope")
if (scope == "core") reg$sources$pediatric_safety <- NULL
validate_registry(reg, a$root, strict = isTRUE(reg$strict_resource_policy)); dir.create(dirname(a$out), recursive = TRUE, showWarnings = FALSE); file.create(a$out); logmsg("registry OK")
