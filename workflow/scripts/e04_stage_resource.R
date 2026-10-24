#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("native", "native_sha256", "out", "source_url", "source_version", "license"))
if (file.exists(a$out) || file.exists(paste0(a$out, ".resource_manifest.json")))
  stop("Refusing existing staged resource")
expected <- tolower(trimws(a$native_sha256))
if (!identical(sha256_file(a$native), expected)) stop("Native checksum mismatch")
dir.create(dirname(a$out), recursive=TRUE, showWarnings=FALSE)
if (!file.copy(a$native, a$out, overwrite=FALSE)) stop("Native staging failed")
if (!identical(sha256_file(a$out), expected) ||
    !identical(sha256_file(a$native), expected)) stop("Staging checksum mismatch")
write_resource_manifest(a$out, "reactome_native", a$source_url,
  a$source_version, a$license,
  notes=paste0("Unmodified native snapshot; source sha256=", expected))
assert_frozen_resource(a$out, "reactome_native")
logmsg("Verified native resource staged")
