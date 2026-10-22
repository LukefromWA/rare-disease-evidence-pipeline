#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# A0.4 manifest: checksums, tool versions, VEP command and flag check, CSQ schema, FORMAT and phase fields, contigs, mode.
a <- parse_args(c("config", "mode", "sample", "vcf", "tbi", "source_vcf", "metadata", "outdir"), c("raw_vcf", "raw_tbi", "fasta", "fai"))
cfg <- read_config(a$config); od <- a$outdir
hdr <- read_vcf_header(a$vcf); schema <- csq_schema_from_header(hdr)
rec1 <- data.table::fread(cmd = paste("gzip -cd", shQuote(a$vcf), "| grep -v '^#' | head -n 2000"), sep = "\t", header = FALSE, colClasses = "character", quote = "", fill = TRUE)
format_keys <- if (nrow(rec1)) unique(unlist(strsplit(rec1$V9, ":", fixed = TRUE))) else character()
contigs <- data.table::fread(cmd = paste("gzip -cd", shQuote(a$vcf), "| grep -v '^#' | cut -f1 | sort | uniq -c | sort -k1,1nr"), header = FALSE, sep = " ", fill = TRUE)
miss_fmt <- setdiff(unlist(cfg$vcf$required_format_fields %||% list("GT", "DP", "GQ", "AD")), format_keys)
if (length(miss_fmt)) stop("A0 STOP: required FORMAT field(s) absent: ", paste(miss_fmt, collapse = ", "), call. = FALSE)
meta_txt <- paste(readLines(a$metadata, warn = FALSE), collapse = "\n")
flags <- unlist(cfg$vep$required_flags %||% list()); flags <- flags[startsWith(flags, "--")]
missing_flags <- flags[!vapply(flags, function(f) grepl(f, meta_txt, fixed = TRUE), TRUE)]
if (length(missing_flags) && isTRUE(cfg$vep$require_flags_in_metadata %||% TRUE) && a$mode == "existing_vep")
  stop("A0 STOP: recorded VEP command/metadata lacks required flag(s): ", paste(missing_flags, collapse = " "), call. = FALSE)
ver <- regmatches(meta_txt, regexpr("(?i)vep[^0-9]{0,20}[0-9]+(\\.[0-9]+)*", meta_txt, perl = TRUE))
if (a$mode == "existing_vep" && !length(ver)) logmsg("A0 warning: no VEP version string found in metadata")
exp <- cfg$vep$version_expected; if (!is.null(exp) && length(ver) && !grepl(as.character(exp), ver[1], fixed = TRUE)) stop("A0 STOP: observed VEP version '", ver[1], "' does not match expected ", exp, call. = FALSE)
inputs <- c(a$source_vcf, if (!is.na(a$raw_vcf %||% NA) && a$raw_vcf != "NA") a$raw_vcf, if (!is.na(a$raw_tbi %||% NA) && a$raw_tbi != "NA") a$raw_tbi, a$metadata, a$fasta %||% NA, a$fai %||% NA)
writeLines(schema, file.path(od, "csq_schema.txt"))
write_node_manifest(od, "a0", cfg, inputs = inputs, outputs = c(a$vcf, a$tbi, file.path(od, "csq_schema.txt")), file = "a0_manifest.json",
  extra = list(sample = a$sample, run_mode = a$mode, vep_version_observed = if (length(ver)) ver[1] else NA, vep_cache_version = cfg$vep$cache_version, required_flags_missing = as.list(missing_flags),
               csq_schema = as.list(schema), format_fields = as.list(format_keys), phase_fields_available = as.list(intersect(format_keys, c("PGT", "PID", "PS"))),
               contigs = lapply(seq_len(nrow(contigs)), function(i) list(n = contigs$V1[i], contig = contigs$V2[i])), n_records_scanned_for_format = nrow(rec1),
               chrM_present = any(grepl("^(chr)?M(T)?$", contigs$V2)), reference_fasta = a$fasta %||% NA, reference_fai = a$fai %||% NA))
logmsg("A0 manifest written: ", file.path(od, "a0_manifest.json"))
