#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1L])
.sf <- normalizePath(.sf, mustWork=TRUE)
.lib <- normalizePath(file.path(dirname(.sf), "..", "lib"), mustWork=TRUE)
source(file.path(.lib, "bootstrap.R"))
a <- parse_args(c("config", "sample", "vcf", "outdir"), c("mode", "chunk_records", "parity"))
cfg <- read_config(a$config)
mode <- a$mode %||% "selftest"
if (!mode %in% c("selftest", "build")) stop("mode must be selftest or build")
chunk_records <- suppressWarnings(as.numeric(a$chunk_records %||% "1000"))
if (length(chunk_records)!=1L || !is.finite(chunk_records) || chunk_records<1 || chunk_records!=floor(chunk_records) || chunk_records>20000) stop("chunk_records must be an integer from 1 to 20000")
chunk_records <- as.integer(chunk_records)
for (tool in c("sort", "bash", "gzip")) if (!nzchar(Sys.which(tool))) stop("Required command: ", tool)
reference_builder <- build_evidence_matrix
helper_paths <- sort(list.files(.lib, pattern="[.]R$", full.names=TRUE))
helper_hashes <- stats::setNames(lapply(helper_paths, sha256_file), helper_paths)
implementation_hash <- sha256_file(.sf)

make_record_builder <- function(fn) {
  pieces <- as.list(body(fn))
  header_hits <- records_hits <- 0L
  for (i in seq_along(pieces)) {
    if (identical(pieces[[i]], quote(hdr <- read_vcf_header(vcf)))) {
      pieces[[i]] <- quote(hdr <- .stream_header)
      header_hits <- header_hits + 1L
    }
    if (identical(pieces[[i]], quote(rec <- read_vcf_records(vcf)))) {
      pieces[[i]] <- quote(rec <- data.table::copy(.stream_records))
      records_hits <- records_hits + 1L
    }
  }
  if (header_hits!=1L || records_hits!=1L) stop("Unexpected builder implementation; refusing automatic adaptation")
  formals(fn) <- c(formals(fn), alist(.stream_header=, .stream_records=))
  body(fn) <- as.call(pieces)
  fn
}
record_builder <- make_record_builder(reference_builder)

qc_counts <- function(em, n_records, qc) {
  c(n_vcf_records=n_records, n_candidates=nrow(em),
    n_pass_core_qc=sum(em$pass_core_qc), n_missing_af=sum(is.na(em$MAX_AF)),
    n_fail_filter=sum(!em$pass_filter), n_fail_dp=sum(!em$pass_dp),
    n_fail_gq=sum(!em$pass_gq), n_fail_ad=sum(!em$pass_ad),
    n_fail_af=sum(!em$pass_af), n_fail_impact=sum(!em$pass_impact),
    n_pass_core_qc_dp_gt=sum(apply_core_qc(em, modifyList(qc, list(dp_rule="gt")))$pass_core_qc),
    n_pass_core_qc_dp_ge=sum(apply_core_qc(em, modifyList(qc, list(dp_rule="ge")))$pass_core_qc))
}

stream_build <- function(vcf, cfg, stage, batch_size) {
  dir.create(stage, recursive=TRUE, showWarnings=FALSE)
  con <- gzfile(vcf, "rt")
  on.exit(close(con), add=TRUE)
  header <- character()
  repeat {
    line <- readLines(con, n=1L, warn=FALSE)
    if (!length(line)) stop("VCF has no #CHROM header")
    header <- c(header, line)
    if (startsWith(line, "#CHROM")) break
  }
  columns <- strsplit(tail(header, 1L), "\t", fixed=TRUE)[[1L]]
  columns[1L] <- "CHROM"
  if (anyDuplicated(columns)) stop("Duplicate VCF column names")
  schema <- csq_schema_from_header(header)
  body_path <- file.path(stage, "body.tsv")
  sorted_path <- file.path(stage, "body.sorted.tsv")
  header_path <- file.path(stage, "header.tsv")
  matrix_path <- file.path(stage, "evidence_matrix.tsv.gz")
  total <- 0
  batch <- 0L
  counts <- NULL
  output_columns <- NULL
  af_cols <- NULL
  format_keys <- character()
  warned_af <- FALSE
  repeat {
    lines <- readLines(con, n=batch_size, warn=FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    if (!length(lines)) next
    fields <- strsplit(paste0(lines, "\t!"), "\t", fixed=TRUE)
    fields <- lapply(fields, function(x) x[-length(x)])
    if (any(lengths(fields)!=length(columns))) stop("VCF row/header field-count mismatch")
    rec <- data.table::as.data.table(data.table::transpose(fields))
    data.table::setnames(rec, columns)
    for (j in seq_along(rec)) data.table::set(rec, which(rec[[j]]==""), j, NA_character_)
    result <- withCallingHandlers(
      record_builder(vcf, cfg$qc, sample_id=cfg$vcf$sample_id,
        sample_index=cfg$vcf$sample_index,
        af_fields=unlist(cfg$qc$af_policy$fields), require_format=character(),
        .stream_header=header, .stream_records=rec),
      warning=function(w) {
        if (startsWith(conditionMessage(w), "no allele-frequency fields found")) {
          if (!warned_af) message("A1 WARNING: ", conditionMessage(w))
          warned_af <<- TRUE
          invokeRestart("muffleWarning")
        }
      })
    em <- result$matrix
    if (nrow(em)!=nrow(rec)) stop("Evidence/VCF row-count mismatch")
    if (!identical(names(em)[1L], "candidate_id")) stop("Expected candidate_id as first output column")
    if (anyNA(em$candidate_id) || any(grepl('[\t\r\n"]', em$candidate_id))) stop("Unsupported candidate ID for external ordering")
    if (is.null(output_columns)) {
      output_columns <- names(em)
      af_cols <- result$af_fields_used
      writeLines(paste(output_columns, collapse="\t"), header_path, useBytes=TRUE)
    } else if (!identical(output_columns, names(em)) || !identical(af_cols, result$af_fields_used)) stop("Chunk output schema or AF policy changed")
    values <- qc_counts(em, nrow(rec), cfg$qc)
    if (is.null(counts)) counts <- values else counts <- counts + values
    format_keys <- unique(c(format_keys, result$format_keys))
    data.table::fwrite(em, body_path, sep="\t", col.names=FALSE,
      append=batch>0L, na="")
    total <- total + nrow(rec)
    batch <- batch + 1L
    if (batch==1L || batch %% 10L==0L) message(sprintf("A1 streaming: batch %d; records %.0f; QC passes %.0f", batch, total, counts[["n_pass_core_qc"]]))
    rm(result, em, rec, fields, lines, values)
    invisible(gc(FALSE))
  }
  if (!batch) stop("No VCF data records")
  missing_format <- setdiff(unlist(cfg$vcf$required_format_fields), format_keys)
  if (length(missing_format)) stop("Required FORMAT fields absent: ", paste(missing_format, collapse=", "))
  status <- system2("sort", c("-s", "-t", shQuote("\t"), "-k1,1", "-S", "128M", "--parallel=1", "-T", shQuote(stage), "-o", shQuote(sorted_path), shQuote(body_path)), env="LC_ALL=C")
  if (status!=0L) stop("External evidence ordering failed")
  check <- file(sorted_path, "rt")
  n_sorted <- 0
  last_id <- NULL
  tryCatch({
    repeat {
      z <- readLines(check, n=10000L, warn=FALSE)
      if (!length(z)) break
      ids <- sub("\t.*$", "", z)
      if (anyDuplicated(ids) || (!is.null(last_id) && identical(ids[1L], last_id))) stop("Duplicate candidate IDs across input records/chunks")
      last_id <- tail(ids, 1L)
      n_sorted <- n_sorted + length(z)
    }
  }, finally=close(check))
  if (n_sorted!=total || counts[["n_candidates"]]!=total) stop("Final evidence row accounting failed")
  command <- paste("cat --", shQuote(header_path), shQuote(sorted_path), "| gzip -n -c >", shQuote(matrix_path))
  if (system2("bash", c("-o", "pipefail", "-c", shQuote(command)))!=0L) stop("Evidence compression failed")
  if (system2("gzip", c("-t", shQuote(matrix_path)))!=0L) stop("Evidence gzip validation failed")
  summary_path <- file.path(stage, "qc_summary.tsv")
  write_table(data.table(metric=names(counts), value=as.numeric(counts)), summary_path)
  list(matrix_file=matrix_path, summary_file=summary_path, counts=counts,
    schema=schema, af_fields_used=af_cols, format_keys=format_keys,
    records=total, batches=batch, output_columns=output_columns)
}

make_fixture <- function(vcf, out, n=128L) {
  input <- gzfile(vcf, "rt")
  on.exit(close(input), add=TRUE)
  hdr <- character()
  repeat {
    z <- readLines(input, n=1L, warn=FALSE)
    if (!length(z)) stop("Missing VCF header")
    hdr <- c(hdr, z)
    if (startsWith(z, "#CHROM")) break
  }
  records <- readLines(input, n=n, warn=FALSE)
  records <- records[nzchar(records) & !startsWith(records, "#")]
  if (!length(records)) stop("No fixture records")
  if (sum(nchar(records, type="bytes"))>16*1024^2) stop("Fixture exceeds 16 MiB; reduce fixture size before testing")
  output <- gzfile(out, "wt", compression=1)
  tryCatch(writeLines(c(hdr, records), output, useBytes=TRUE), finally=close(output))
  list(header=hdr, records=records)
}

if (mode=="selftest") {
  if (file.exists(a$outdir)) stop("Selftest output directory must be new")
  dir.create(a$outdir, recursive=TRUE)
  fixture <- file.path(a$outdir, "fixture.vcf.gz")
  fixture_data <- make_fixture(a$vcf, fixture)
  reference <- reference_builder(fixture, cfg$qc,
    sample_id=cfg$vcf$sample_id, sample_index=cfg$vcf$sample_index,
    af_fields=unlist(cfg$qc$af_policy$fields),
    require_format=unlist(cfg$vcf$required_format_fields))
  reference_file <- file.path(a$outdir, "reference.tsv.gz")
  write_table(reference$matrix, reference_file)
  expected <- read_gz(reference_file)
  expected_counts <- qc_counts(reference$matrix, reference$n_vcf_records, cfg$qc)
  for (size in c(17L, 53L)) {
    result <- stream_build(fixture, cfg, file.path(a$outdir, paste0("chunk_", size)), size)
    actual <- read_gz(result$matrix_file)
    comparison <- all.equal(as.data.frame(expected), as.data.frame(actual), tolerance=1e-12, check.attributes=TRUE)
    if (!isTRUE(comparison)) stop("Parity failure at chunk size ", size, ": ", paste(comparison, collapse="; "))
    if (!identical(names(expected_counts), names(result$counts)) || !isTRUE(all.equal(as.numeric(expected_counts), as.numeric(result$counts), tolerance=0))) stop("QC-count parity failure")
    message("PASS: matrix order, columns, values, missingness and QC counts; chunk size ", size)
  }
  duplicate_fixture <- file.path(a$outdir, "duplicate.vcf.gz")
  output <- gzfile(duplicate_fixture, "wt", compression=1)
  tryCatch(writeLines(c(fixture_data$header, fixture_data$records, fixture_data$records[1L]), output, useBytes=TRUE), finally=close(output))
  duplicate_error <- tryCatch({stream_build(duplicate_fixture, cfg, file.path(a$outdir,"duplicate_test"),17L); NULL}, error=function(e) conditionMessage(e))
  if (is.null(duplicate_error) || !grepl("duplicate", duplicate_error, ignore.case=TRUE)) stop("Duplicate rejection test failed")
  report <- list(status="passed", test_scope="bounded prefix against active modular builder, not historical full-run parity",
    fixture_records=length(fixture_data$records), fixture_sha256=sha256_file(fixture),
    chunk_sizes_tested=c(17L,53L), numeric_tolerance=1e-12,
    cross_chunk_duplicate_test=TRUE, implementation_sha256=implementation_hash,
    helpers=helper_hashes, config_sha256=cfg$config_sha256,
    source_variant_filter_added=FALSE, production_A1_written=FALSE)
  jsonlite::write_json(report, file.path(a$outdir,"parity_report.json"), pretty=TRUE, auto_unbox=TRUE)
  message("A1 STREAMING SELFTEST PASSED. No full-WGS build executed.")
} else {
  if (is.null(a$parity) || !file.exists(a$parity)) stop("A passed parity report is required")
  p <- jsonlite::read_json(a$parity)
  if (!identical(p$status,"passed") || !identical(p$implementation_sha256,implementation_hash) || !identical(p$config_sha256,cfg$config_sha256)) stop("Parity report/code/config mismatch")
  if (!identical(sort(names(p$helpers)), sort(names(helper_hashes)))) stop("Helper inventory differs from tested inventory")
  for (path in names(helper_hashes)) if (!identical(p$helpers[[path]],helper_hashes[[path]])) stop("A tested helper changed")
  od <- a$outdir
  dir.create(od, recursive=TRUE, showWarnings=FALSE)
  final_matrix <- file.path(od,"evidence_matrix.tsv.gz")
  final_summary <- file.path(od,"qc_summary.tsv")
  if (any(file.exists(c(final_matrix, final_summary, file.path(od,"manifest.json"))))) stop("Refusing existing A1 products; inspect incomplete outputs before restarting")
  stage <- tempfile(".a1_stream_", tmpdir=od)
  dir.create(stage)
  on.exit_cleanup <- function() unlink(stage, recursive=TRUE)
  input_hash <- sha256_file(a$vcf)
  result <- tryCatch(stream_build(a$vcf,cfg,stage,chunk_records), error=function(e) {on.exit_cleanup(); stop(e)})
  if (!identical(input_hash,sha256_file(a$vcf))) {on.exit_cleanup(); stop("VCF changed during execution")}
  if (!identical(implementation_hash,sha256_file(.sf))) {on.exit_cleanup(); stop("Streaming implementation changed")}
  for (path in names(helper_hashes)) if (!identical(helper_hashes[[path]],sha256_file(path))) {on.exit_cleanup(); stop("A helper changed during execution")}
  if (!file.rename(result$matrix_file,final_matrix) || !file.rename(result$summary_file,final_summary)) stop("Could not publish staged A1 files")
  write_node_manifest(od,"a1",cfg,
    inputs=c(a$vcf,a$parity,.sf,helper_paths), outputs=c(final_matrix,final_summary),
    extra=list(io_backend="single_pass_chunks_disk_spool_external_order",
      chunk_records=chunk_records, sort_memory="128M", n_batches=result$batches,
      n_vcf_records=result$records, af_fields_used=as.list(result$af_fields_used),
      dp_rule=cfg$qc$dp_rule, n_pass_core_qc=unname(result$counts[["n_pass_core_qc"]]),
      phase_fields_present=as.list(intersect(result$format_keys,c("PGT","PID","PS"))),
      cross_chunk_duplicates_checked=TRUE, added_variant_filter=FALSE,
      output_order="candidate_id; LC_ALL=C; checked against active-builder fixture",
      historical_full_run_parity_established=FALSE))
  on.exit_cleanup()
  message("A1 STREAMING BUILD COMPLETE: records=",result$records)
}
