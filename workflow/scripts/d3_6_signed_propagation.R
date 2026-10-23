#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

# Explicitly unavailable until the rescue engine and its controls are validated.
a <- parse_args(c("config", "sample", "targets", "outdir"))
cfg <- read_config(a$config)
tu <- read_gz(a$targets)
need_cols(tu, c("scenario_id", "drug_id"), "target universe")

od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)

out <- unique(tu[, .(scenario_id, drug_id)])
out[, `:=`(
    rescue_mean = NA_real_,
    rescue_lo = NA_real_,
    rescue_hi = NA_real_,
    p_rescue_gt0 = NA_real_,
    propagation_status = "disabled_pending_rescue_validation"
)]

write_table(out, file.path(od, "method6_propagation.tsv.gz"))
write_method_status(
    od,
    "d3_6_signed_propagation",
    "disabled_pending_rescue_validation",
    paste(
        "Requires distinct disease-gene/intervention-target wiring,",
        "validated readout and convergence, documented positive control,",
        "and rescue-statistic comparison against appropriate nulls.",
        "A degree-preservation check alone is not sufficient."
    ),
    cfg,
    inputs = a$targets,
    node = "d3_6"
)
