#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

# Structural checks are not statistical validation of the rescue method.
a <- parse_args(c("graph", "sentinel"), c("report"))

if (file.exists(a$sentinel))
    stop("An old propagation sentinel exists; archive it before this check.")

report <- data.table::data.table(
    check = c(
        "production_target_gene_wiring",
        "rescue_readout_and_convergence",
        "documented_positive_control",
        "rescue_statistic_null_comparison"
    ),
    pass = FALSE,
    status = "not_validated",
    interpretation = "Enablement criteria not established; no sentinel written"
)

if (!is.null(a$report)) {
    dir.create(dirname(a$report), recursive = TRUE, showWarnings = FALSE)
    write_table(report, a$report)
}
logmsg("D3.6 remains disabled: rescue validation is not complete")
