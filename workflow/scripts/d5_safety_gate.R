#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

a <- parse_args(
    c("config", "sample", "evidence", "outdir", "safety", "targets")
)

cfg <- read_config(a$config)
b <- require_bundle(cfg)
d <- read_gz(a$evidence)
sf <- read_gz(a$safety, colClasses = "character")

need_cols(sf, c(
    "drug_id", "pediatric_evidence_status", "contraindication_status",
    "route_feasible", "review_state"
), "case safety")
if (anyDuplicated(sf$drug_id)) stop("Duplicate safety drug IDs")

mf <- file.path(dirname(a$safety), "run_manifest.json")
sums <- file.path(dirname(a$safety), "checksums.sha256")
m <- jsonlite::read_json(mf)

if (!identical(m$drugs_input_sha256, sha256_file(a$targets)))
    stop("Safety evidence was not acquired from THIS D2 target table")

pairs <- strsplit(readLines(sums, warn = FALSE), "  ", fixed = TRUE)
hit <- vapply(pairs, function(x) {
    length(x) == 2L && identical(x[2L], basename(a$safety))
}, logical(1))

if (sum(hit) != 1L ||
    !identical(pairs[[which(hit)]][1L], sha256_file(a$safety)))
    stop("Case safety checksum mismatch")

sf[is.na(review_state) | review_state != "reviewed", `:=`(
    pediatric_evidence_status = "unknown",
    contraindication_status = "unknown",
    route_feasible = NA_character_
)]

old <- intersect(names(d), c(
    "pediatric_evidence_status", "contraindication_status",
    "route_feasible", "safety_gate", "eligibility", "exclusion_reason"
))
if (length(old)) d[, (old) := NULL]

d <- merge(d, sf[, .(
    drug_id, pediatric_evidence_status,
    contraindication_status, route_feasible
)], by = "drug_id", all.x = TRUE)

d[, safety_gate := safety_gate(
    pediatric_evidence_status, contraindication_status, route_feasible
)]
d[, eligibility := final_eligibility(
    direction_gate, safety_gate, isTRUE(cfg$features$pediatric_mode)
)]
d[, exclusion_reason := fcase(
    direction_gate %in% "fail", "direction_incompatible",
    safety_gate %in% "fail", "safety_or_feasibility",
    safety_gate %in% "not_configured",
        "safety_unknown_unreviewed_or_incomplete",
    !direction_gate %in% "pass", "direction_unknown",
    default = NA_character_
)]

od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)

write_table(d, file.path(od, "gated.tsv.gz"))
write_table(d[eligibility == "excluded", .(
    scenario_id, drug_id, direction_gate, safety_gate, exclusion_reason
)], file.path(od, "excluded.tsv"))
write_table(d[eligibility != "eligible_for_research_tiering", .(
    scenario_id, drug_id, direction_gate, safety_gate,
    eligibility, exclusion_reason
)], file.path(od, "held_or_excluded.tsv"))

write_node_manifest(
    od, "d5", cfg,
    inputs = c(a$evidence, a$safety, a$targets, mf, sums),
    outputs = file.path(od, c(
        "gated.tsv.gz", "excluded.tsv", "held_or_excluded.tsv"
    )),
    extra = list(eligibility_counts = as.list(table(d$eligibility)))
)

stopifnot(
    all(d[direction_gate %in% "fail" | safety_gate %in% "fail",
          eligibility] == "excluded"),
    !any(d[safety_gate != "pass",
           eligibility] %in% "eligible_for_research_tiering")
)
