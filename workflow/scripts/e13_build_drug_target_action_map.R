#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

a <- parse_args(
    c("registry", "raw"),
    c("lookup", "out", "native_dir", "prepare_only")
)

if (identical(a$prepare_only, "true")) {
    if (is.null(a$native_dir)) stop("Preparation needs --native_dir")
    prepare_chembl_drug_actions(a$native_dir, a$raw)
} else {
    if (is.null(a$lookup) || is.null(a$out))
        stop("Mapping needs --lookup and --out")

    reg <- yaml::read_yaml(a$registry)$track2_resources
    s <- reg$sources$drug_targets
    assert_frozen_resource(a$raw, "drug_targets_raw")
    assert_frozen_resource(a$lookup, "uniprot_gene_lookup")

    mapped <- build_drug_target_action_map(
        read_gz(a$raw, colClasses = "character"),
        read_gz(a$lookup),
        s$version %||% "unknown"
    )

    # Opposing mapped actions for the same molecule/gene are held,
    # rather than passed to D2 as unambiguous direct actions.
    conflicts <- mapped[
        mapping_status == "mapped" & !is.na(drug_effect_sign),
        .(n_signs = uniqueN(drug_effect_sign)),
        by = .(drug_id, target_gene)
    ][n_signs > 1L]

    if (nrow(conflicts)) {
        mapped[conflicts, on = .(drug_id, target_gene),
            `:=`(
                mapping_status = "conflicting_actions",
                action_confidence = "conflicting"
            )
        ]
    }

    write_table(mapped, a$out)
    write_resource_manifest(
        a$out, "drug_target_action_map",
        s$url, s$version, s$license,
        notes = paste0(
            "Raw input sha256=", sha256_file(a$raw),
            "; rejected native mechanisms retained beside raw input; ",
            "opposing mapped molecule/gene actions flagged, not ranked."
        )
    )
    logmsg(
        "Drug mapping rows: ", nrow(mapped),
        "; mapped: ", sum(mapped$mapping_status == "mapped"),
        "; conflicting molecule/gene pairs: ", nrow(conflicts)
    )
}
