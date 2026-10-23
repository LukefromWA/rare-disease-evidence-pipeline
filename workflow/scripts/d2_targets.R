#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

# BEGIN d2_model_views_v1
# Preserve discovery associations and weights; add separate model provenance.

a <- parse_args(c(
    "config", "sample", "mechanism", "model_views", "hierarchy_dir", "outdir"
))
cfg <- read_config(a$config)
b <- require_bundle(cfg)
sc <- read_gz(a$mechanism)
mv <- read_gz(a$model_views)

need_cols(sc, c(
    "hypothesis_id", "scenario_id", "gene", "variant"
), "D1 mechanism scenarios")

value_columns <- c(
    "view_status", "score_type", "score", "rank",
    "p_rank1", "p_top3", "p_top5", "p_top10",
    "source_field", "source_sha256"
)
need_cols(mv, c(
    "hypothesis_id", "gene", "variant", "model_view", value_columns
), "D0 model-view ledger")

if (anyDuplicated(mv[, .(hypothesis_id, model_view)]))
    stop("Duplicate hypothesis/model-view records")

expected_views <- c(
    "stan", "a10_informed", "a10_agnostic", "dependence_ensemble"
)
if (!setequal(unique(mv$model_view), expected_views))
    stop("Unexpected model-view set")

h <- unique(sc[, .(hypothesis_id, gene, variant)])
if (anyDuplicated(h$hypothesis_id))
    stop("Inconsistent hypothesis identities in D1")

view_keys <- unique(mv[, .(hypothesis_id, gene, variant)])
if (anyDuplicated(view_keys$hypothesis_id))
    stop("Inconsistent hypothesis identities in model-view ledger")

if (!setequal(h$hypothesis_id, view_keys$hypothesis_id))
    stop("D1 and model-view hypothesis universes differ")

vi <- match(h$hypothesis_id, view_keys$hypothesis_id)
if (!all(h$gene == view_keys$gene[vi]) ||
    !all(h$variant == view_keys$variant[vi]))
    stop("D1 and model-view gene/variant identities differ")

view_counts <- mv[, .N, by = hypothesis_id]
if (any(view_counts$N != length(expected_views)))
    stop("Incomplete model-view coverage per hypothesis")

# One wide row per hypothesis: never multiply associations by model count.
wide <- data.table::dcast(
    mv,
    hypothesis_id ~ model_view,
    value.var = value_columns
)
data.table::setnames(
    wide,
    setdiff(names(wide), "hypothesis_id"),
    paste0("mv_", setdiff(names(wide), "hypothesis_id"))
)
if (anyDuplicated(wide$hypothesis_id))
    stop("Model-view pivot produced duplicate hypothesis keys")

dr <- read_gz(b$drugs)
need_cols(dr, c(
    "drug_id", "drug_name", "target_gene", "action_normalized",
    "mapping_status", "drug_effect_sign"
), "drug_target_action_map")

n_un <- sum(dr$mapping_status != "mapped", na.rm = TRUE)

# Original D2 association construction: no score threshold or reweighting.
tu <- merge(
    sc,
    dr[
        mapping_status == "mapped" & !is.na(target_gene),
        .(
            drug_id, drug_name, target_gene,
            action_normalized, drug_effect_sign
        )
    ],
    by.x = "gene",
    by.y = "target_gene",
    allow.cartesian = TRUE
)

# Retain an exact reference for association/weight preservation checks.
before <- data.table::copy(tu)
tu[, .mv_original_row := seq_len(.N)]
tu <- merge(
    tu, wide,
    by = "hypothesis_id",
    all.x = TRUE,
    sort = FALSE
)
data.table::setorder(tu, .mv_original_row)
tu[, .mv_original_row := NULL]

if (nrow(tu) != nrow(before))
    stop("Model annotation changed the number of association rows")
for (column in names(before)) {
    if (!identical(tu[[column]], before[[column]]))
        stop("Model annotation changed original column: ", column)
}

# Include hypotheses with zero mapped direct associations.
counts <- tu[, .(
    n_direct_drugs = data.table::uniqueN(drug_id),
    n_association_rows = .N
), by = hypothesis_id]

audit <- merge(h, wide, by = "hypothesis_id", all.x = TRUE, sort = FALSE)
audit <- merge(
    audit, counts,
    by = "hypothesis_id", all.x = TRUE, sort = FALSE
)
audit[is.na(n_direct_drugs), n_direct_drugs := 0L]
audit[is.na(n_association_rows), n_association_rows := 0L]
audit[, direct_association_status := data.table::fifelse(
    n_direct_drugs > 0,
    "mapped_direct_association_found",
    "no_mapped_direct_association"
)]
stopifnot(
    nrow(audit) == nrow(h),
    !anyDuplicated(audit$hypothesis_id)
)


# BEGIN d2_pathway_discovery_v1
# Exploratory shared-pathway discovery, separate from direct associations.

# BEGIN d2_hierarchy_annotations_v1
hierarchy_dir <- a$hierarchy_dir
hierarchy_files <- c(
    "ReactomePathways.txt",
    "ReactomePathwaysRelation.txt",
    "human_pathway_hierarchy.tsv",
    "manifest.json"
)
checksum_file <- file.path(hierarchy_dir, "checksums.sha256")
checksum_rows <- strsplit(
    readLines(checksum_file, warn = FALSE), "  ", fixed = TRUE
)
if (!length(checksum_rows) ||
    any(lengths(checksum_rows) != 2L))
    stop("Malformed supplemental checksums")

checksum_names <- vapply(checksum_rows, `[`, character(1), 2L)
checksum_values <- vapply(checksum_rows, `[`, character(1), 1L)
if (anyDuplicated(checksum_names) ||
    !setequal(checksum_names, hierarchy_files))
    stop("Unexpected supplemental checksum file list")

for (j in seq_along(checksum_names)) {
    path <- file.path(hierarchy_dir, checksum_names[j])
    if (!file.exists(path) ||
        sha256_file(path) != checksum_values[j])
        stop("Supplemental checksum mismatch: ", checksum_names[j])
}

hierarchy_meta <- jsonlite::fromJSON(
    file.path(hierarchy_dir, "manifest.json")
)
if (!identical(
    sha256_file(b$pathways),
    hierarchy_meta$existing_gene_set_sha256
))
    stop("Hierarchy snapshot references different gene sets")

hierarchy_table <- data.table::fread(
    file.path(hierarchy_dir, "human_pathway_hierarchy.tsv")
)
need_cols(hierarchy_table, c(
    "pathway_id", "pathway_name", "hierarchy_status",
    "is_hierarchy_root"
), "supplemental hierarchy")
if (anyDuplicated(hierarchy_table$pathway_id))
    stop("Duplicate hierarchy pathway IDs")

hierarchy_inputs <- file.path(
    hierarchy_dir, c("checksums.sha256", hierarchy_files)
)
# END d2_hierarchy_annotations_v1


pw_membership <- read_gz(b$pathways)
need_cols(
    pw_membership,
    c("pathway_id", "gene_symbol"),
    "pathway discovery membership"
)
pw_membership <- unique(pw_membership[
    !is.na(pathway_id) & !is.na(gene_symbol) &
        nzchar(pathway_id) & nzchar(gene_symbol),
    .(pathway_id, gene_symbol)
])

pw_sizes <- pw_membership[, .(
    pathway_n_genes = uniqueN(gene_symbol)
), by = pathway_id]
pw_membership <- merge(
    pw_membership, pw_sizes,
    by = "pathway_id", all.x = TRUE, sort = FALSE
)

mapped_targets <- unique(dr[
    mapping_status == "mapped" &
        !is.na(target_gene) & nzchar(target_gene),
    .(drug_id, drug_name, target_gene)
])
target_genes <- unique(mapped_targets$target_gene)
pathways_by_gene <- split(
    pw_membership$pathway_id, pw_membership$gene_symbol
)

witness_parts <- vector("list", nrow(h))
for (i in seq_len(nrow(h))) {
    current_gene <- h$gene[i]
    current_pathways <- pathways_by_gene[[current_gene]]
    if (!length(current_pathways)) next

    shared <- pw_membership[
        pathway_id %in% current_pathways &
            gene_symbol %in% target_genes
    ]
    if (!nrow(shared)) next

    w <- merge(
        shared, mapped_targets,
        by.x = "gene_symbol", by.y = "target_gene",
        allow.cartesian = TRUE, sort = FALSE
    )
    setnames(w, "gene_symbol", "drug_target_gene")
    w[, `:=`(
        hypothesis_id = h$hypothesis_id[i],
        gene = current_gene
    )]

    # Same-gene targeting belongs to the existing direct-association branch.
    w <- w[drug_target_gene != gene]
    if (!nrow(w)) next

    witness_parts[[i]] <- unique(w[, .(
        hypothesis_id, gene, drug_id, drug_name,
        drug_target_gene, pathway_id, pathway_n_genes
    )])
}

nonempty <- Filter(Negate(is.null), witness_parts)
if (length(nonempty)) {
    pathway_witnesses <- rbindlist(nonempty, use.names = TRUE)
} else {
    pathway_witnesses <- data.table(
        hypothesis_id = character(),
        gene = character(),
        drug_id = character(),
        drug_name = character(),
        drug_target_gene = character(),
        pathway_id = character(),
        pathway_n_genes = integer()
    )
}


witness_rows_before_hierarchy <- nrow(pathway_witnesses)
pathway_witnesses <- merge(
    pathway_witnesses,
    hierarchy_table[, .(
        pathway_id,
        pathway_name,
        pathway_hierarchy_status = hierarchy_status,
        pathway_is_hierarchy_root = as.integer(is_hierarchy_root)
    )],
    by = "pathway_id", all.x = TRUE, sort = FALSE
)
if (nrow(pathway_witnesses) != witness_rows_before_hierarchy ||
    anyNA(pathway_witnesses$pathway_name) ||
    anyNA(pathway_witnesses$pathway_is_hierarchy_root) ||
    any(!pathway_witnesses$pathway_is_hierarchy_root %in% c(0L, 1L)) ||
    any(pathway_witnesses$pathway_hierarchy_status != "in_hierarchy"))
    stop("Incomplete or inconsistent witness hierarchy annotation")


if (nrow(pathway_witnesses)) {
    pathway_pairs <- pathway_witnesses[, .(
        gene = gene[1L],
        drug_name = paste(sort(unique(drug_name)), collapse = ";"),
        n_shared_target_genes = uniqueN(drug_target_gene),
        n_shared_pathways = uniqueN(pathway_id),
        n_shared_root_pathways = uniqueN(
            pathway_id[pathway_is_hierarchy_root == 1L]
        ),
        n_shared_nonroot_pathways = uniqueN(
            pathway_id[pathway_is_hierarchy_root == 0L]
        ),
        smallest_shared_pathway_n_genes = min(pathway_n_genes),
        largest_shared_pathway_n_genes = max(pathway_n_genes)
    ), by = .(hypothesis_id, drug_id)]
} else {
    pathway_pairs <- data.table(
        hypothesis_id = character(),
        drug_id = character(),
        gene = character(),
        drug_name = character(),
        n_shared_target_genes = integer(),
        n_shared_pathways = integer(),
        n_shared_root_pathways = integer(),
        n_shared_nonroot_pathways = integer(),
        smallest_shared_pathway_n_genes = integer(),
        largest_shared_pathway_n_genes = integer()
    )
}

# A multi-target drug can have both direct and indirect associations.
direct_keys <- unique(merge(
    h[, .(hypothesis_id, gene)],
    mapped_targets[, .(gene = target_gene, drug_id)],
    by = "gene", allow.cartesian = TRUE
)[, .(hypothesis_id, drug_id)])

pathway_pairs[, hierarchy_support := fcase(
    n_shared_nonroot_pathways > 0,
        "non_root_pathway_support",
    n_shared_root_pathways > 0,
        "root_only_pathway_support",
    default = "unclassified"
)]
if (any(pathway_pairs$hierarchy_support == "unclassified"))
    stop("Unclassified pathway pair")

pathway_pairs[, also_has_direct_association := FALSE]
pathway_pairs[
    direct_keys,
    on = .(hypothesis_id, drug_id),
    also_has_direct_association := TRUE
]
pathway_pairs[, `:=`(
    association_type = "indirect_shared_pathway_exploratory",
    direction_status = "not_assessed",
    rescue_status = "not_assessed",
    safety_status = "not_assessed"
)]

pathway_annotation_counts <- pw_membership[, .(
    n_annotated_pathways = uniqueN(pathway_id)
), by = .(gene = gene_symbol)]

pathway_pair_counts <- pathway_pairs[, .(
    n_exploratory_pathway_drugs = uniqueN(drug_id),
    n_root_only_pathway_drugs = sum(
        hierarchy_support == "root_only_pathway_support"
    ),
    n_nonroot_pathway_drugs = sum(
        hierarchy_support == "non_root_pathway_support"
    )
), by = hypothesis_id]

pathway_audit <- merge(
    audit, pathway_annotation_counts,
    by = "gene", all.x = TRUE, sort = FALSE
)
pathway_audit <- merge(
    pathway_audit, pathway_pair_counts,
    by = "hypothesis_id", all.x = TRUE, sort = FALSE
)
pathway_audit[
    is.na(n_annotated_pathways), n_annotated_pathways := 0L
]
pathway_audit[
    is.na(n_exploratory_pathway_drugs),
    n_exploratory_pathway_drugs := 0L
]
for (column in c(
    "n_root_only_pathway_drugs", "n_nonroot_pathway_drugs"
)) {
    set(pathway_audit,
        i = which(is.na(pathway_audit[[column]])),
        j = column, value = 0L)
}
stopifnot(all(
    pathway_audit$n_exploratory_pathway_drugs ==
    pathway_audit$n_root_only_pathway_drugs +
    pathway_audit$n_nonroot_pathway_drugs
))

pathway_audit[, pathway_discovery_status := fcase(
    n_annotated_pathways == 0,
        "no_pathway_annotation",
    n_exploratory_pathway_drugs == 0,
        "annotated_no_mapped_indirect_targets",
    default = "exploratory_shared_pathway_candidates"
)]

stopifnot(
    nrow(pathway_audit) == nrow(h),
    !anyDuplicated(pathway_audit$hypothesis_id),
    !anyDuplicated(pathway_pairs[, .(hypothesis_id, drug_id)])
)
# END d2_pathway_discovery_v1


od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(tu, file.path(od, "target_universe.tsv.gz"))
write_table(audit, file.path(od, "hypothesis_discovery_audit.tsv.gz"))
write_table(pathway_pairs, file.path(od, "pathway_discovery.tsv.gz"))
write_table(pathway_witnesses, file.path(od, "pathway_discovery_witnesses.tsv.gz"))
write_table(pathway_audit, file.path(od, "pathway_discovery_audit.tsv.gz"))

write_node_manifest(
    od, "d2", cfg,
    inputs = c(a$mechanism, a$model_views, b$drugs, b$pathways, hierarchy_inputs),
    outputs = file.path(od, c(
        "target_universe.tsv.gz",
        "hypothesis_discovery_audit.tsv.gz",
        "pathway_discovery.tsv.gz",
        "pathway_discovery_witnesses.tsv.gz",
        "pathway_discovery_audit.tsv.gz"
    )),
    extra = list(
        n_hypotheses = nrow(h),
        n_hypotheses_with_direct_association =
            sum(audit$n_direct_drugs > 0),
        n_hypotheses_without_direct_association =
            sum(audit$n_direct_drugs == 0),
        n_scenarios_with_target = uniqueN(tu$scenario_id),
        n_rows = nrow(tu),
        n_unmapped_targets_dropped = n_un,
        model_views = expected_views,
        hierarchy_snapshot = hierarchy_dir,
        hierarchy_annotations_applied = TRUE,
        hierarchy_numeric_release_verified = isTRUE(
            hierarchy_meta$compatibility_checks$numeric_release_verified
        ),
        n_exploratory_pathway_pairs = nrow(pathway_pairs),
        n_pathway_witnesses = nrow(pathway_witnesses),
        n_hypotheses_with_pathway_annotation =
            sum(pathway_audit$n_annotated_pathways > 0),
        pathway_size_filter_applied = FALSE,
        pathway_hierarchy_filter_applied = FALSE,
        pathway_overlap_is_rescue_evidence = FALSE,
        pathway_pairs_added_to_direct_target_universe = FALSE,
        association_rows_and_original_columns_preserved = TRUE,
        score_threshold_applied = FALSE,
        new_consensus_probability_created = FALSE,
        scalar_weight_scope =
            "unchanged_input_scenario_weights_not_cross_model_consensus"
    )
)
logmsg(
    "D2: hypotheses=", nrow(h),
    "; with direct associations=", sum(audit$n_direct_drugs > 0),
    "; association rows=", nrow(tu),
    "; original association columns and weights preserved"
)
# END d2_model_views_v1
