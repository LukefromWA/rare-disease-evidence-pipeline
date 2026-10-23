#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

a <- parse_args(
    c("ready", "manifest", "checksums", "root", "graph", "edges",
      "conflicts", "audit", "genes", "lookup", "pathways", "drugs"),
    c("safety", "scope")
)

scope <- a$scope %||% "full"
if (!scope %in% c("core", "full")) stop("Invalid scope")

# A failed revalidation must not leave old readiness behind.
if (file.exists(a$ready)) unlink(a$ready)

P <- a[c(
    "graph", "edges", "conflicts", "audit",
    "genes", "lookup", "pathways", "drugs"
)]

if (scope == "full") {
    if (is.null(a$safety)) stop("Full scope requires --safety")
    P$safety <- a$safety
}

for (id in names(P)) assert_frozen_resource(P[[id]], id)

edge <- read_gz(P$edges)
gene <- read_gz(P$genes)
lookup <- read_gz(P$lookup)
drug <- read_gz(P$drugs)
graph <- readRDS(P$graph)

need_cols(edge, c(
    "source_gene", "target_gene", "sign", "edge_probability"
), "signed edges")
need_cols(gene, "gene_symbol", "gene map")
need_cols(lookup, "uniprot_id", "UniProt lookup")
need_cols(drug, c(
    "drug_id", "target_uniprot_id", "action_normalized",
    "mapping_status"
), "drug map")

stopifnot(
    all(edge$sign %in% c(-1L, 1L)),
    all(is.finite(edge$edge_probability)),
    all(edge$edge_probability > 0 & edge$edge_probability <= 1),
    !anyDuplicated(edge[, .(source_gene, target_gene)]),
    !anyDuplicated(gene$gene_symbol),
    !anyDuplicated(lookup$uniprot_id),
    igraph::is_directed(graph),
    igraph::ecount(graph) == nrow(edge)
)

if (scope == "full") {
    safety <- read_gz(P$safety)
    need_cols(safety, c(
        "drug_id", "pediatric_evidence_status",
        "contraindication_status", "route_feasible"
    ), "safety map")
    stopifnot(
        !anyDuplicated(safety$drug_id),
        all(safety$pediatric_evidence_status %in% PED_STATUS),
        all(safety$contraindication_status %in% CONTRA)
    )
}

manifest_text <- unlist(lapply(P, function(path) {
    readLines(paste0(path, ".resource_manifest.json"), warn = FALSE)
}))

if (any(grepl("data/private|results/private|C4", manifest_text)))
    stop("Private/C4 reference in public resource manifests")

dir.create(dirname(a$manifest), recursive = TRUE, showWarnings = FALSE)
write_bundle_checksums(unlist(P), a$checksums, root = a$root)

jsonlite::write_json(
    list(
        scope = scope,
        safety_scope = if (scope == "core")
            "case_level_not_validated_here" else "bundle_level",
        artifacts = lapply(P, function(path) {
            list(path = normalizePath(path), sha256 = sha256_file(path))
        }),
        timestamp_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    a$manifest, auto_unbox = TRUE, pretty = TRUE
)

dir.create(dirname(a$ready), recursive = TRUE, showWarnings = FALSE)
writeLines(paste0("validated_scope=", scope), a$ready)
logmsg("Validated resource scope: ", scope)
