# workflow/lib/track2_resources.R -- Track 2 resource builders (corrected) and signed propagation checks
suppressPackageStartupMessages({ library(data.table); library(igraph) })

build_gene_protein_map <- function(raw) {                 # returns BOTH tables (secondary UniProt IDs survive)
  x <- as.data.table(copy(raw)); if (!"reviewed" %in% names(x)) x[, reviewed := NA]
  x[, `:=`(uniprot_id = trimws(uniprot_id), gene_symbol = trimws(gene_symbol), reviewed = as.logical(reviewed))]
  x <- unique(x[nzchar(uniprot_id) & nzchar(gene_symbol)])
  lookup <- x[order(uniprot_id, gene_symbol)][, .(gene_symbol = gene_symbol[1L], n_genes = .N), by = uniprot_id]
  canon  <- x[order(gene_symbol, -as.integer(reviewed %in% TRUE), uniprot_id)][, .SD[1L], by = gene_symbol]
  stopifnot(!anyDuplicated(lookup$uniprot_id), !anyDuplicated(canon$gene_symbol))
  list(gene_protein_map = canon[, .(gene_symbol, uniprot_id, reviewed)], uniprot_lookup = lookup)
}

normalize_effect_sign <- function(effect) {
  x <- tolower(trimws(as.character(effect)))
  fcase(x %in% c("up-regulates", "activates", "increases activity", "positive regulation"), 1L,
        x %in% c("down-regulates", "inhibits", "decreases activity", "negative regulation"), -1L, default = NA_integer_)
}
edge_prob <- function(confidence, evidence_count, directness) {
  c0 <- suppressWarnings(as.numeric(confidence)); n0 <- suppressWarnings(as.numeric(evidence_count))
  cpart <- fifelse(is.finite(c0), pmin(pmax(c0, 0), 1), 0.50); npart <- fifelse(is.finite(n0), 1 - exp(-pmax(n0, 0)), 0.50)
  dpart <- fcase(directness == "direct", 1.00, directness == "indirect", 0.75, default = 0.60)
  pmin(pmax(dpart * (0.60 * cpart + 0.40 * npart), 0.05), 0.99)
}

build_signed_edges <- function(raw, lookup, release_id) {
  x <- as.data.table(copy(raw))
  for (col in c("direct", "confidence", "evidence_count")) if (!col %in% names(x)) x[, (col) := NA]
  x <- merge(x, lookup[, .(source_entity = uniprot_id, source_gene = gene_symbol)], by = "source_entity", all.x = TRUE)
  x <- merge(x, lookup[, .(target_entity = uniprot_id, target_gene = gene_symbol)], by = "target_entity", all.x = TRUE)
  x[, sign := normalize_effect_sign(effect)]
  x[, directness := fcase(
    toupper(trimws(as.character(direct))) %in% c("YES", "TRUE", "DIRECT"), "direct",
    toupper(trimws(as.character(direct))) %in% c("NO", "FALSE", "INDIRECT"), "indirect",
    default = "unknown"
  )]
  x[, p := edge_prob(confidence, evidence_count, directness)]
  x[, reason := fcase(is.na(source_gene) | is.na(target_gene), "unmapped_entity", source_gene == target_gene, "self_loop",
                      is.na(sign), "unsigned_effect", default = "kept")]
  e <- x[reason == "kept", .(edge_probability = 1 - prod(1 - p), n_records = .N), by = .(source_gene, target_gene, sign)]
  e[, `:=`(edge_weight = edge_probability, source_release = release_id)]
  list(edges = e[], drop_audit = x[, .N, by = reason])
}
split_edge_conflicts <- function(edges) {
  e <- as.data.table(copy(edges)); f <- e[, .(conflict = any(sign == 1L) & any(sign == -1L)), by = .(source_gene, target_gene)]
  e <- merge(e, f, by = c("source_gene", "target_gene"))
  list(usable_edges = e[conflict == FALSE][, conflict := NULL][], conflict_audit = e[conflict == TRUE])
}
build_signed_graph <- function(edges) {
  e <- as.data.table(copy(edges)); stopifnot(all(e$sign %in% c(-1L, 1L)), all(e$edge_probability > 0 & e$edge_probability <= 1))
  if (anyDuplicated(e[, .(source_gene, target_gene)])) stop("parallel or conflicting edges remain")
  igraph::graph_from_data_frame(e[, .(from = source_gene, to = target_gene, sign, weight = edge_weight, edge_probability)], directed = TRUE)
}

# deterministic signed propagation: x <- tanh(alpha * A^T x + x0), A[u,v] = sign * (weight if expected)
propagate_signed <- function(g, x0, alpha = 0.5, iterations = 50L, use_probability = FALSE, include = NULL) {
  nodes <- V(g)$name; A <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  el <- as_edgelist(g); w <- if (use_probability) E(g)$edge_probability else rep(1, ecount(g))
  if (!is.null(include)) w <- w * include                # Monte Carlo draw: 0/1 edge inclusion
  A[cbind(el[, 1], el[, 2])] <- E(g)$sign * w
  x0v <- stats::setNames(numeric(length(nodes)), nodes); x0v[names(x0)] <- x0; x <- x0v
  for (i in seq_len(iterations)) x <- tanh(alpha * as.numeric(t(A) %*% x) + x0v)
  stats::setNames(x, nodes)
}
sign_convention_audit <- function(g, expected) {          # expected: named signs after perturbing the root
  x <- propagate_signed(g, c(A = 1)); got <- sign(x[names(expected)])
  data.table(node = names(expected), expected = expected, observed = as.integer(got), agree = got == expected)
}
degree_matched_null <- function(g, n_null = 20L, seed) {
  set.seed(seed); lapply(seq_len(n_null), function(i) igraph::rewire(g, igraph::keeping_degseq(niter = 10 * ecount(g))))
}

# ---- resource-bundle builders (Track 2 resource DAG) ----
download_resource <- function(url, output_path, allow_live_downloads = FALSE, overwrite = FALSE, timeout_seconds = 300, expected_sha256 = NULL) {
  if (!isTRUE(allow_live_downloads)) stop("Live downloads disabled; supply a frozen local snapshot.", call. = FALSE)
  if (is.null(url) || !nzchar(url)) stop("Verified URL is not configured.", call. = FALSE)
  if (file.exists(output_path) && !overwrite) return(normalizePath(output_path)); dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  partial <- paste0(output_path, ".partial"); if (file.exists(partial)) unlink(partial); old <- getOption("timeout"); options(timeout = timeout_seconds); on.exit(options(timeout = old), add = TRUE)
  utils::download.file(url, partial, mode = "wb", quiet = FALSE); if (!file.exists(partial) || file.info(partial)$size <= 0) stop("Empty or failed download: ", url)
  if (!file.rename(partial, output_path)) stop("Could not finalize download")
  if (!is.null(expected_sha256) && !identical(sha256_file(output_path), expected_sha256)) stop("Downloaded file checksum mismatch: ", url, call. = FALSE); normalizePath(output_path)
}
validate_registry <- function(reg, root, strict = TRUE) {
  txt <- paste(unlist(reg), collapse = " "); if (grepl("data/private|results/private|C4|hypothesis", txt)) stop("Private/C4 reference in resource registry", call. = FALSE)
  for (id in names(reg$sources)) { s <- reg$sources[[id]]; if (!isTRUE(s$enabled)) next
    if (strict) { bad <- c("url", "version", "license")[vapply(c("url", "version", "license"), function(k) is.null(s[[k]]) || !nzchar(s[[k]]), TRUE)]; if (length(bad)) stop(id, ": registry missing ", paste(bad, collapse = ", "), call. = FALSE) }
    raw <- file.path(root, s$raw_path); if (!file.exists(raw) && !isTRUE(reg$allow_live_downloads) && isTRUE(s$required)) stop(id, ": frozen raw snapshot missing and live downloads disabled: ", raw, call. = FALSE)
    if (file.exists(raw) && !is.null(s$expected_sha256) && !identical(sha256_file(raw), s$expected_sha256)) stop(id, ": raw checksum differs from registry", call. = FALSE) }
  invisible(TRUE)
}
build_reactome_gene_sets <- function(reactome_raw, lookup, release_id) {
  r <- as.data.table(copy(reactome_raw)); need_cols(r, c("uniprot_id", "pathway_id", "pathway_name"), "Reactome raw"); if ("species" %in% names(r)) r <- r[species == "Homo sapiens"]
  out <- merge(r, lookup[, .(uniprot_id, gene_symbol)], by = "uniprot_id"); out[, `:=`(species = "Homo sapiens", source_database = "Reactome", source_release = release_id)]
  out <- unique(out[, .(pathway_id, pathway_name, gene_symbol, uniprot_id, species, source_database, source_release)]); out[, n_genes := uniqueN(gene_symbol), by = pathway_id][]
}
normalize_action <- function(x) { x <- tolower(trimws(as.character(x))); fcase(x %in% c("inhibitor", "antagonist", "degrader"), "inhibitory", x %in% c("agonist", "activator", "stabilizer"), "activating", default = "uncertain") }
build_drug_target_action_map <- function(raw, lookup, release_id) {
  d <- as.data.table(copy(raw)); need_cols(d, c("drug_id", "drug_name", "target_uniprot_id", "action_raw"), "Drug target raw")
  d[, action_normalized := normalize_action(action_raw)][, drug_effect_sign := fcase(action_normalized == "inhibitory", -1L, action_normalized == "activating", 1L, default = NA_integer_)]
  d <- merge(d, lookup[, .(target_uniprot_id = uniprot_id, target_gene = gene_symbol)], by = "target_uniprot_id", all.x = TRUE)
  d[, `:=`(mapping_status = fifelse(is.na(target_gene), "unmapped_target", "mapped"), action_confidence = fifelse(is.na(drug_effect_sign), "uncertain", "explicit"), source_release = release_id)]
  unique(d[, .(drug_id, drug_name, target_uniprot_id, target_gene, action_raw, action_normalized, drug_effect_sign, action_confidence, mapping_status, source_release)])
}
PED_STATUS <- c("supportive", "limited", "unknown", "contraindicated"); CONTRA <- c("none", "relative", "absolute", "unknown")
build_pediatric_safety_map <- function(raw, release_id) {
  s <- as.data.table(copy(raw)); need_cols(s, c("drug_id", "pediatric_evidence_status", "contraindication_status", "route_feasible"), "Safety raw")
  s[, `:=`(pediatric_evidence_status = tolower(trimws(pediatric_evidence_status)), contraindication_status = tolower(trimws(contraindication_status)), route_feasible = as.logical(route_feasible))]
  stopifnot(all(s$pediatric_evidence_status %in% PED_STATUS), all(s$contraindication_status %in% CONTRA))
  out <- s[order(-(pediatric_evidence_status == "contraindicated"), -(contraindication_status == "absolute"))][, .SD[1L], by = drug_id]      # most restrictive row per drug
  out[, source_release := release_id][]
}
audit_signed_graph <- function(g, edges, drugs = NULL) {
  e <- as.data.table(edges); sc <- igraph::components(g, mode = "strong")
  data.table(n_nodes = igraph::vcount(g), n_edges = igraph::ecount(g), directed = igraph::is_directed(g), positive_fraction = mean(e$sign == 1L), n_strong_components = sc$no,
             largest_scc_fraction = max(sc$csize) / igraph::vcount(g), zero_indegree_fraction = mean(igraph::degree(g, mode = "in") == 0), median_edge_probability = median(e$edge_probability),
             drug_targets_in_graph = if (is.null(drugs)) NA_real_ else mean(unique(drugs$target_gene[drugs$mapping_status == "mapped"]) %in% igraph::V(g)$name))
}


# BEGIN integrated_ChEMBL_preparation_v1
prepare_chembl_drug_actions <- function(native_dir, out) {
    mf_path <- file.path(native_dir, "download_manifest.json")
    sums_path <- file.path(native_dir, "native.sha256")
    mf <- jsonlite::read_json(mf_path)

    if (!identical(mf$chembl_db_version, "ChEMBL_37"))
        stop("Unexpected ChEMBL version")
    if (!identical(sha256_file(sums_path),
                   mf$native_sha256_file_sha256))
        stop("Native checksum-list hash mismatch")

    pairs <- strsplit(readLines(sums_path, warn = FALSE),
                      "  ", fixed = TRUE)
    checked <- character()

    for (pair in pairs) {
        if (length(pair) != 2L) stop("Malformed native checksum entry")
        path <- normalizePath(file.path(native_dir, pair[2L]),
                              mustWork = TRUE)
        prefix <- paste0(normalizePath(native_dir), "/")
        if (!startsWith(path, prefix))
            stop("Native checksum entry escapes snapshot directory")
        if (!identical(sha256_file(path), pair[1L]))
            stop("Native checksum mismatch: ", pair[2L])
        checked <- c(checked, path)
    }

    read_collection <- function(pattern, key, id_key) {
        files <- sort(list.files(
            file.path(native_dir, "pages"),
            pattern = pattern, full.names = TRUE
        ))
        if (!length(files)) stop("Missing collection: ", key)

        records <- list()
        for (path in files) {
            if (!normalizePath(path) %in% checked)
                stop("Unmanifested input page: ", path)
            page <- jsonlite::read_json(path)
            if (!is.list(page[[key]]))
                stop("Unexpected collection schema: ", key)
            records <- c(records, page[[key]])
        }

        ids <- vapply(records, function(x) {
            value <- x[[id_key]]
            if (is.null(value)) "" else as.character(value)
        }, character(1))

        if (any(!nzchar(ids)) || anyDuplicated(ids))
            stop("Missing or duplicate collection IDs: ", key)

        stats::setNames(records, ids)
    }

    mechanisms <- read_collection(
        "^mechanism_[0-9]+[.]json$", "mechanisms", "mec_id"
    )
    targets <- read_collection(
        "^target_batch[0-9]+_[0-9]+[.]json$",
        "targets", "target_chembl_id"
    )
    molecules <- read_collection(
        "^molecule_batch[0-9]+_[0-9]+[.]json$",
        "molecules", "molecule_chembl_id"
    )

    if (length(mechanisms) != mf$counts$mechanisms ||
        length(targets) != mf$counts$referenced_targets ||
        length(molecules) != mf$counts$referenced_molecules_and_parents)
        stop("Acquisition record-count mismatch")

    txt <- function(x) {
        if (is.null(x) || !length(x)) "" else as.character(x)[1L]
    }
    yes <- function(x) {
        txt(x) %in% c("1", "TRUE", "true")
    }
    no <- function(x) {
        txt(x) %in% c("0", "FALSE", "false")
    }

    admitted <- list()
    audit <- list()
    allowed <- c("INHIBITOR", "ANTAGONIST", "ACTIVATOR", "AGONIST")

    for (id in names(mechanisms)) {
        m <- mechanisms[[id]]
        tid <- txt(m$target_chembl_id)
        did <- txt(m$molecule_chembl_id)

        t <- if (nzchar(tid)) targets[[tid]] else NULL
        drug <- if (nzchar(did)) molecules[[did]] else NULL

        components <- t$target_components
        if (is.null(components)) components <- list()
        if (!is.list(components))
            stop("Malformed target components: ", tid)

        accession <- if (length(components) == 1L)
            txt(components[[1L]]$accession) else ""
        component_type <- if (length(components) == 1L)
            txt(components[[1L]]$component_type) else ""

        action <- toupper(trimws(txt(m$action_type)))
        reasons <- character()

        if (!nzchar(tid) || is.null(t))
            reasons <- c(reasons, "missing_target")
        if (txt(t$tax_id) != "9606")
            reasons <- c(reasons, "not_explicit_human_target")
        if (!no(t$species_group_flag))
            reasons <- c(reasons, "species_group_or_unknown")
        if (txt(t$target_type) != "SINGLE PROTEIN")
            reasons <- c(reasons, "not_single_protein")
        if (length(components) != 1L ||
            component_type != "PROTEIN" || !nzchar(accession))
            reasons <- c(reasons, "not_one_explicit_protein_accession")
        if (!yes(m$direct_interaction))
            reasons <- c(reasons, "direct_interaction_not_confirmed")
        if (!yes(m$molecular_mechanism))
            reasons <- c(reasons, "molecular_mechanism_not_confirmed")
        if (!action %in% allowed)
            reasons <- c(reasons, "action_outside_signed_policy")
        if (!is.null(m$variant_sequence) &&
            length(m$variant_sequence) &&
            !identical(m$variant_sequence, ""))
            reasons <- c(reasons, "variant_specific_mechanism")
        if (!nzchar(did) || is.null(drug))
            reasons <- c(reasons, "missing_molecule")

        audit[[length(audit) + 1L]] <- data.table::data.table(
            mechanism_id = id,
            drug_id = did,
            parent_drug_id = txt(m$parent_molecule_chembl_id),
            target_chembl_id = tid,
            target_type = txt(t$target_type),
            target_tax_id = txt(t$tax_id),
            action_raw = txt(m$action_type),
            target_uniprot_id = accession,
            status = if (length(reasons)) "excluded" else "admitted",
            reasons = paste(reasons, collapse = ";"),
            mechanism_refs_json = as.character(jsonlite::toJSON(
                m$mechanism_refs, auto_unbox = TRUE, null = "null"
            ))
        )

        if (!length(reasons)) {
            admitted[[length(admitted) + 1L]] <-
                data.table::data.table(
                    drug_id = did,
                    drug_name = txt(drug$pref_name),
                    target_uniprot_id = accession,
                    action_raw = action,
                    mechanism_id = id,
                    parent_drug_id = txt(m$parent_molecule_chembl_id),
                    target_chembl_id = tid
                )
        }
    }

    raw <- if (length(admitted)) {
        data.table::rbindlist(admitted)
    } else {
        data.table::data.table(
            drug_id = character(), drug_name = character(),
            target_uniprot_id = character(), action_raw = character(),
            mechanism_id = character(), parent_drug_id = character(),
            target_chembl_id = character()
        )
    }

    au <- data.table::rbindlist(audit, fill = TRUE)
    audit_path <- paste0(out, ".preparation_audit.tsv.gz")
    summary_path <- paste0(out, ".preparation_summary.json")

    dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
    write_table(raw, out)
    write_table(au, audit_path)

    release <- paste0(
        mf$chembl_db_version, "; release ", mf$chembl_release_date
    )
    notes <- paste0(
        "Prepared offline from native acquisition manifest sha256=",
        sha256_file(mf_path),
        "; exact component accessions; human single-protein targets; ",
        "direct molecular mechanisms; four-action policy; ",
        "variant-specific records excluded. ",
        "Rejected records retained in preparation audit. ",
        "No gene mapping, approval filtering, safety or efficacy inference."
    )

    write_resource_manifest(
        out, "drug_targets_raw", mf$api_base, release,
        mf$license, notes = notes
    )

    summary <- list(
        source_release = release,
        native_manifest_sha256 = sha256_file(mf_path),
        native_checksums_sha256 = sha256_file(sums_path),
        raw_sha256 = sha256_file(out),
        audit_sha256 = sha256_file(audit_path),
        counts = list(
            native_mechanisms = length(mechanisms),
            admitted_mechanisms = nrow(raw),
            excluded_mechanisms = nrow(au) - nrow(raw),
            distinct_drugs = data.table::uniqueN(raw$drug_id),
            distinct_accessions = data.table::uniqueN(
                raw$target_uniprot_id
            )
        ),
        exclusion_reasons = au[status == "excluded", .N, by = reasons]
    )
    jsonlite::write_json(
        summary, summary_path, auto_unbox = TRUE,
        pretty = TRUE, dataframe = "rows"
    )

    cat("OFFLINE CHEMBL PREPARATION COMPLETE\n")
    print(summary$counts)
    print(summary$exclusion_reasons)
    if (!nrow(raw))
        warning("No mechanisms admitted; inspect component schema and policy")
    invisible(summary)
}
# END integrated_ChEMBL_preparation_v1
