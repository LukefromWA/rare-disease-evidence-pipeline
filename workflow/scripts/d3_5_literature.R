#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))

# BEGIN mechanism_literature_search_v1
# Literature retrieval only: no efficacy, safety or causal probabilities.

a <- parse_args(c("config", "sample", "scenarios", "outdir"))
cfg <- read_config(a$config)
sc <- read_gz(a$scenarios)

need_cols(sc, c(
    "hypothesis_id", "gene", "variant", "consequence"
), "mechanism-search hypotheses")

h <- unique(sc[, .(
    hypothesis_id, gene, variant, consequence
)])
if (anyDuplicated(h$hypothesis_id))
    stop("Inconsistent hypothesis annotations")

h[, stop_gained := !is.na(consequence) &
    grepl("stop_gained", consequence, fixed = TRUE)]

settings <- cfg$track2$mechanism_search
enabled <- isTRUE(settings$enabled)
allow_live <- isTRUE(settings$allow_live_downloads) &&
    !identical(Sys.getenv("MVA_MECHANISM_OFFLINE"), "1")
limit <- settings$max_results_per_query %||% 50L
if (length(limit) != 1L || !is.numeric(limit) ||
    !is.finite(limit) || limit < 1 || limit > 1000 ||
    limit != as.integer(limit))
    stop("max_results_per_query must be an integer from 1 to 1000")
limit <- as.integer(limit)

queries <- data.table(
    query_id = c(
        "ptc_readthrough",
        "nmd_transcript_availability",
        "ribosome_translation_modulation"
    ),
    query = c(
        paste(
            '("premature termination codon" OR "nonsense mutation")',
            'AND (readthrough OR "read-through")',
            'AND (drug OR compound OR pharmacological)'
        ),
        paste(
            '("nonsense-mediated decay" OR "nonsense mediated mRNA decay")',
            'AND (inhibitor OR inhibition OR modulator OR compound)',
            'AND (nonsense OR "premature termination")'
        ),
        paste(
            '("ribosome stabilization" OR "ribosome rescue"',
            'OR "ribosome stalling" OR "translation termination")',
            'AND (drug OR compound OR pharmacological)',
            'AND (nonsense OR "premature termination")'
        )
    )
)

od <- a$outdir
dir.create(od, recursive = TRUE, showWarnings = FALSE)
cache_dir <- file.path(od, "mechanism_search_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

fetch_query <- function(query) {
    key <- digest::digest(
        paste(query, limit, "core-json-v1", sep = "\n"),
        algo = "sha256", serialize = FALSE
    )
    path <- file.path(cache_dir, paste0(key, ".json"))

    if (!file.exists(path)) {
        if (!allow_live)
            stop("Uncached query and live mechanism search is disabled")

        url <- paste0(
            "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            "?query=", utils::URLencode(query, reserved = TRUE),
            "&format=json&resultType=core&pageSize=", limit
        )
        partial <- tempfile(pattern = ".query_", tmpdir = cache_dir)
        on.exit(unlink(partial), add = TRUE)

        old_timeout <- getOption("timeout")
        options(timeout = 90)
        on.exit(options(timeout = old_timeout), add = TRUE)

        Sys.sleep(0.7)
        utils::download.file(
            url, partial, mode = "wb", quiet = TRUE
        )
        candidate <- jsonlite::read_json(partial)
        if (is.null(candidate$hitCount) ||
            !is.list(candidate$resultList))
            stop("Unexpected Europe PMC response; no evidence invented")

        if (!file.rename(partial, path))
            stop("Could not finalize query cache")
    }

    response <- jsonlite::read_json(path)
    if (is.null(response$hitCount) ||
        !is.list(response$resultList))
        stop("Invalid cached Europe PMC response")

    list(response = response, path = path)
}

txt <- function(x) {
    if (is.null(x) || !length(x)) "" else as.character(x)[1L]
}

publications <- data.table(
    query_id = character(),
    record_source = character(),
    record_id = character(),
    pmid = character(),
    doi = character(),
    title = character(),
    abstract_html = character(),
    publication_year = character(),
    article_url = character()
)
query_audit <- copy(queries)
query_audit[, `:=`(
    returned_records = 0L,
    reported_total_hits = NA_real_,
    truncated = NA,
    cache_sha256 = NA_character_,
    retrieval_status = "not_run"
)]
cache_inputs <- character()

if (!enabled) {
    state <- "disabled_pending_enablement"
    reason <- "Mechanism search disabled in configuration"
} else if (!any(h$stop_gained)) {
    state <- "not_applicable"
    reason <- "Initial module targets stop_gained hypotheses; none present"
} else {
    for (i in seq_len(nrow(queries))) {
        fetched <- fetch_query(queries$query[i])
        records <- fetched$response$resultList$result
        if (is.null(records)) records <- list()
        total <- as.numeric(fetched$response$hitCount)

        if (length(total) != 1L || !is.finite(total) || total < 0)
            stop("Invalid Europe PMC hit count")

        parts <- lapply(records, function(record) {
            source <- txt(record$source)
            id <- txt(record$id)
            if (!nzchar(source) || !nzchar(id))
                stop("Publication lacks a stable source/id")

            data.table(
                query_id = queries$query_id[i],
                record_source = source,
                record_id = id,
                pmid = if (source == "MED") id else txt(record$pmid),
                doi = txt(record$doi),
                title = txt(record$title),
                abstract_html = txt(record$abstractText),
                publication_year = txt(record$pubYear),
                article_url = paste0(
                    "https://europepmc.org/article/", source, "/", id
                )
            )
        })

        if (length(parts)) {
            batch <- rbindlist(parts)
            if (anyDuplicated(batch[, .(record_source, record_id)]))
                stop("Duplicate publication IDs within a query response")
            publications <- rbindlist(list(publications, batch))
        }

        query_audit[i, `:=`(
            returned_records = length(records),
            reported_total_hits = total,
            truncated = total > length(records),
            cache_sha256 = sha256_file(fetched$path),
            retrieval_status = "retrieved"
        )]
        cache_inputs <- c(cache_inputs, fetched$path)
    }

    state <- "ok"
    reason <- paste(
        "Mechanism literature retrieved for screening;",
        "bounded search; unreviewed seeded abstract extraction, no rescue inference"
    )
}

links <- data.table(
    hypothesis_id = character(),
    query_id = character(),
    applicability_status = character()
)
if (enabled && any(h$stop_gained)) {
    links <- data.table::CJ(
        hypothesis_id = h[stop_gained == TRUE, hypothesis_id],
        query_id = queries$query_id,
        unique = TRUE
    )
    links[, applicability_status :=
        "screening_rationale_only_variant_applicability_unreviewed"]
}

hypothesis_audit <- copy(h)
hypothesis_audit[, search_scope := fifelse(
    stop_gained,
    "stop_gained_mechanism_screen",
    "not_targeted_by_initial_module"
)]
hypothesis_audit[, search_enabled := enabled]

worksheet <- copy(publications[, .(
    query_id, record_source, record_id, pmid, doi,
    title, article_url
)])
worksheet[, `:=`(
    compound_name = "",
    experimental_model = "",
    reported_effect = "",
    evidence_locator = "",
    gene_specific_evidence = "not_assessed",
    variant_specific_evidence = "not_assessed",
    readthrough_or_functional_rescue = "not_assessed",
    safety_status = "not_assessed",
    review_state = "unreviewed",
    reviewer_notes = ""
)]


# BEGIN compound_evidence_ledger_v1
registry_path <- settings$compound_seed_registry
if (is.null(registry_path) || !nzchar(registry_path))
    stop("Configure mechanism_search.compound_seed_registry")
registry_manifest_path <- paste0(registry_path, ".manifest.json")
registry_meta <- jsonlite::read_json(registry_manifest_path)
if (!identical(sha256_file(registry_path), registry_meta$sha256))
    stop("Compound seed registry checksum mismatch")

# BEGIN seed_registry_csv_reader_v1
read_seed_registry <- function(path) {
    data.table::as.data.table(utils::read.delim(
        file = path,
        header = TRUE,
        sep = "\t",
        quote = "\"",
        comment.char = "",
        colClasses = "character",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        na.strings = character(),
        fill = FALSE,
        allowEscapes = FALSE
    ))
}
seed_registry <- read_seed_registry(registry_path)
# END seed_registry_csv_reader_v1
need_cols(seed_registry, c(
    "assertion_id", "compound_name", "aliases_json", "pmid",
    "reported_role", "experimental_gene", "experimental_model",
    "reported_outcome", "curation_state"
), "compound seed registry")
if (anyDuplicated(seed_registry$assertion_id))
    stop("Duplicate seed assertion IDs")
if (any(seed_registry$curation_state != "unreviewed_abstract_seed"))
    stop("Unexpected seed curation state")

clean_abstract <- function(text) {
    text <- gsub("&lt;", "<", text, fixed = TRUE)
    text <- gsub("&gt;", ">", text, fixed = TRUE)
    text <- gsub("<[^>]+>", " ", text)
    text <- gsub("&amp;", "&", text, fixed = TRUE)
    trimws(gsub("[[:space:]]+", " ", text))
}

# Literal matching with alphanumeric boundaries, restricted to the seed PMID.
find_alias <- function(aliases, text) {
    lower_text <- tolower(text)
    for (alias in aliases) {
        positions <- gregexpr(
            tolower(alias), lower_text, fixed = TRUE
        )[[1L]]
        if (positions[1L] < 0L) next

        for (position in positions) {
            last <- position + nchar(alias) - 1L
            before <- if (position > 1L)
                substr(lower_text, position - 1L, position - 1L) else ""
            after <- if (last < nchar(lower_text))
                substr(lower_text, last + 1L, last + 1L) else ""

            if (!grepl("[[:alnum:]]", before) &&
                !grepl("[[:alnum:]]", after))
                return(alias)
        }
    }
    NA_character_
}

compound_evidence <- data.table(
    assertion_id = character(), compound_name = character(),
    matched_alias = character(), publication_id = character(),
    publication_records = character(), query_ids = character(),
    reported_role = character(), experimental_gene = character(),
    experimental_model = character(), reported_outcome = character(),
    evidence_location = character(), supporting_quote = character(),
    raw_abstract_sha256 = character(),
    source_abstract_versions = integer(),
    applicability_to_current_hypothesis = character(),
    intervention_direction = character(),
    safety_status = character(), review_state = character(),
    annotation_origin = character()
)

extraction_audit <- seed_registry[, .(
    assertion_id, compound_name, pmid
)]
extraction_audit[, extraction_status := "source_publication_not_in_cache"]

for (i in seq_len(nrow(seed_registry))) {
    seed <- seed_registry[i]
    aliases <- jsonlite::fromJSON(seed$aliases_json)
    if (!is.character(aliases) || !length(aliases) ||
        anyNA(aliases) || any(!nzchar(aliases)))
        stop("Invalid seed alias list")

    source_rows <- publications[
        pmid == seed$pmid |
            (record_source == "MED" & record_id == seed$pmid)
    ]
    if (!nrow(source_rows)) next

    source_rows <- copy(source_rows)
    source_rows[, clean_text := vapply(
        abstract_html, clean_abstract, character(1)
    )]
    source_rows[, matched_alias := vapply(
        clean_text, function(text) find_alias(aliases, text), character(1)
    )]
    matched <- source_rows[!is.na(matched_alias)]

    if (!nrow(matched)) {
        extraction_audit[i, extraction_status :=
            "publication_present_alias_not_found_in_abstract"]
        next
    }

    matched[, abstract_length := nchar(clean_text)]
    setorderv(
        matched,
        c("abstract_length", "record_source", "record_id", "query_id"),
        c(-1L, 1L, 1L, 1L)
    )
    chosen <- matched[1L]

    record <- data.table(
        assertion_id = seed$assertion_id,
        compound_name = seed$compound_name,
        matched_alias = chosen$matched_alias,
        publication_id = paste0("PMID:", seed$pmid),
        publication_records = paste(sort(unique(paste(
            source_rows$record_source, source_rows$record_id, sep = ":"
        ))), collapse = ";"),
        query_ids = paste(
            sort(unique(source_rows$query_id)), collapse = ";"
        ),
        reported_role = seed$reported_role,
        experimental_gene = seed$experimental_gene,
        experimental_model = seed$experimental_model,
        reported_outcome = seed$reported_outcome,
        evidence_location = "cached_abstract_basic_HTML_cleanup",
        supporting_quote = chosen$clean_text,
        raw_abstract_sha256 = digest::digest(
            chosen$abstract_html, algo = "sha256", serialize = FALSE
        ),
        source_abstract_versions = uniqueN(source_rows$abstract_html),
        applicability_to_current_hypothesis = "not_established",
        intervention_direction = "not_assessed",
        safety_status = "not_assessed",
        review_state = "unreviewed",
        annotation_origin = "seeded_abstract_screen_not_full_text_review"
    )
    compound_evidence <- rbindlist(list(compound_evidence, record))
    extraction_audit[i, extraction_status := "abstract_mention_matched"]
}

# Indexed subjects remain a separate table: not assumed to be medications.
chemical_mentions <- data.table(
    record_source = character(), record_id = character(),
    indexed_name = character(), registry_number_raw = character(),
    source_cache_sha256 = character()
)

for (path in unique(cache_inputs)) {
    raw <- jsonlite::read_json(path)
    records <- raw$resultList$result
    if (is.null(records)) records <- list()

    for (record in records) {
        chemicals <- record$chemicalList$chemical
        if (is.null(chemicals) || !length(chemicals)) next
        if (!is.list(chemicals))
            stop("Unexpected indexed chemical schema")
        if (!is.null(chemicals$name)) chemicals <- list(chemicals)

        for (chemical in chemicals) {
            if (!is.list(chemical))
                stop("Unexpected indexed chemical entry")
            name <- txt(chemical$name)
            if (!nzchar(name)) next

            chemical_mentions <- rbindlist(list(
                chemical_mentions,
                data.table(
                    record_source = txt(record$source),
                    record_id = txt(record$id),
                    indexed_name = name,
                    registry_number_raw = txt(chemical$registryNumber),
                    source_cache_sha256 = sha256_file(path)
                )
            ))
        }
    }
}
chemical_mentions <- unique(chemical_mentions)
chemical_mentions[, `:=`(
    entity_type = "unclassified_indexed_subject",
    intervention_evidence = "not_established",
    review_state = "unreviewed"
)]

stopifnot(!anyDuplicated(compound_evidence$assertion_id))
# END compound_evidence_ledger_v1


files <- c(
    "d3_5_literature_status.tsv",
    "mechanism_search_queries.tsv",
    "mechanism_publications.tsv.gz",
    "hypothesis_mechanism_search_links.tsv.gz",
    "mechanism_search_hypothesis_audit.tsv.gz",
    "mechanism_compound_review_worksheet.tsv",
    "mechanism_compound_evidence.tsv.gz",
    "mechanism_chemical_index_mentions.tsv.gz",
    "mechanism_compound_extraction_audit.tsv"
)
write_table(data.table(
    method = "d3_5_literature",
    status = state,
    reason = reason
), file.path(od, files[1L]))
write_table(query_audit, file.path(od, files[2L]))
write_table(publications, file.path(od, files[3L]))
write_table(links, file.path(od, files[4L]))
write_table(hypothesis_audit, file.path(od, files[5L]))
write_table(worksheet, file.path(od, files[6L]))
write_table(compound_evidence, file.path(od, "mechanism_compound_evidence.tsv.gz"))
write_table(chemical_mentions, file.path(od, "mechanism_chemical_index_mentions.tsv.gz"))
write_table(extraction_audit, file.path(od, "mechanism_compound_extraction_audit.tsv"))

write_node_manifest(
    od, "d3_5", cfg,
    inputs = c(a$scenarios, unique(cache_inputs), registry_path, registry_manifest_path),
    outputs = file.path(od, files),
    file = "manifest_d3_5_literature.json",
    extra = list(
        status = state,
        provider = "Europe PMC",
        n_stop_gained_hypotheses = sum(h$stop_gained),
        max_results_per_query = limit,
        any_query_truncated = any(query_audit$truncated %in% TRUE),
        private_case_fields_sent = FALSE,
        automatic_compound_extraction = FALSE,
        seeded_abstract_assertions_emitted = TRUE,
        seed_registry_exhaustive = FALSE,
        indexed_subjects_are_validated_drugs = FALSE,
        live_access_allowed = allow_live,
        automatic_efficacy_score = FALSE,
        mechanism_priors_updated = FALSE,
        safety_review_performed = FALSE
    )
)
logmsg("D3.5: ", state, "; publication-query records=", nrow(publications))
# END mechanism_literature_search_v1
