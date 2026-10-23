#!/usr/bin/env Rscript
.sf <- gsub("~+~", " ", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]), fixed=TRUE)
source(file.path(dirname(normalizePath(.sf, mustWork=TRUE)), "..", "lib", "bootstrap.R"))
a <- parse_args(c("native", "native_sha256", "lookup", "outdir", "source_url", "source_version", "license"))
if (file.exists(a$outdir)) stop("Refusing to reuse output directory: ",a$outdir)
if (!identical(sha256_file(a$native),tolower(trimws(a$native_sha256)))) stop("Native Reactome checksum mismatch")
assert_frozen_resource(a$lookup,"uniprot_gene_lookup")
lib <- file.path(dirname(normalizePath(.sf)),"..","lib","track2_resources.R")
helper_sha <- sha256_file(lib)
lookup_sha <- sha256_file(a$lookup)
script_sha <- sha256_file(.sf)

# Synthetic regression fixture: repeated accession and two accessions per gene.
fixture <- data.table::data.table(
  uniprot_id=c("fixture_U1","fixture_U1","fixture_U2","fixture_U3","fixture_MOUSE"),
  pathway_id=c(rep("R-HSA-FIXTURE",4),"R-MMU-FIXTURE"),
  pathway_name=c(rep("Fixture human pathway",4),"Fixture mouse pathway"),
  species=c(rep("Homo sapiens",4),"Mus musculus")
)
fixture_lookup <- data.table::data.table(uniprot_id=c("fixture_U1","fixture_U2","fixture_U3","fixture_MOUSE"),gene_symbol=c("fixture_G1","fixture_G1","fixture_G2","fixture_G3"))
ft <- build_reactome_gene_sets(fixture,fixture_lookup,"SYNTHETIC_REGRESSION_ONLY")
if (nrow(ft)!=3L || !all(ft$species=="Homo sapiens") || !all(ft$n_genes==2L)) stop("Reactome helper regression failed: n_genes must count DISTINCT genes; duplicate mappings and nonhuman rows must not inflate it. No real outputs written.")
cat("Human filtering, duplicate handling and distinct-gene regression passed.\n")

raw <- data.table::fread(a$native,header=FALSE,sep="\t",quote="",colClasses="character",na.strings=NULL,fill=FALSE)
if (ncol(raw)!=6L) stop("Expected six native Reactome columns; found ",ncol(raw))
data.table::setnames(raw,c("uniprot_id","pathway_id","url","pathway_name","evidence_code","species"))
for (nm in c("uniprot_id","pathway_id","species")) data.table::set(raw,j=nm,value=trimws(raw[[nm]]))
raw[, native_source_row := .I]
human <- data.table::copy(raw[species=="Homo sapiens"])
if (!nrow(human) || any(!grepl("^R-HSA-",human$pathway_id))) stop("Missing human mappings or unexpected human pathway identifiers")
lookup <- read_gz(a$lookup,colClasses="character",na.strings=NULL)
need_cols(lookup,c("uniprot_id","gene_symbol"),"UniProt lookup")
if (anyDuplicated(lookup$uniprot_id) || any(!nzchar(lookup$uniprot_id)) || any(!nzchar(lookup$gene_symbol))) stop("Nonunique or blank lookup keys")
mi <- match(human$uniprot_id,lookup$uniprot_id)
data.table::set(human,j="gene_symbol",value=lookup$gene_symbol[mi])
human[, mapping_status := ifelse(is.na(gene_symbol),"unmapped_accession","mapped")]
unmapped <- human[mapping_status=="unmapped_accession", .(n_native_membership_rows=.N),by=uniprot_id]

# Reuse the canonical scientific helper; retain annotation-level evidence separately.
result <- build_reactome_gene_sets(raw,lookup,a$source_version)
if (!nrow(result)) stop("No mapped pathway memberships")
need_cols(result,c("pathway_id","gene_symbol","uniprot_id","species","n_genes"),"Reactome output")
if (anyDuplicated(result[, .(pathway_id,gene_symbol,uniprot_id)])) stop("Duplicate pathway/gene/accession memberships")
sizes <- result[, .(actual_distinct_genes=data.table::uniqueN(gene_symbol), n_reported_values=data.table::uniqueN(n_genes), reported_n_genes=n_genes[1L]),by=pathway_id]
if (any(sizes$n_reported_values!=1L) || any(sizes$actual_distinct_genes!=sizes$reported_n_genes)) stop("Incorrect pathway gene counts")
if (!identical(helper_sha,sha256_file(lib)) || !identical(lookup_sha,sha256_file(a$lookup)) || !identical(tolower(trimws(a$native_sha256)),sha256_file(a$native))) stop("An input changed during processing")

dir.create(a$outdir,recursive=TRUE,showWarnings=FALSE)
paths <- c(gene_sets=file.path(a$outdir,"reactome_gene_sets.tsv.gz"),
           membership_evidence=file.path(a$outdir,"human_membership_evidence.tsv.gz"),
           unmapped=file.path(a$outdir,"unmapped_human_accessions.tsv"),
           pathway_sizes=file.path(a$outdir,"pathway_gene_counts.tsv"))
write_table(result,paths[["gene_sets"]]);write_table(human,paths[["membership_evidence"]]);write_table(unmapped,paths[["unmapped"]]);write_table(sizes,paths[["pathway_sizes"]])
chk <- read_gz(paths[["gene_sets"]],colClasses="character")
if (nrow(chk)!=nrow(result) || anyNA(suppressWarnings(as.integer(chk$n_genes)))) stop("Gene-set read-back check failed")
notes <- paste0("UniProt2Reactome_All_Levels native sha256=",a$native_sha256,"; human-only memberships; exact accession lookup; n_genes counts distinct mapped genes, not rows/accessions. Full evidence codes and URLs retained in the separate human_membership_evidence table. All-level hierarchical memberships are overlapping annotations, not independent evidence. No pathway-size exclusion, enrichment test, causal sign inference or case-specific selection performed. Helper sha256=",helper_sha,"; adapter sha256=",script_sha,".")
for (nm in names(paths)) {
  write_resource_manifest(paths[[nm]],paste0("reactome_",nm),a$source_url,a$source_version,a$license,notes=notes)
  assert_frozen_resource(paths[[nm]],nm)
}
native_mf <- paste0(a$native,".resource_manifest.json")
if (file.exists(native_mf)) {
  assert_frozen_resource(a$native,"reactome_native")
} else {
  write_resource_manifest(a$native,"reactome_all_levels_native",a$source_url,a$source_version,a$license,
    notes="Unmodified six-column headerless multi-species UniProt2Reactome_All_Levels export. Retrieval label is not a claimed release number. Derived transformations are separately manifested.")
}
summary <- list(adapter="e04_build_reactome_scratch_v1",timestamp_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
  native=list(path=normalizePath(a$native),sha256=sha256_file(a$native),rows=nrow(raw)),
  lookup=list(path=normalizePath(a$lookup),sha256=lookup_sha),helper_sha256=helper_sha,adapter_sha256=script_sha,
  counts=list(human_membership_rows=nrow(human),human_accessions=data.table::uniqueN(human$uniprot_id),
    human_pathway_ids=data.table::uniqueN(human$pathway_id),mapped_human_membership_rows=sum(human$mapping_status=="mapped"),
    unmapped_human_accessions=nrow(unmapped),output_gene_accession_memberships=nrow(result),output_pathways=data.table::uniqueN(result$pathway_id),output_genes=data.table::uniqueN(result$gene_symbol)),
  outputs=lapply(paths,function(p) list(path=normalizePath(p),sha256=sha256_file(p))),
  source_url=a$source_url,source_version=a$source_version,license=a$license,
  regression_tests=list(human_filter=TRUE,duplicate_handling=TRUE,distinct_gene_counts=TRUE),
  interpretation="pathway annotation resource only; not enrichment, target direction, clinical efficacy or readiness validation")
jsonlite::write_json(summary,file.path(a$outdir,"build_summary.json"),pretty=TRUE,auto_unbox=TRUE,null="null")
cat("\nREACTOME SCRATCH BUILD COMPLETE\n")
cat("Native rows:",nrow(raw),"\nHuman membership rows:",nrow(human),"\nHuman accessions:",data.table::uniqueN(human$uniprot_id),"\n")
cat("Mapped human membership rows:",sum(human$mapping_status=="mapped"),"\nUnmapped human accessions:",nrow(unmapped),"\n")
cat("Output gene/accession memberships:",nrow(result),"\nOutput genes:",data.table::uniqueN(result$gene_symbol),"\nOutput pathways:",data.table::uniqueN(result$pathway_id),"\n")
cat("\nPathway sizes (distinct mapped genes):\n");print(base::summary(sizes$actual_distinct_genes))
cat("\nUnmapped accession examples:\n");print(head(unmapped[order(-n_native_membership_rows)],15))
cat("\nOutputs:",a$outdir,"\nNo .ready or propagation-enable sentinel created.\n")
