# workflow/lib/ids.R -- candidate identifiers and missing-value helpers
suppressPackageStartupMessages(library(data.table))
# One id form per table ("chrom:pos:ref:alt", no chr prefix); the "GRCh38:" prefix is added only when joining.
make_candidate_id <- function(chrom, pos, ref, alt) paste0(sub("^chr", "", chrom), ":", pos, ":", ref, ":", alt)
normalize_candidate_id <- function(x) {
  x <- trimws(as.character(x)); x[x == ""] <- NA_character_
  has_build <- grepl("^GRCh38:", x, ignore.case = TRUE)
  x[!is.na(x) & !has_build] <- paste0("GRCh38:", x[!is.na(x) & !has_build]); x
}
missing_value <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA", ".", "-")
