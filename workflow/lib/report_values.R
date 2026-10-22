# workflow/lib/report_values.R -- the only path from a result to LaTeX: generated macros and tables (never hand-typed numbers)
suppressPackageStartupMessages(library(data.table))
tex_escape <- function(x) gsub("([#$%&_{}])", "\\\\\\1", as.character(x))
write_values_tex <- function(values, path) {
  lines <- sprintf("\\expandafter\\def\\csname mvaval@%s\\endcsname{%s}", gsub("_", "-", names(values)), vapply(values, tex_escape, ""))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); writeLines(c("% GENERATED. Do not edit.", lines), path); invisible(path)
}
write_latex_table <- function(dt, path, caption = NULL, digits = 3, max_rows = 40) {
  d <- data.table::as.data.table(utils::head(dt, max_rows)); for (j in names(d)) if (is.numeric(d[[j]])) data.table::set(d, j = j, value = signif(d[[j]], digits))
  hdr <- paste(tex_escape(names(d)), collapse = " & "); rows <- apply(d, 1, function(r) paste(tex_escape(ifelse(is.na(r), "", r)), collapse = " & "))
  body <- c(sprintf("\\begin{longtable}{%s}", paste(rep("l", ncol(d)), collapse = "")), "\\toprule", paste0(hdr, " \\\\"), "\\midrule", paste0(rows, " \\\\"), "\\bottomrule",
            if (!is.null(caption)) sprintf("\\caption{%s}", tex_escape(caption)), "\\end{longtable}")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE); writeLines(c("% GENERATED. Do not edit.", body), path); invisible(path)
}
