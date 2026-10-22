# workflow/lib/bootstrap.R -- sourced first by every node script: loads all lib files in dependency order
.lib_dir <- local({ f <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1]))) else getwd() })
.lib_dir <- if (nzchar(Sys.getenv("MVA_LIB_DIR"))) Sys.getenv("MVA_LIB_DIR") else if (basename(.lib_dir) == "lib") .lib_dir else normalizePath(file.path(.lib_dir, "..", "lib"))   # MVA_LIB_DIR overrides (tests)
for (.f in c("io.R", "manifest.R", "ids.R", "qc.R", "csq.R", "hpo.R", "network.R", "integration.R", "stability.R", "external_features.R", "gates.R",
             "track2_resources.R", "classes.R", "report_values.R")) source(file.path(.lib_dir, .f))
