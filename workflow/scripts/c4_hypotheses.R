#!/usr/bin/env Rscript
.sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]); source(file.path(dirname(normalizePath(.sf)), "..", "lib", "bootstrap.R"))
# C4: full weighted hypothesis set. q_j over EVERY candidate; p_m from the declared mechanism prior; never collapsed to rank 1.
a <- parse_args(c("config", "sample", "a10", "m0m1", "frozen_scores", "outdir"), c("evidence_matrix")); cfg <- read_config(a$config)
a10 <- read_gz(a$a10); m <- read_gz(a$m0m1); fs <- read_gz(a$frozen_scores)
d <- merge(merge(a10[, .(candidate_id = best_candidate_id, gene = gene_symbol, impact = best_IMPACT, consequence = if ("best_consequence" %in% names(a10)) best_consequence else NA_character_)], m[, .(candidate_id, m0, m1, rank_m0, rank_m1)], by = "candidate_id"),
           fs[, .(candidate_id, p_frozen, rank_frozen, in_conformal_set)], by = "candidate_id")
q <- 0.5 * (0.5 * d$m0 + 0.5 * d$m1) + 0.5 * d$p_frozen; d[, q_j := q / sum(q)]
d[, cclass := consequence_class(consequence, impact)]
mp <- mechanism_prior_table(cfg$mechanism_prior, unique(d$cclass))
fh <- merge(d, mp, by = "cclass", allow.cartesian = TRUE); setnames(fh, "gene", "gene"); fh[, cclass := NULL]
obj <- new("HypothesisSet", data = fh[, .(gene, variant = candidate_id, q_j, p_m, mechanism)])        # S4 validity: sum(q_j)=1 over hypotheses, p_m in [0,1]
fh[, variant := candidate_id]; setorder(fh, -q_j, -p_m)
od <- a$outdir; dir.create(od, recursive = TRUE, showWarnings = FALSE)
write_table(d[order(-q_j)], file.path(od, "hypothesis_weights.tsv.gz")); write_table(fh, file.path(od, "final_hypothesis_set.tsv.gz"))
write_node_manifest(od, "c4", cfg, inputs = c(a$a10, a$m0m1, a$frozen_scores), outputs = file.path(od, c("hypothesis_weights.tsv.gz", "final_hypothesis_set.tsv.gz")),
  extra = list(n_hypotheses = nrow(d), sum_q_j = sum(d$q_j), mechanism_prior_mode = cfg$mechanism_prior$mode, mechanism_prior_is_biological = identical(cfg$mechanism_prior$mode, "declared")))
logmsg("C4: ", nrow(d), " hypotheses; sum q_j = ", round(sum(d$q_j), 6), "; top: ", paste(head(d[order(-q_j), gene], 3), collapse = ", "))
