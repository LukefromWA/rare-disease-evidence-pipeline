# A11 synthetic DAG execution

Tested source commit: 3c84e71ae15aa03b2b755b5bf86e476bdffee02a

## Observed results
- PASS: synthetic regression suite regenerated its fixtures.
- PASS: DAG dry-run planned calculation, report and attachment.
- PASS: all processing rules executed successfully.
- PASS: producer/report/attachment provenance binding.
- PASS: 35 union hypotheses retained; 34 current; 33 A11-scored.
- PASS: unavailable A11 coverage retained rather than zero-imputed.
- PASS: second dry-run required no work.
- No real case inputs or retained Stan fits executed.
- Full A0-F pipeline NOT tested.

DAG test: 32 draws per specification; 18 specifications.
Separate regression suite: 64 draws per specification.
Production default: 10000 draws per specification.

## Environment
Snakemake: 9.27.1
Python: 3.14.8
R: R version 4.3.3 (2024-02-29)
R package versions: {"data.table": "1.14.10", "digest": "0.6.34", "jsonlite": "1.8.8"}
RNG kind: ["Mersenne-Twister", "Inversion", "Rejection"]

Resolved Python packages:
environment/a11-snakemake-test.requirements.txt

The package list records this tested Python environment. It is not a
complete cross-platform environment lock; restoration remains to be tested.

## Scope and privacy
Raw logs and runtime manifests are retained locally, not copied here.
No hostnames, absolute run paths, case evidence or report tables are included.
Synthetic execution is not biological or clinical validation.
Reporting findings remain open in docs/review/A11_FIRST_PASS.tsv.
