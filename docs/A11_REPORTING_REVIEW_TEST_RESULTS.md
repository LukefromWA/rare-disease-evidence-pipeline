# A11 reporting/CLI review test results

Base source checkpoint: 695af5a172afda1755df69ab91acf67d277f786a.
PASS: new R syntax parsing.
PASS: all 16 coverage-presence combinations plus low-score three-view case.
PASS: rank-range boundaries; rank one does not imply robustness.
PASS: clean rejection of empty, odd-length and unknown CLI options.
PASS: original calculation regression and fixed-seed repeatability.
PASS: legacy report columns, candidate membership and all non-label values.
PASS: new report label expectations, input-binding and overwrite rejection.
PASS: source separation and retained union coverage.
PASS: standalone synthetic DAG and no-op second dry-run.
PASS: DAG label/provenance assertions; 35 union, 33 A11-scored.
Producer source unchanged; priors, seed and calculation not edited.
Regression used 64 draws/specification; DAG used 32; production default 10000.
No real case inputs, real predictor scoring or full A0-F execution.

Snakemake: 9.27.1
Local execution logs were not copied into the repository.
A11-R04 environment restoration remains open.
