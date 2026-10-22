# Local synthetic core test results

Executed UTC: 2026-10-22T23:20:14.543029+00:00

Command: `Rscript --vanilla tests/layer_a/test_core.R`

Exit status: 0

```text
R version 4.3.3 (2024-02-29)
data.table 1.14.10
PASS: QC boundaries, informative VAF and copy semantics
PASS: name-aligned weights and missing evidence
PASS: empty-panel informed/agnostic agreement
PASS: Shapley efficiency with missing evidence
```

Scope: synthetic QC/integration core only. No patient inputs, resource downloads or model fitting were used by this harness.

This is not full A0-A11 validation, B training/deployment validation, Stan compilation, or comprehensive privacy clearance.
