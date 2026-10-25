# Shared library — consolidation placeholder

> Directory distinction: this top-level `lib/` directory is a consolidation placeholder. Existing shared workflow source files are maintained in [`workflow/lib/` on the `workflow` branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib). The checklist below concerns this placeholder directory, not the implementation or validation status of the existing workflow library.


[Project overview](../README.md) · [Workflow](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md) · [Model specifications](../models/README.md)

This directory is reserved for the reviewed reusable functions that support the workflow. This README is a documentation placeholder; no library functions or package are created by this commit.

## Intended responsibilities

- Configuration, schema/key validation and explicit missing-input checks.
- Shared QC, annotation parsing and exact evidence joins.
- Phenotype similarity, network propagation and coverage-aware integration.
- Dependence diagnostics, attribution and ranking sensitivity.
- External-cohort preparation, split/feature contracts and frozen prediction interfaces.
- Candidate reconciliation, mechanism scenarios, evidence-family integration and hard gates.
- Provenance, manifests and synthetic test helpers.

Scientific functions should have explicit inputs and outputs and be testable without invoking the whole workflow. Node scripts should coordinate configuration, I/O, validation and artifact writing rather than duplicate the core calculations.

## Import checklist

- [ ] Import the existing reviewed implementations; do not replace them with placeholder calculations.
- [ ] Document each function's inputs, output schema, assumptions and missing-data behavior.
- [ ] Keep private inputs, secrets, absolute workspace paths and generated case results out of library code.
- [ ] Use the same feature construction and frozen state for B9 and C1.
- [ ] Keep normalized evidence scores, sensitivity frequencies and model probabilities distinct.
- [ ] Add synthetic unit/regression tests and dependency records.
- [ ] Check names and call sites against the imported workflow before documenting a runnable command.

No dummy scientific functions, synthetic performance claims or completed-validation claims are supplied here. Source consolidation and verification remain separate work.
