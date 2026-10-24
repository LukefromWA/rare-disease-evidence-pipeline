# Shared library

[Project overview](../README.md) · [Workflow](../workflow/README.md) · [Model specifications](../models/README.md)

The shared workflow library is maintained on the `workflow` branch under
`workflow/lib`. This directory on `main` provides navigation and describes
the library's role in the analysis architecture.

## Implementation

| Looking for | Start with |
|---|---|
| Shared calculations and validation | [Workflow library](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib) |
| Analysis scripts and library call sites | [Workflow scripts](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/scripts) |
| Stage dependencies and commands | [Rule modules](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/rules) |
| Analysis stages and execution instructions | [Workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md) |
| Runtime configuration and supported execution routes | [CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md) |
| External-model specifications and retained-fit records | [External-training branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training) |

## Role in the workflow

The shared library supports reusable calculations and validation across analysis
stages. Node scripts coordinate configuration, input/output handling, and artifact
writing; shared functions keep common logic consistent across their call sites.

The workflow's shared-function responsibilities include:

- Configuration, schema and key validation, and explicit missing-input handling.
- QC, annotation parsing, and exact evidence joins.
- Phenotype similarity, network propagation, and coverage-aware integration.
- Dependence diagnostics, attribution, and ranking sensitivity.
- External-cohort preparation, split and feature contracts, and frozen
  prediction interfaces.
- Candidate reconciliation, mechanism scenarios, evidence integration, and gates.
- Provenance and manifest handling.

The implementation and its call sites provide the reference for individual
function interfaces, supported methods, and dependency requirements.

## Scientific contracts

Shared calculations preserve the distinctions between evidence types:

- Missing or unavailable evidence remains distinguishable from an observed
  negative result.
- Normalized evidence scores, ranking-sensitivity frequencies, and model
  probabilities retain their separate meanings.
- B9 training features and C1 inference features follow the matching feature
  definitions and retained preprocessing state.
- Frozen inference applies retained external-model artifacts without
  case-specific refitting.
- Candidate reconciliation preserves identifiers, source membership, coverage,
  and the distinct A and C evidence views.

## Development conventions

Reusable functions should expose explicit inputs and outputs, document their
assumptions and missing-data behavior, and support testing independently of a
complete workflow run.


