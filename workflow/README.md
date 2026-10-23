# Analysis workflow

[Project overview](../README.md) · [Model specifications](../models/README.md) · [Shared functions](lib/)

The workflow separates case evidence, external model development, frozen inference, and Track 2 research. Snakemake connects the stages; the calculations live in the scripts and shared libraries.

## Start here

- [Master workflow](Snakefile): runs the A, C, E, and D stages using existing external model artifacts. It does not run external training.
- [Stage rules](rules/): input/output dependencies and commands for individual stages.
- [Analysis scripts](scripts/): executable preparation, ranking, inference, resource, and Track 2 scripts.
- [Shared functions](lib/): annotation parsing, QC, phenotype and network methods, integration, gates, and provenance helpers.
- [Standalone A11 workflow](layer_a/Snakefile): ranking-sensitivity analysis with its own [configuration and reproduction guide](../docs/A11_REPRODUCIBILITY.md).
- [Master workflow checkpoint](../docs/MASTER_WORKFLOW_CHECKPOINT.md): completed checks and remaining setup work.

## Stages and handoffs

| Stage | Role | Handoff |
|---|---|---|
| A | Prepare case evidence, apply QC, compare phenotype/network views, and rank candidates | Candidate tables with genomic evidence, coverage, and ranking summaries |
| A11 | Test ranking sensitivity under specified perturbations | Sensitivity summaries, separate from causal probabilities |
| B | Prepare external public cases, audit groups and splits, fit features and models, and freeze artifacts | External model and preprocessing artifacts |
| C | Apply the frozen feature/model contract without refitting on the case | Model-based candidate results |
| A ∪ C | Reconcile candidate identities while retaining both evidence views | Candidate registry for Track 2 |
| D | Examine mechanism scenarios, targets, compound evidence, gates, and sensitivity | Research hypotheses and their supporting evidence |
| E | Prepare and validate case-independent reference resources | Versioned resource bundle and validation records |
| F | Assemble saved outputs and provenance into reports | Reviewable tables, figures, and reports |

The master `all` target requests A10 ranked candidates, C4 hypotheses, and the D7 report table. PDF reporting and optional Track 2 targets are separate. A11 also has a separate entrypoint; its presence in the repository does not mean the master runs it automatically.

## External training and frozen inference

The B rule module and training sources are available for inspection, but the no-training master does not import that module. It reuses existing external artifacts for C. The [external-training branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training) contains the retained-fit documentation and diagnostic material.

The training sequence is B8 → B10 → B9: construct candidate sets, assign group-aware splits, then fit features on the training partition. C1 must use the same feature definitions and frozen preprocessing state. Held-out evaluation also requires the recorded deduplication, exclusion, and split-overlap audits.

Track 2 uses the retained A ∪ C registry, preserving source membership, genomic evidence, and the distinct score meanings. Missing C predictions are not evidence against an A-only candidate.

## Execution status

The repository master passed a local dry run with a 39-job plan through A10, C4, and D7. That check used an existing private configuration, resource snapshots, frozen artifacts, and a local R launcher. It checked workflow planning, not end-to-end execution or a fresh-checkout setup.

Portable setup still needs a complete master configuration example, environment restoration instructions, and resource preparation guidance. The imported R lock records one environment; it has not yet been verified as a complete, restorable environment for every stage. See the [checkpoint](../docs/MASTER_WORKFLOW_CHECKPOINT.md) for the detailed scope.

## Method status and data boundaries

Methods report unavailable, disabled, or not-applicable states when their inputs are missing; these states should not be filled with toy production evidence. In the current implementation, D1 carries mechanism priors forward without an evidence update, and D4 combines direction, pathway, and propagation evidence where available. Other research outputs should not be described as eight active aggregate-scoring channels.

Keep populated private configurations, case inputs, and generated case results outside the source release. Source review, runtime checks, and export review remain separate tasks. The presence of a rule or a successful dry run does not establish that its scientific results have been validated.
