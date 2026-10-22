# Workflow — consolidation placeholder

[Project overview](../README.md) · [Model specifications](../models/README.md) · [Shared library](../lib/README.md)

This directory is reserved for the reviewed workflow implementation. This README is a documentation placeholder, not a runnable entrypoint or evidence of a completed reproduction.

## Layer responsibilities

| Layer | Role and handoff |
|---|---|
| A | Case evidence preparation, QC, phenotype/network views, integration and candidate handoff |
| A11 | Ranking sensitivity under declared perturbations; distinct from causal probability |
| B | External public cohort preparation, deduplication/split audits, feature fitting, model development and freezing |
| C | Apply the frozen feature/model contract without case-specific refitting |
| A ∪ C | Reconcile candidate identities and retain both source views for Track 2 |
| D | Mechanism scenarios, complementary discovery methods, evidence integration, hard gates and sensitivity |
| E | Case-independent resource preparation, versioning and validation |
| F | Run-specific tables, figures, provenance and reviewed export |

## Import checklist

- [ ] Add reviewed node scripts, orchestration rules and declared input/output contracts.
- [ ] Add matching shared functions in `lib` without duplicating scientific logic across nodes.
- [ ] Preserve the B8 → B10 → B9 ordering: candidate sets, group-aware splits, train-fitted features.
- [ ] Verify B9/C1 feature parity and frozen inference behavior.
- [ ] Preserve the A ∪ C candidate union and distinct score meanings.
- [ ] Link the actual deduplication rules, exclusions and split-overlap audit records.
- [ ] Add synthetic tests, dependency records and a tested reproduction command.
- [ ] Record unavailable/disabled methods explicitly instead of supplying toy production evidence.
- [ ] Complete the export review before adding data-derived artifacts.

Only this README is added here by the documentation commit. It does not add a Snakefile, model source, private configuration, patient data or executable workflow. The existing Stan definitions remain available on the [external-training branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/workflow/stan).
