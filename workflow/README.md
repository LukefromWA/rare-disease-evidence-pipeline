# Analysis workflow

[Project overview](../README.md) · [CLI guide](../docs/CLI.md)

## Entrypoints

- `../bin/mva`: interactive Track 1 / Track 2 / connected-run selection.
- `Snakefile`: canonical configured no-training workflow.
- `track2_genes.Snakefile`: independent unweighted gene-only Track 2 discovery.
- `layer_a/Snakefile`: standalone A11 sensitivity/reporting workflow only.

Rules live in `rules/`; calculations live in `scripts/` and `lib/`.
Layer directories are not six interchangeable independent runners.

## Master targets

- `track1`: A10, frozen Stan predictions and A11 in fresh-Stan mode.
- `layer_a`: through A10.
- `layer_c`: selected C4 handoff.
- `layer_d`: through D7.
- `resources`: native preparation and core resource validation.
- `all`: resource readiness, A10, selected C4 handoff and D7.

`layer_b` requests an existing artifact; it is not a training command.
The master does not import external-training rules.
E21 writes a separate validation-status report and does not imply enablement.

## Execution and reporting scope

A configured existing-VEP resource/case workflow completed locally after repair
and resume. Output integrity and cross-stage consistency were checked. This does
not establish fresh-checkout reproducibility, clinical validity or a clean
uninterrupted runtime benchmark.

The CLI gene-only route was executed, and all five weighted D2 discovery tables
remained content-identical in a regression run.
The expanded interactive/raw-annotation paths require their own integration tests.

D7's safety-gated research table is separate from optional PDF reporting.
Legacy PDF asset wiring still requires alignment with the fresh-Stan handoff;
do not describe PDF generation as part of the demonstrated all target.

See the CLI guide for input requirements, historical-support limitations, private
configuration handling, path support and the current validation boundary.
