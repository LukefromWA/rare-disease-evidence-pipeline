# Shared workflow library

[Project overview](../README.md) · [Workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md) · [Shared library source](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib)

## Source location

The shared workflow source is maintained in
[`workflow/lib/` on the `workflow` branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow/workflow/lib).

This top-level `lib/` directory on `main` is a documentation and navigation
location. It is not a separate library package or the runtime source
directory. Its contents do not describe the implementation or validation
status of the existing workflow library.

## Organization

Shared calculations and validation belong in the workflow library.
Executable stages in `workflow/scripts/` coordinate configuration,
input/output handling and artifact writing, using shared functions
where implemented.

The workflow separates case-evidence preparation, external model
development, frozen inference and downstream research analysis.
Function-level behavior and availability should be read from the
corresponding source and documentation, not inferred from this directory.

## Running the workflow

Use the
[workflow guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/workflow/README.md)
and
[CLI guide](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/CLI.md)
for entrypoints, dependencies, configuration and required inputs.

No source consolidation or relocation is required by this README.
