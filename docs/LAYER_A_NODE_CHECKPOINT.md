# Layer A node-source checkpoint — 2026-10-22

## Scope
A0–A10 node scripts, supporting manifest scripts and a_private.smk
copied unchanged from Written Final Scripts/workflow.
Shared libraries were imported in the preceding source checkpoint.

The local cmdstan-backend checkout helped identify the node filenames;
it was not substituted for the selected source tree.

## Important distinctions
This is the Written Final Scripts implementation.
It is not a claim that the separate private-run A implementation has
been fully imported or that its scientific behavior is identical.

Existing workflow/layer_a/R extracted helpers remain unchanged.
Their test results do not establish tests of these newly imported nodes.

## Not included
- Separate private-run runner and A11 implementation.
- Populated configs, sample tables, phenotype lists or genomic inputs.
- Results, logs, fitted objects or runtime manifests.
- Master Snakefile and Layer D rules.

## Execution status
Source copied only. No node execution, VEP run, fitting or new test pass.
a_private.smk depends on variables supplied by a future entrypoint.
It is not a standalone layer runner.
Configuration, environment and resource setup remain integration tasks.

## Next work
Connect the reviewed A11 path after removing private runtime settings.
Add layer/master entrypoints, simple case inputs and connected tests.
Preserve the distinction between the two A implementation families.
Original working files are not modified.
