# Local reviewed A/B import

Included:
- Private-run QC and integration helpers, copied without scientific edits.
- Four unchanged Stan definitions from external-training commit
  a31490ed5c3793bf73ab93bfeff46b3353f07ce0.
- An allowlisted numeric Layer A parameter record.
- Four synthetic core test groups.

This is a partial import, not the complete A0-A11 or B training workflow.
The numeric YAML is not a complete executable configuration.

Excluded:
- Original populated configs and private overlays.
- Patient inputs, identifiers, phenotype lists and case-derived results.
- Runtime manifests and commands containing private paths.
- The unparameterized A11 script and private runner.

The imported integration module expects a logger from its execution layer;
the unit tests supply a no-op logger. Legacy posterior column names do not
establish calibrated causal probabilities. Allele-fraction inference and
mosaic-review labels do not establish biological mosaicism.

A11 remains a separate Dirichlet/Beta sensitivity component. Its recorded
seed is 20260912; the inspected Layer A configuration seed is 20261014.
Do not silently replace either seed or conflate A11 with the A9 helper.

No complete workflow, Stan compilation, model fitting or privacy clearance
is established by these core tests. Record actual test output separately.
