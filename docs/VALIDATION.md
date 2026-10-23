# Validation scope

The reusable implementation is published without the internal synthetic
test scripts. Those scripts and execution logs are retained locally.

Completed synthetic checks:
- Original A11 calculation and non-label report-value regression.
- Fixed-seed repeatability and candidate coverage preservation.
- Input/provenance mismatch and overwrite rejection.
- Separate Stan and A11 source weights.
- Coverage combinations, neutral rank labels and CLI edge cases.
- Standalone A11 calculation -> report -> union attachment DAG.
- No-work second dry-run.

Tests used reduced simulation counts. Production defaults are unchanged.
No real case inputs, real retained-predictor scoring, full A0-F execution
or biological/clinical validation is established by these checks.

Detailed existing validation records describe execution history.
The exact internal test suite is not distributed in this source tree;
independent reproduction of that suite requires obtaining it separately.
Environment restoration and the full master workflow remain pending.
