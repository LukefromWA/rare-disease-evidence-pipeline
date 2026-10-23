# Master workflow source checkpoint

## Scope

Imports the local no-training master workflow, Track 2 and resource rules,
supporting source scripts, and the recorded no-training R dependency lock.

The existing integration helper is preserved, including its
best_genomic_score export. External training is not included in the master
workflow's rule imports; existing external artifacts are reused.

## Checks actually performed

- Imported 63 new source/environment files from a checksum-recorded local bundle.
- Checked changed/new source for selected private identifiers, machine-specific
  paths, and credential patterns. Two flags were URL parameter labels,
  not embedded credentials. This is not comprehensive privacy certification.
- Parsed new Python source successfully with ast.parse.
- git diff --cached --check returned exit code 0 before this document was added.
- Repository workflow/Snakefile dry-run returned exit code 0.
- The dry-run planned 39 jobs, including the all target, reaching A10, C4, and D7.
- No analysis jobs were executed by that dry-run.

## Reproduction boundary

The successful dry-run used an existing private configuration, sample sheet,
local resource snapshots, frozen external artifacts, and a local R launcher.
It does not establish fresh-checkout or end-to-end scientific reproduction.

Still required for portable execution:
- A complete unpopulated master configuration template and input contracts.
- Documented restoration and verification of the R/Python/system environment.
- Configurable resource snapshot paths and acquisition/preparation instructions.
- Documentation of optional Track 2 targets outside the default all target.
- Runtime tests and a clean-checkout verification with appropriately
  distributable inputs.

The R lock records one environment; its presence alone does not establish
that all workflow dependencies are covered or successfully restorable.

Private configurations, case data, resource bundles, and generated case
outputs are not part of this source import. F report/PDF targets are separate
from the master's default all target. Repository status documentation may
describe earlier checkpoints; this document records this import's scope.
