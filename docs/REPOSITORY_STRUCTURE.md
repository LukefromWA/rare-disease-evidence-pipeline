# Repository organization

main is the accepted baseline; workflow is the integration branch.
Layer development branches: layer-a through layer-f.

Target directories: workflow/layer_a through workflow/layer_f.
Shared dependencies remain under workflow/lib.

The A/B/C source checkpoint has been preserved on workflow.
Functional source relocation requires corresponding bootstrap, rule,
helper/model path, documentation and checksum updates.
Existing checkpoint documents describe their original import commits.

B/C completeness, A11 integration, remaining D/E/F imports and runnable
entrypoints remain explicitly tracked work, not assumed completed.
No scientific implementation or repository visibility was changed by
this scaffold commit.
