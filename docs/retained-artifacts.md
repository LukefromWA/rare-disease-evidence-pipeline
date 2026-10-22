# Retained fit artifacts

The `retained-folds-v1` release asset set contains twelve original retained objects: three folds each for M2a, M2b, M3a and M3b. It also contains `fold_assets.csv` and `fold_assets.sha256`.

The asset set is uploaded and checksum-verified.
[Open the retained-fold release](https://github.com/LukefromWA/rare-disease-evidence-pipeline/releases/tag/retained-folds-v1).
Downloads become available when the draft is published; repository
access permissions apply.

## Contents and verification

Each model/fold page records its artifact name, byte count and SHA-256. Verify downloads using `sha256sum -c fold_assets.sha256` from the download directory. The compact asset inventory is also available as [CSV](retained-fold-assets.csv).

The three-component M3b predictor and across-specification averaging are different constructions. The retained downstream adapter uses the former. These fitted objects retain their recorded diagnostic qualifications; see the model pages and diagnostic tables.
