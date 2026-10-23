# A11 rewire test results

PASS: R syntax parsing of new R sources.

PASS: legacy calculation regression; fixed-seed repeatability; eligibility ledger; legacy report regression; input-mismatch rejection; overwrite rejection; union coverage; source separation

PASS: standalone synthetic DAG executed with Snakemake 9.27.1; second dry-run required no work. See A11_SYNTHETIC_DAG_RESULTS.md for scope and environment.

Synthetic regression used 64 draws per specification; production default remains 10000.
No real case data or retained Stan fits were executed.

Environment:

```text
R version 4.3.3 (2024-02-29) 
data.table 1.14.10 
digest 0.6.34 
jsonlite 1.8.8
```
