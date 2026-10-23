# Software and citations

This record separates the software observed in the local execution environments from the historical tools that produced annotations and external model fits. Installed versions alone do not establish that a package was used by every stage.

## Recorded environments

The following versions were collected locally on October 22, 2026. These observations describe the checkpoint environment, not a verified fresh-checkout installation.

### R runtime

| Software | Version | Role or scope |
|---|---|---|
| R | 4.3.3 (2024-02-29) | Analysis runtime |
| data.table | 1.14.10 | Table processing and joins |
| yaml | 2.3.8 | Configuration parsing |
| jsonlite | 1.8.8 | JSON input/output and manifests |
| digest | 0.6.34 | Checksums |
| igraph | 1.6.0 | Graph construction and network analysis |
| Matrix | 1.6-5 | Sparse matrix operations |
| renv | 1.0.3 | Environment management |

CmdStanR and posterior were not visible to the lean R launcher's library search paths during collection. Their historical training versions, along with the CmdStan version, still need to be verified from the external-training environment or saved provenance. Do not infer those versions from the no-training lock.

### Python environment

| Software | Version | Observation |
|---|---|---|
| Python | 3.14.8 | Interpreter in the Snakemake environment |
| Snakemake | 9.27.1 | Workflow orchestration |
| PyYAML | 6.0.3 | Installed; YAML is imported by the master workflow |
| requests | 2.34.2 | Installed; direct application use remains to be audited |

The recorded requirements include transitive dependencies. They should not all be described as independently selected scientific methods.

### Command-line tools

These versions were observed on the shell's PATH.

| Tool | Version |
|---|---|
| GNU Bash | 5.2.21 |
| Git | 2.43.0 |
| gzip | 1.12 |
| curl | 8.5.0 |
| BCFtools | 1.19 |
| tabix / HTSlib | 1.19 |
| bgzip / HTSlib | 1.19 |
| sha256sum / GNU coreutils | 9.4 |

The existing-VEP workflow consumes an annotated VCF, its index, and annotation metadata. The VEP software/cache version that produced those annotations must be taken from that metadata. Querying a tool's version for a manifest is not evidence that its analysis command ran.

## R and package references

The references below come from citation() and citation(package) in the observed R environment. Package versions are recorded above, including where a recommended paper predates the installed release.

- R Core Team (2024). R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/
- Barrett T, Dowle M, Srinivasan A (2023). data.table: Extension of data.frame. R package version 1.14.10. https://CRAN.R-project.org/package=data.table
- Garbett SP and contributors (2023). yaml: Methods to Convert R Data to YAML and Back. R package version 2.3.8. https://CRAN.R-project.org/package=yaml . The full author list is preserved in the locally generated BibTeX collection.
- Ooms J (2014). The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects. arXiv:1403.2805. https://arxiv.org/abs/1403.2805
- Eddelbuettel D and contributors (2024). digest: Create Compact Hash Digests of R Objects. R package version 0.6.34. https://CRAN.R-project.org/package=digest . The full contributor list is preserved in the locally generated BibTeX collection.
- Csardi G, Nepusz T (2006). The igraph software package for complex network research. InterJournal, Complex Systems, 1695. https://igraph.org
- Bates D, Maechler M, Jagan M (2024). Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.6-5. https://CRAN.R-project.org/package=Matrix
- Ushey K, Wickham H (2023). renv: Project Environments. R package version 1.0.3. https://CRAN.R-project.org/package=renv

The installed igraph citation output also supplied a manual entry dated 2026 for version 1.6.0. Its version-specific bibliographic metadata remains under review; the observed package version is not in dispute. Raw BibTeX entries also need unique citation keys before publication.

## Additional software references to verify

- Snakemake: follow the project's [citation guidance](https://snakemake.readthedocs.io/en/stable/project_info/citations.html).
- Stan, CmdStan and CmdStanR: identify the historical training versions and use the corresponding [CmdStanR citation information](https://mc-stan.org/cmdstanr/authors.html) and Stan project references.
- BCFtools: Danecek P et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience 10(2), giab008. https://doi.org/10.1093/gigascience/giab008 . This shared paper does not imply that SAMtools was used in this workflow.
- HTSlib: Bonfield JK et al. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience 10(2), giab007. https://doi.org/10.1093/gigascience/giab007
- VEP: verify the annotation provenance and add the software reference appropriate to the tool that produced the retained annotations.
- Python and directly imported third-party packages: record their role and version, then verify available project-recommended references.

Reference datasets, annotation snapshots, ontologies, pathways and interaction resources need their own citations, versions, licenses and checksums. They are not covered by citing the R/Python package that reads them.

## Source and validation status

The current implementation is on the [workflow branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/workflow). Its [master workflow checkpoint](https://github.com/LukefromWA/rare-disease-evidence-pipeline/blob/workflow/docs/MASTER_WORKFLOW_CHECKPOINT.md) records the successful local 39-job dry-run and its limitations.

This is an initial software inventory, not a completed repository-wide dependency audit. Remaining work includes checking all direct R/Python imports and subprocess calls, separating optional and training dependencies, reviewing the bibliography, and testing portable environment restoration. No private paths, sample identifiers, populated configurations or case results are included here.
