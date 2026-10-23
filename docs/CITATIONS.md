# Literature, data-source and software citations

Scope: public implementation and methods attribution for the rare-disease evidence pipeline. This register contains no patient identifiers, alleles, genomic coordinates, private paths, case results or literature-to-patient links. A citation establishes attribution, not scientific validation, case applicability, permission to redistribute a dataset, or proof that every described feature was used.

## How to use this register

- Cite the method paper and the actual data release separately where appropriate.
- Record resource releases, access dates, checksums and preparation steps in the applicable run manifest. A paper's year is not a resource snapshot version.
- Use the citation returned by `citation("package")` in the environment that actually ran the analysis. Current online package pages may describe a different version.
- Observed software versions are recorded in [software_versions.observed.tsv](../environment/software_versions.observed.tsv). The [R lock](../environment/no-training.renv.lock) and [Python requirements](../environment/a11-snakemake-test.requirements.txt) remain the fuller dependency records, including transitive packages. They are not proof of restored or historical training environments.
- Detailed case-specific allele justification stays in the authorized submission rather than being copied into this public register.

## Readthrough model and exported data

1. Toledano I, Supek F, Lehner B. Genome-scale quantification and prediction of pathogenic stop codon readthrough by small molecules. *Nature Genetics*. 2024;56:1914-1924. DOI: [10.1038/s41588-024-01878-5](https://doi.org/10.1038/s41588-024-01878-5).
2. Toledano I, Supek F, Lehner B. *Readthrough predictions*. Figshare. DOI: [10.6084/m9.figshare.23708901](https://doi.org/10.6084/m9.figshare.23708901). This is the verified base dataset identifier. Exact version-specific citation must be checked against the pinned dataset metadata; no version-specific DOI is asserted here.

RTDetective concerns stop-codon sequence-context readthrough, not a direct nonsense-mediated-decay model or a clinical response probability. Dataset lookups, published model performance and functional experiments are separate evidence. This register does not verify an exported score scale or infer treatment rankings.

## Experimental comparator literature

3. Trzaska C, Amand S, Bailly C, et al. 2,6-Diaminopurine as a highly potent corrector of UGA nonsense mutations. *Nature Communications*. 2020;11:1509. DOI: [10.1038/s41467-020-15140-z](https://doi.org/10.1038/s41467-020-15140-z). PMID: 32198346.
4. Friesen WJ, Trotta CR, Tomizawa Y, et al. The nucleoside analog clitocine is a potent and efficacious readthrough agent. *RNA*. 2017;23:567-577. DOI: [10.1261/rna.060236.116](https://doi.org/10.1261/rna.060236.116). PMID: 28096517.
5. Sharma J, Du M, Wong E, et al. A small molecule that induces translational readthrough of CFTR nonsense mutations by eRF1 depletion. *Nature Communications*. 2021;12:4358. DOI: [10.1038/s41467-021-24575-x](https://doi.org/10.1038/s41467-021-24575-x). PMID: 34272367.
6. Coelho JPL, Yip MCJ, Oltion K, Taunton J, Shao S. The eRF1 degrader SRI-41315 acts as a molecular glue at the ribosomal decoding center. *Nature Chemical Biology*. 2024;20:877-884. DOI: [10.1038/s41589-023-01521-0](https://doi.org/10.1038/s41589-023-01521-0). PMID: 38172604.

These sources concern their studied experimental systems, not established clinical suitability. DAP's reported low experimental toxicity must not be replaced by a blanket statement that all alternative readthrough compounds are more toxic. SRI-41315 structural findings must not be assigned to ataluren.

## Ataluren literature

7. Ng MY, Li H, Ghelfi MD, Goldman YE, Cooperman BS. Ataluren and aminoglycosides stimulate read-through of nonsense codons by orthogonal mechanisms. *PNAS*. 2021;118:e2020599118. DOI: [10.1073/pnas.2020599118](https://doi.org/10.1073/pnas.2020599118). PMID: 33414181.
8. Bezzerri V, Bardelli D, Morini J, et al. Ataluren-driven restoration of Shwachman-Bodian-Diamond syndrome protein function in Shwachman-Diamond syndrome bone marrow cells. *American Journal of Hematology*. 2018;93:527-536. DOI: [10.1002/ajh.25025](https://doi.org/10.1002/ajh.25025). PMID: 29285795.
9. Bezzerri V, Pegoraro A, Hristodor AM, et al. Ataluren improves hematopoietic and pancreatic disorders in Shwachman-Diamond syndrome patients: a compassionate program case-series. *Nature Communications*. 2025;16:8189. DOI: [10.1038/s41467-025-63137-3](https://doi.org/10.1038/s41467-025-63137-3). PMID: 40897730.
10. Vlodavets D, Wu S, Kostera-Pruszczyk A, et al. Confirmatory long-term efficacy and safety results of ataluren in patients with nmDMD from Study 041, an international, randomized, double-blind, placebo-controlled, Phase III trial. *Journal of Comparative Effectiveness Research*. 2025;14:e240238. DOI: [10.57264/cer-2024-0238](https://doi.org/10.57264/cer-2024-0238). PMID: 40843507.
11. Kerem E, Konstan MW, De Boeck K, et al. Ataluren for the treatment of nonsense-mutation cystic fibrosis: a randomised, double-blind, placebo-controlled phase 3 trial. *Lancet Respiratory Medicine*. 2014;2:539-547. DOI: [10.1016/S2213-2600(14)70100-6](https://doi.org/10.1016/S2213-2600(14)70100-6). PMID: 24836205.
12. Konstan MW, VanDevanter DR, Rowe SM, et al. Efficacy and safety of ataluren in patients with nonsense-mutation cystic fibrosis not receiving chronic inhaled aminoglycosides: The international, randomized, double-blind, placebo-controlled Ataluren Confirmatory Trial in Cystic Fibrosis (ACT CF). *Journal of Cystic Fibrosis*. 2020;19:595-601. DOI: [10.1016/j.jcf.2020.01.007](https://doi.org/10.1016/j.jcf.2020.01.007). PMID: 31983658.

Publications, outcome claims and overlapping subgroups are not independent study counts. Source-anchor matching is not completed scientific adjudication. Preserve unfavorable findings and endpoint hierarchy. Study 041's favorable ITT analyses do not replace its nonsignificant primary mITT result. Patient-cell findings and the small uncontrolled SDS series do not establish benefit in another disorder.

## Regulatory context (checked 23 October 2026)

- [UK Translarna SmPC](https://www.medicines.org.uk/emc/product/2576/smpc), updated 31 March 2026: conditional authorization for specified ambulatory nonsense-mutation Duchenne patients aged 2 years and older; not approval for another indication.
- [EMA Translarna EPAR](https://www.ema.europa.eu/en/medicines/human/EPAR/translarna): EU authorization not renewed on 28 March 2025 because effectiveness was not confirmed.
- [PTC regulatory update, 12 February 2026](https://ir.ptcbio.com/news-releases/news-release-details/ptc-therapeutics-provides-regulatory-update-translarnatm-2): US NDA resubmission withdrawn following effectiveness feedback; not withdrawal of a previously approved US product.

Authorization, availability, case-specific access and safety are different questions. No treatment-initiation or dosing recommendation is made.

## Public annotation, phenotype and network resources

| Resource | Citation or official citation guidance | Reproduction detail to retain |
|---|---|---|
| Ensembl VEP | McLaren et al. *The Ensembl Variant Effect Predictor*. Genome Biology 2016;17:122. DOI: [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4) | VEP/cache version, reference build, flags and plugin versions |
| HPO | Gargano et al. *The Human Phenotype Ontology in 2024: phenotypes around the world*. NAR 2024;52:D1333-D1346. DOI: [10.1093/nar/gkad1005](https://doi.org/10.1093/nar/gkad1005) | Ontology and gene-annotation snapshot |
| MANE | Morales et al. *A joint NCBI and EMBL-EBI transcript set for clinical genomics and research*. Nature 2022;604:310-315. DOI: [10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8) | MANE release and exact transcript versions |
| Phenopacket Store | Danis et al. *A corpus of GA4GH phenopackets: Case-level phenotyping for genomic diagnostics and discovery*. HGG Advances 2025;6:100371. DOI: [10.1016/j.xhgg.2024.100371](https://doi.org/10.1016/j.xhgg.2024.100371) | Actual release, inclusion policy, grouping and split audit |
| STRING | [Official citation list](https://string-db.org/cgi/about); 2023 resource paper DOI [10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000); 2025 paper DOI [10.1093/nar/gkae1113](https://doi.org/10.1093/nar/gkae1113) | Exact release, network file/type and confidence threshold; citing a paper does not enable all its features |
| Reactome | [Official citation instructions](https://reactome.org/cite) | Release, mappings and stable pathway IDs; pathway-specific citation where applicable |
| UniProt | [Official publications/citation guidance](https://www.uniprot.org/help/publications) | Release, reviewed-status policy and mapping file |
| SIGNOR | Lo Surdo et al. *SIGNOR 3.0, the SIGnaling network open resource 3.0: 2022 update*. NAR 2023;51:D631-D637. DOI: [10.1093/nar/gkac883](https://doi.org/10.1093/nar/gkac883) | Actual snapshot/schema, species, signs and evidence mapping; paper is not proof of snapshot version |
| ChEMBL | Zdrazil et al. *The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods*. NAR 2024;52:D1180-D1192. DOI: [10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004). Web services: Davies et al. 2015, DOI [10.1093/nar/gkv352](https://doi.org/10.1093/nar/gkv352) | Actual release, native acquisition and action normalization; applicable source license |
| gnomAD | [Resource](https://gnomad.broadinstitute.org/) and [maintainer citation guidance](https://github.com/broadinstitute/gnomad-browser/issues/1832) | Exact major/minor release, genome/exome subset and access date; do not substitute a current browser release for the historical lookup |
| ClinVar | Landrum et al. *ClinVar: public archive of relationships among sequence variation and human phenotype*. NAR 2014;42:D980-D985. DOI: [10.1093/nar/gkt1113](https://doi.org/10.1093/nar/gkt1113); [citation guidance](https://www.ncbi.nlm.nih.gov/clinvar/docs/faq/) | Snapshot/query date and accession/version where a record was used |
| REVEL | Ioannidis et al. *REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants*. AJHG 2016;99:877-885. DOI: [10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016) | Score release and actual retrieval source |
| SIFT | Ng and Henikoff. *SIFT: predicting amino acid changes that affect protein function*. NAR 2003;31:3812-3814. DOI: [10.1093/nar/gkg509](https://doi.org/10.1093/nar/gkg509) | Annotation/cache version; do not count correlated predictors as independent evidence |
| PolyPhen-2 | Adzhubei et al. *A method and server for predicting damaging missense mutations*. Nature Methods 2010;7:248-249. DOI: [10.1038/nmeth0410-248](https://doi.org/10.1038/nmeth0410-248) | Score model and annotation/cache source |

These are resource attributions, not a complete locked-data manifest. Dynamic source websites and generic resource papers do not certify exact historical releases or redistribution rights.

## Software papers and recommended citations

- R Core Team (2024). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna. [R project](https://www.R-project.org/). Observed R version: 4.3.3; not a measurement of the original training environment.
- Moelder F, Jablonski KP, Letcher B, et al. *Sustainable data analysis with Snakemake*. F1000Research 10:33. [Version 3, 2025](https://doi.org/10.12688/f1000research.29032.3); first published in 2021. The observed executable version is recorded separately.
- Danecek P, Bonfield JK, Liddle J, et al. *Twelve years of SAMtools and BCFtools*. GigaScience 2021;10:giab008. DOI: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008). Citing this joint paper does not assert that SAMtools was used by every stage.
- Bonfield JK, Marshall J, Danecek P, et al. *HTSlib: C library for reading/writing high-throughput sequencing data*. GigaScience 2021;10:giab007. DOI: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007). Relevant to HTSlib tooling including bgzip/tabix; log actual tool versions.
- Carpenter B, Gelman A, Hoffman MD, et al. *Stan: A Probabilistic Programming Language*. JSS 2017;76(1):1-32. DOI: [10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01).
- Csardi G, Nepusz T (2006). *The igraph software package for complex network research*. InterJournal, Complex Systems:1695. [igraph](https://igraph.org/).
- Ooms J (2014). *The jsonlite Package: A Practical and Consistent Mapping Between JSON Data and R Objects*. [arXiv:1403.2805](https://arxiv.org/abs/1403.2805).

### R package attribution

The following direct/support package citations were supplied by the participant's R citation export or identified from package documentation. Version-specific installed-package metadata takes precedence over today's website. Locked transitive packages are listed in `environment/no-training.renv.lock` and should retain their own package attribution when used; this is not a claim that every installed dependency was called directly.

| Package | Attribution / official citation metadata |
|---|---|
| data.table | Barrett T, Dowle M, Srinivasan A. *data.table: Extension of data.frame*. Participant export: 2023, version 1.14.10. [Package](https://cran.r-project.org/package=data.table) |
| yaml | Garbett SP, Stephens J, Simonov K, et al. *yaml: Methods to Convert R Data to YAML and Back*. Participant export: 2023, version 2.3.8. [Package](https://cran.r-project.org/package=yaml) |
| jsonlite | Ooms (2014), paper above; participant observed version 1.8.8. [Package](https://cran.r-project.org/package=jsonlite) |
| digest | Eddelbuettel D, et al. *digest: Create Compact Hash Digests of R Objects*. Participant export: 2024, version 0.6.34. [Package](https://cran.r-project.org/package=digest) |
| igraph | Csardi and Nepusz (2006), above. Participant observed version 1.6.0; generate release-specific citation in the actual environment. |
| Matrix | Bates D, Maechler M, Jagan M. *Matrix: Sparse and Dense Matrix Classes and Methods*. Participant export: 2024, version 1.6-5. The version probe rendered 1.6.5; lock/citation spelling is preserved separately. [Package](https://cran.r-project.org/package=Matrix) |
| renv | Ushey K, Wickham H. *renv: Project Environments*. Participant export: 2023, version 1.0.3. [Package](https://cran.r-project.org/package=renv) |
| R6 | [Package and citation metadata](https://cran.r-project.org/package=R6); lock records 2.5.1. |
| loo | [Package and citation metadata](https://cran.r-project.org/package=loo); lock records 2.6.0. Installation does not establish that stacking or LOO was exercised in a particular run. |
| matrixStats | [Package and citation metadata](https://cran.r-project.org/package=matrixStats); lock records 1.5.0. |
| cmdstanr | Gabry J, Cesnovar R, Johnson A, Bronder S. *cmdstanr: R Interface to CmdStan*. [Official citation](https://mc-stan.org/cmdstanr/authors.html). Training version not established by the no-training probe. |
| posterior | Buerkner PC, Gabry J, Kay M, Vehtari A. *posterior: Tools for Working with Posterior Distributions*. [Official citation metadata](https://cloud.r-project.org/web/packages/posterior/citation.html). Training version not established by the no-training probe; do not assign the current online version to historical fits. |

R base/recommended namespaces (for example stats, utils, methods, graphics, grDevices and parallel) are covered by R plus any package-specific citation returned by the installed environment. The lock also records backports, checkmate, cli, cpp11, lattice, lifecycle, magrittr, pkgconfig and rlang; use each installed package's `citation()` output if explicitly used. No new versions or authors are inferred for those entries here.

### Python and command-line attribution

- [Python](https://www.python.org/): cite the interpreter version used.
- [Snakemake](https://snakemake.readthedocs.io/): paper above; workflow planning is distinct from executed scientific analysis.
- [PyYAML](https://pypi.org/project/PyYAML/) and [Requests](https://pypi.org/project/requests/): record package versions and their actual role. An installed package is not proof of direct execution by project code.
- [Bash](https://www.gnu.org/software/bash/), [Git](https://git-scm.com/), [GNU gzip](https://www.gnu.org/software/gzip/), [curl](https://curl.se/) and [GNU coreutils](https://www.gnu.org/software/coreutils/): version records are in the observed table; this register does not invent scientific papers for these utilities.
- VEP: citation above; its version was not collected by the supplied CLI probe. CmdStan is separate from CmdStanR and needs its own observed version in the training environment.

## Capture and interpretation boundaries

The observed version table is transcribed from participant-supplied terminal output and the committed lock; no environment was probed by the assistant. The supplied no-training R probe could not resolve cmdstanr/posterior versions. That does not show that the historical training environment lacked them. A directory name or current package website is not a version measurement.

This register does not license third-party source/data for redistribution, alter claim-review status, or validate biological causality, functional rescue, clinical benefit or safety. Keep private configurations, genotype-scale outputs, case-specific manifests and local inventory paths outside the public source release. Do not upload model/resource artifacts merely to supply their citations.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

See [AI assistance disclosure](AI_ASSISTANCE.md). The formal challenge-dataset publication citation is separate from these resource/software citations; do not invent a DOI when the organizer-prescribed reference is unavailable.
