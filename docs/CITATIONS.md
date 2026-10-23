# Literature, data-source and software citations

Scope: attribution for the public implementation and methods of the rare-disease evidence pipeline. This register contains no case identifiers, alleles, coordinates, private paths, case results or literature-to-patient links. A citation documents attribution; it does not establish scientific validity, case applicability or redistribution permission.

## Use and versioning

- Cite method papers and the exact data release separately where relevant.
- Run manifests should capture resource releases, access dates, checksums and transformations. A paper's publication year is not a data snapshot version.
- Use `citation("package")` in the environment that actually ran the analysis. Current package websites may show a different release.
- The [observed software-version table](../environment/software_versions.observed.tsv) scopes participant-observed versions. [R lock](../environment/no-training.renv.lock) and [Python requirements](../environment/a11-snakemake-test.requirements.txt) remain the fuller dependency records. Neither proves a historical training environment was restored or used.
- Keep case-specific allele justification in authorized submission materials, not this public register.

## Readthrough model and dataset

1. Toledano I, Supek F, Lehner B. Genome-scale quantification and prediction of pathogenic stop codon readthrough by small molecules. *Nature Genetics*. 2024;56:1914–1924. DOI: [10.1038/s41588-024-01878-5](https://doi.org/10.1038/s41588-024-01878-5).
2. Toledano I, Supek F, Lehner B. *Readthrough predictions*. Figshare. Base DOI: [10.6084/m9.figshare.23708901](https://doi.org/10.6084/m9.figshare.23708901). The workflow declares dataset version 2. This base DOI was verified; confirm the exact version-specific citation against pinned inventory metadata before a later publication. No version-specific DOI is asserted here.

RTDetective predicts stop-codon sequence-context readthrough; it is not a direct NMD or clinical-response model. Published model performance, experimental drug activity and an imported exact-coordinate lookup are distinct evidence. This register does not verify exported score units or treatment rankings.

## Experimental readthrough comparators

3. Trzaska C, Amand S, Bailly C, et al. 2,6-Diaminopurine as a highly potent corrector of UGA nonsense mutations. *Nature Communications*. 2020;11:1509. DOI: [10.1038/s41467-020-15140-z](https://doi.org/10.1038/s41467-020-15140-z). PMID: 32198346.
4. Friesen WJ, Trotta CR, Tomizawa Y, et al. The nucleoside analog clitocine is a potent and efficacious readthrough agent. *RNA*. 2017;23(4):567–577. DOI: [10.1261/rna.060236.116](https://doi.org/10.1261/rna.060236.116). PMID: 28096517.
5. Sharma J, Du M, Wong E, et al. A small molecule that induces translational readthrough of CFTR nonsense mutations by eRF1 depletion. *Nature Communications*. 2021;12:4358. DOI: [10.1038/s41467-021-24575-x](https://doi.org/10.1038/s41467-021-24575-x). PMID: 34272367.
6. Coelho JPL, Yip MCJ, Oltion K, Taunton J, Shao S. The eRF1 degrader SRI-41315 acts as a molecular glue at the ribosomal decoding center. *Nature Chemical Biology*. 2024;20(7):877–884. DOI: [10.1038/s41589-023-01521-0](https://doi.org/10.1038/s41589-023-01521-0). PMID: 38172604.

These papers report findings in their studied systems, not established clinical suitability. DAP's reported low experimental toxicity is not proof of pediatric clinical safety or a reason to label every alternative as more toxic. SRI-41315 structural findings should not be assigned to ataluren.

## Ataluren literature

7. Ng MY, Li H, Ghelfi MD, Goldman YE, Cooperman BS. Ataluren and aminoglycosides stimulate read-through of nonsense codons by orthogonal mechanisms. *Proceedings of the National Academy of Sciences*. 2021;118(2):e2020599118. DOI: [10.1073/pnas.2020599118](https://doi.org/10.1073/pnas.2020599118). PMID: 33414181.
8. Bezzerri V, Bardelli D, Morini J, et al. Ataluren-driven restoration of Shwachman-Bodian-Diamond syndrome protein function in Shwachman-Diamond syndrome bone marrow cells. *American Journal of Hematology*. 2018;93(4):527–536. DOI: [10.1002/ajh.25025](https://doi.org/10.1002/ajh.25025). PMID: 29285795.
9. Bezzerri V, Pegoraro A, Hristodor AM, et al. Ataluren improves hematopoietic and pancreatic disorders in Shwachman-Diamond syndrome patients: a compassionate program case-series. *Nature Communications*. 2025;16:8189. DOI: [10.1038/s41467-025-63137-3](https://doi.org/10.1038/s41467-025-63137-3). PMID: 40897730.
10. Vlodavets D, Wu S, Kostera-Pruszczyk A, et al. Confirmatory long-term efficacy and safety results of ataluren in patients with nmDMD from Study 041, an international, randomized, double-blind, placebo-controlled, Phase III trial. *Journal of Comparative Effectiveness Research*. 2025;14:e240238. DOI: [10.57264/cer-2024-0238](https://doi.org/10.57264/cer-2024-0238). PMID: 40843507.
11. Kerem E, Konstan MW, De Boeck K, et al. Ataluren for the treatment of nonsense-mutation cystic fibrosis: a randomised, double-blind, placebo-controlled phase 3 trial. *The Lancet Respiratory Medicine*. 2014;2(7):539–547. DOI: [10.1016/S2213-2600(14)70100-6](https://doi.org/10.1016/S2213-2600(14)70100-6). PMID: 24836205.
12. Konstan MW, VanDevanter DR, Rowe SM, et al. Efficacy and safety of ataluren in patients with nonsense-mutation cystic fibrosis not receiving chronic inhaled aminoglycosides: The international, randomized, double-blind, placebo-controlled Ataluren Confirmatory Trial in Cystic Fibrosis (ACT CF). *Journal of Cystic Fibrosis*. 2020;19(4):595–601. DOI: [10.1016/j.jcf.2020.01.007](https://doi.org/10.1016/j.jcf.2020.01.007). PMID: 31983658.

The 15 claims in the supplied augmented ledger come from six publication IDs; they are not 15 independent studies. Source-anchor matching is not completed scientific adjudication. Preserve unfavorable findings and analysis hierarchy. Study 041's favorable ITT analyses do not replace its nonsignificant primary mITT result. Cell findings and a small uncontrolled SDS series do not establish benefit in another disorder.

## Regulatory sources (checked 23 October 2026)

- [UK Translarna SmPC](https://www.medicines.org.uk/emc/product/2576/smpc), updated 31 March 2026: conditional authorization for the specified ambulatory nonsense-mutation Duchenne population aged 2 years and older; not an authorization for another indication.
- [EMA Translarna EPAR](https://www.ema.europa.eu/en/medicines/human/EPAR/translarna): EU authorization not renewed on 28 March 2025 because effectiveness was not confirmed.
- [PTC regulatory update, 12 February 2026](https://ir.ptcbio.com/news-releases/news-release-details/ptc-therapeutics-provides-regulatory-update-translarnatm-2): the US NDA resubmission was withdrawn after FDA feedback; ataluren had not been approved for this indication in the US.

Authorization, availability, access, safety and case suitability are different questions. No treatment-initiation or dosing recommendation is made.

## Public annotation, phenotype and network resources

| Resource | Citation or official guidance | Reproduction detail to retain |
|---|---|---|
| Ensembl VEP | McLaren W, Gil L, Hunt SE, et al. *The Ensembl Variant Effect Predictor*. *Genome Biology*. 2016;17:122. DOI: [10.1186/s13059-016-0974-4](https://doi.org/10.1186/s13059-016-0974-4). PMID: 27268795. | VEP/cache version, reference build, flags and plugins. |
| HPO | Gargano MA, Matentzoglu N, Coleman B, et al. *The Human Phenotype Ontology in 2024: phenotypes around the world*. *Nucleic Acids Research*. 2024;52(D1):D1333–D1346. DOI: [10.1093/nar/gkad1005](https://doi.org/10.1093/nar/gkad1005). | Ontology and gene-annotation snapshot. |
| MANE | Morales J, Pujar S, Loveland JE, et al. *A joint NCBI and EMBL-EBI transcript set for clinical genomics and research*. *Nature*. 2022;604:310–315. DOI: [10.1038/s41586-022-04558-8](https://doi.org/10.1038/s41586-022-04558-8). | MANE release and exact transcript versions. |
| Phenopacket Store | Danis D, et al. *A corpus of GA4GH phenopackets: Case-level phenotyping for genomic diagnostics and discovery*. *HGG Advances*. 2025;6(1):100371. DOI: [10.1016/j.xhgg.2024.100371](https://doi.org/10.1016/j.xhgg.2024.100371). | Actual release, inclusion policy, grouping and split audit. |
| STRING | [Official citation guidance](https://string-db.org/cgi/about); Szklarczyk D, Kirsch R, Koutrouli M, et al. *The STRING database in 2023*. *NAR*. 2023;51(D1):D638–D646. DOI: [10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000); 2025 update DOI [10.1093/nar/gkae1113](https://doi.org/10.1093/nar/gkae1113). | Exact release, network file/type and confidence threshold. A citation does not assert every STRING feature was used. |
| Reactome | [Official citation instructions](https://reactome.org/cite). | Release, mappings and stable pathway IDs; cite pathway-specific source where appropriate. |
| UniProt | [Official publication guidance](https://www.uniprot.org/help/publications). | Release, reviewed-status policy and mapping file. |
| SIGNOR | Lo Surdo P, Iannuccelli M, Contino S, et al. *SIGNOR 3.0, the SIGnaling network open resource 3.0: 2022 update*. *NAR*. 2023;51(D1):D631–D637. DOI: [10.1093/nar/gkac883](https://doi.org/10.1093/nar/gkac883). | Actual snapshot/schema, species, signs and evidence mapping. The paper is not proof of snapshot version. |
| ChEMBL | Zdrazil B, Felix E, Hunter F, et al. *The ChEMBL Database in 2023*. *NAR*. 2024;52(D1):D1180–D1192. DOI: [10.1093/nar/gkad1004](https://doi.org/10.1093/nar/gkad1004). Web-service paper: Davies M, Nowotka M, Papadatos G, et al. *NAR*. 2015;43(W1):W612–W620. DOI: [10.1093/nar/gkv352](https://doi.org/10.1093/nar/gkv352). | Actual release, native acquisition, action normalization and applicable data license. |
| gnomAD | [Resource](https://gnomad.broadinstitute.org/) and [maintainer citation guidance](https://github.com/broadinstitute/gnomad-browser/issues/1832). | Exact release, exome/genome subset and access date; do not substitute a current version for a historical lookup. |
| ClinVar | Landrum MJ, Lee JM, Benson M, et al. *ClinVar: public archive of relationships among sequence variation and human phenotype*. *NAR*. 2014;42(D1):D980–D985. DOI: [10.1093/nar/gkt1113](https://doi.org/10.1093/nar/gkt1113); [official citation guidance](https://www.ncbi.nlm.nih.gov/clinvar/docs/faq/). | Snapshot/query date and accession/version for records used. |
| REVEL | Ioannidis NM, Rothstein JH, Pejaver V, et al. *REVEL: An Ensemble Method for Predicting the Pathogenicity of Rare Missense Variants*. *AJHG*. 2016;99(4):877–885. DOI: [10.1016/j.ajhg.2016.08.016](https://doi.org/10.1016/j.ajhg.2016.08.016). | Score release and retrieval source. |
| SIFT | Ng PC, Henikoff S. *SIFT: Predicting amino acid changes that affect protein function*. *NAR*. 2003;31(13):3812–3814. DOI: [10.1093/nar/gkg509](https://doi.org/10.1093/nar/gkg509). | Annotation/cache version; correlated predictors are not independent evidence. |
| PolyPhen-2 | Adzhubei IA, Schmidt S, Peshkin L, et al. *A method and server for predicting damaging missense mutations*. *Nature Methods*. 2010;7:248–249. DOI: [10.1038/nmeth0410-248](https://doi.org/10.1038/nmeth0410-248). | Score model and annotation/cache source. |

These citations attribute methods/resources; they do not prove the exact historical release used or permission to redistribute the underlying data.

## Software and environment citations

- R Core Team. (2024). *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing, Vienna. [Project](https://www.R-project.org/). The supplied no-training launcher reported R 4.3.3; this is not a verified historical training environment.
- Moelder F, Jablonski KP, Letcher B, et al. *Sustainable data analysis with Snakemake*. *F1000Research*. 2021;10:33. Version 3, 2025. DOI: [10.12688/f1000research.29032.3](https://doi.org/10.12688/f1000research.29032.3). Snakemake 9.27.1 was observed in the supplied Python environment; planning a DAG is not execution of its scientific jobs.
- Danecek P, Bonfield JK, Liddle J, et al. *Twelve years of SAMtools and BCFtools*. *GigaScience*. 2021;10(2):giab008. DOI: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008). Supplied CLI probes observed bcftools 1.19; cite tools only for stages that actually used them.
- Bonfield JK, Marshall J, Danecek P, et al. *HTSlib: C library for reading/writing high-throughput sequencing data*. *GigaScience*. 2021;10(2):giab007. DOI: [10.1093/gigascience/giab007](https://doi.org/10.1093/gigascience/giab007). Supplied probes observed tabix/bgzip 1.19.
- Carpenter B, Gelman A, Hoffman MD, et al. *Stan: A Probabilistic Programming Language*. *JSS*. 2017;76(1):1–32. DOI: [10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01).
- Gabry J, Češnovar R, Johnson A, Bronder S. *cmdstanr: R Interface to CmdStan*. [Official package citation](https://mc-stan.org/cmdstanr/authors.html). The supplied no-training R probe did not resolve cmdstanr; its use/version must be stated from the historical training environment.
- Bürkner PC, Gabry J, Kay M, Vehtari A. *posterior: Tools for Working with Posterior Distributions*. [Official package citation metadata](https://cloud.r-project.org/web/packages/posterior/citation.html). The supplied no-training R probe did not resolve posterior; do not assign a current web version to historical fits.
- Csardi G, Nepusz T. (2006). *The igraph software package for complex network research*. *InterJournal*, Complex Systems:1695. The supplied runtime probe observed igraph 1.6.0.
- Ooms J. (2014). *The jsonlite Package: A Practical and Consistent Mapping Between JSON Objects and R Classes*. *arXiv:1403.2805*.

### Direct/support R packages recorded in the supplied citation export

| Package | Version observed | Citation note |
|---|---:|---|
| data.table | 1.14.10 | Barrett T, Dowle M, Srinivasan A. *data.table: Extension of data.frame*. Use installed `citation("data.table")` output. |
| yaml | 2.3.8 | Garbett SP, Stephens J, Simonov K, et al. *yaml: Methods to Convert R Data to YAML and Back*. Use installed citation output. |
| jsonlite | 1.8.8 | Ooms 2014 paper above. |
| digest | 0.6.34 | Eddelbuettel D, et al. *digest: Create Compact Hash Digests of R Objects*. Use installed citation output. |
| igraph | 1.6.0 | Csardi & Nepusz 2006; also use release-specific package citation if generated. |
| Matrix | 1.6.5 in version probe; `1.6-5` in lock | Bates D, Maechler M, Jagan M. *Matrix: Sparse and Dense Matrix Classes and Methods*. Reconcile package-version spelling to the exact run’s `packageVersion()` output. |
| renv | 1.0.3 | Ushey K, Wickham H. *renv: Project Environments*. Use installed citation output. |
| R6 | 2.5.1 in lock | Record/cite if directly used. |
| loo | 2.6.0 in lock | Installation does not establish that LOO/stacking was executed in a given run. |
| matrixStats | 1.5.0 in lock | Record/cite if directly used. |

The committed no-training `renv` lock also contains transitive packages (`backports`, `checkmate`, `cli`, `cpp11`, `lattice`, `lifecycle`, `magrittr`, `pkgconfig`, `rlang`). The Python requirements file records Snakemake/PyYAML/requests and transitive dependencies. These files are complete environment records, while this table is not a citation for every transitive dependency.

### Other observed tools

The participant supplied version output for Python 3.14.8, PyYAML 6.0.3, requests 2.34.2, Bash 5.2.21, Git 2.43.0, gzip 1.12, curl 8.5.0, bcftools 1.19, tabix/bgzip 1.19 and coreutils/sha256sum 9.4. These are observations from the current Snakemake environment/PATH, not proof every tool contributed to every run. VEP and CmdStan historical versions were not established by that probe.

## Capture and interpretation boundaries

Version observations are transcribed from participant-supplied terminal output and the committed no-training lock; the assistant did not probe the user’s machine. The no-training launcher’s failure to resolve cmdstanr/posterior does not establish absence from a different training environment. Current online package pages are not historical version records.

Citations do not establish data licenses, exact local snapshots, feature usage, scientific adjudication, biological causality, functional rescue, clinical benefit or safety. Keep private configs, genotype-scale outputs, case manifests and source-root mappings outside public documentation. Do not upload resource/model artifacts merely to provide citations.

## Acknowledgement

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.

See [AI assistance disclosure](AI_ASSISTANCE.md). The formal challenge-dataset publication citation is separate from these resource/software citations; do not invent a DOI when the organizer-prescribed reference is unavailable.
