# Rare-Disease Evidence Pipeline

An evidence-preserving research framework for rare-disease candidate prioritization, externally trained Bayesian ranking, and mechanism-aware drug-repurposing exploration.

I separate case evidence, external model development, frozen inference, and downstream research hypotheses so that each result can be traced to its inputs, assumptions and limitations. The aim is not just to produce a ranked list, but to explain what supports it, what is missing, and how it changes under alternative assumptions.

## Repository status and navigation

The CLI provides `bin/mva` for Track 1, Track 2 and independent candidate-gene exploration. Start with [How to run](#how-to-run), or see the [CLI guide](docs/CLI.md) for detailed inputs and historical handoff preparation.

This README describes the workflow, its mathematical contracts and the steps for running the configured analysis.

- `main` currently contains the portfolio overview and development-policy documentation.
- The [external-training branch](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training) contains the imported Stan definitions, retained-fit diagnostic material and model documentation. These files are not merged by this README update.
- [Model documentation](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/models) explains the four fitted specifications.
- [Stan definitions](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/workflow/stan) are the executable reference for the equations below.
- [Environment records](https://github.com/LukefromWA/rare-disease-evidence-pipeline/tree/external-training/environment) record the imported tooling environment.
- [AI assistance disclosure](docs/AI_ASSISTANCE.md) records the current disclosure.

The configured Linux/WSL workflow has completed end to end, including Track 1, historical candidate handoff and Track 2, with output checks and 15 passing CLI tests. Setup requires the runtime environment, model artifacts and resource bundles described below. Optional methods record their availability explicitly.

## How to run

### 1. Prepare the environment and resources

Use Linux or WSL with the configured Python, R, Snakemake and supporting
command-line tools available. Activate the runtime environment and provision
the resource/model/analysis profile, frozen model artifacts, core resource
bundle and matching Reactome hierarchy.

From the repository root:

```bash
./bin/mva --help
```

Set `MVA_PYTHON` if the launcher needs an explicit Python interpreter.
Keep the matching Snakemake and supporting tools on `PATH`.

Setup time depends on which dependencies and resource bundles are already
available. The execution times below exclude installation, resource preparation
and raw VCF annotation.

### 2. Run Track 1

Start the interactive menu:

```bash
./bin/mva
```

Choose `1`, then enter:

1. The number of cores.
2. The complete profile and project/data root.
3. Confirmation of the profile's QC thresholds, priors and mechanism assumptions.
4. The variant inputs.
5. The selected sample label and observed HPO text file.
6. A new private run directory outside the repository.

For an existing VEP-annotated VCF, answer `n` to raw input and `y` to existing
VEP output. Supply these files at their respective prompts:

| Prompt | Input |
|---|---|
| Existing VEP-annotated VCF | Annotated `.vcf.gz` |
| Annotated VCF index | Matching `.vcf.gz.tbi` |
| Original VEP command/provenance | Original command or provenance text file |
| Observed HPO text file | One `HP:nnnnnnn` term per line |

The wizard saves `run.yaml`, a sample sheet and phenotype snapshots in the
new workspace. To execute, answer `n` to “Plan only” and confirm execution.

A saved configuration can also be run directly:

```bash
./bin/mva track1 \
  --config /path/to/workspace/run.yaml \
  --root /path/to/project-data \
  --cores 1
```

Add `--dry-run` to preview the job plan without executing the analysis.

Estimated runtime: about 22 minutes for the observed single-core,
existing-VEP run. Raw annotation time is additional and was not measured
in this run.

### 3. Prepare the historical candidate handoff

The current connected route combines freshly scored candidates with explicitly
configured same-case historical C4/A10 candidates.

Follow [Historical handoff preparation](docs/CLI.md#existing-case-historical-handoff-preparation)
to verify the historical manifest and input checksums, save a separate
`track2_handoff.yaml`, and build the new C4 union using the completed Track 1
artifacts.

Track 1 alone does not build this historical union. The existing-run Track 2
command expects it to be prepared first.

Estimated runtime: about 5 seconds for the observed handoff build after
configuration and verification.

### 4. Run Track 2

Use the prepared handoff configuration:

```bash
./bin/mva track2 \
  --config /path/to/workspace/track2_handoff.yaml \
  --root /path/to/project-data \
  --cores 1
```

Add `--dry-run` to preview execution.

The interactive existing-run route is also available through menu choice `2`,
then `1`. It resolves `run.yaml` before legacy `config.yaml`. Use `--config`
when selecting a separately prepared handoff configuration.

Estimated runtime: about 1 minute for the observed core-resource run.
Available methods, candidate counts, resource sizes and live-query latency
affect runtime.

For independent candidate-gene exploration instead of weighted case analysis:

```bash
./bin/mva track2 --genes TP53 BRCA1 \
  --resources /path/to/validated-core-resources \
  --hierarchy /path/to/matching-reactome-hierarchy \
  --out /path/to/new-gene-run \
  --cores 1 --dry-run
```

Remove `--dry-run` to execute. The gene-only route produces unweighted
exploration; it does not create variant evidence, Stan probabilities or
mechanism weights.

### 5. Inspect the results

Case outputs are stored under the configured output root and sample label:

| Output | Location |
|---|---|
| Track 1 ranking | `a10/ranked_candidates.tsv.gz` |
| Frozen prediction summary | `c3_stan_cv/stan_cv_predictions_summary.tsv` |
| Ranking-sensitivity results | `a11_review/calculation/` |
| Evidence-method status | `d3/*status.tsv` |
| Held/excluded scenario-drug rows | `d5/held_or_excluded.tsv` |
| Sensitivity and scenario coverage | `d6/status.tsv` |
| Track 2 final report | `d7/final_report.tsv` |
| Execution logs | `logs/` |

Optional methods report unavailable, disabled or not-applicable status
explicitly. A report with no eligible drugs is a valid result; the eligibility
audit records why rows were held or excluded. Patient HPO inputs and the
optional drug-indication HPO resource serve different purposes.

Outputs support computational research and human review, not treatment
recommendations. Keep private case inputs, configurations and generated results
outside the public repository.

### 6. End-to-end run and reproducibility

The configured Linux/WSL run completed Track 1, the historical candidate
handoff and Track 2. Output checks covered the final report, method-status
files, eligibility audit and scenario coverage. All 15 existing CLI tests passed.

| Stage | Completed steps | Observed elapsed time |
|---|---:|---:|
| Track 1 | 22/22 | 22 minutes, 7 seconds |
| Historical handoff | 2/2 | About 5 seconds |
| Track 2 | 17/17 | 56 seconds |

Allow approximately 23 minutes of execution time for a comparable configured
single-core run, plus setup and manual handoff preparation. These are observed
timings, not benchmarks for other machines or datasets.

Run the CLI tests with the configured Python:

```bash
python -m unittest discover -s tests -p 'test_cli*.py' -v
```

The tested path used existing VEP output, prepared local resources and an
explicit historical handoff. Raw VCF-to-VEP execution and optional methods
require their own dependencies and checks.

## Workflow overview

A processes case evidence. B develops models using external public training data and does not train on the hackathon child's genome. C applies the frozen B artifacts to case features without refitting. Track 2 receives the union of the retained A and C candidate outputs, preserving their distinct evidence views. D evaluates mechanism-aware research hypotheses; E supplies validated reference resources; F reports the results and provenance.

```mermaid
flowchart TD
    PP[External public Phenopackets] --> B[B: external cohort preparation and training]
    B --> P[Retained M3b predictor: three equal-weight components]
    IN[Case variant and phenotype inputs] --> A[A: evidence preparation and ranking]
    A --> DEP[Dependence and ranking-sensitivity analysis]
    A --> C[C: frozen inference]
    P --> C
    A --> U[Candidate union: A and C]
    C --> U
    U --> D[D: Track 2 mechanism-aware research]
    DEP --> D
    E[E: versioned reference resources] --> D
    D --> F[F: reporting and provenance]
```

The union is a candidate-registry operation, not an addition of scores:

```math
\mathcal H_{\mathrm{Track\ 2}}=\mathcal H_A\cup\mathcal H_C.
```

Align candidate identifiers and preserve source membership, underlying genomic evidence, coverage and the separate A/C scores. Do not silently drop A-only candidates or treat different score types as interchangeable. The diagram is an overview, not an exhaustive list of resource dependencies.

## A — Case evidence preparation and ranking

### A0–A3: annotation, evidence extraction and QC

Inputs are authorized variant data, annotation records, phenotype information and configuration. Normalize alleles and genome-build conventions, read annotation schemas from their headers, and join transcript consequences using exact allele/transcript identifiers rather than gene symbols alone.

Retain individual QC flags and exclusion reasons. Thresholds belong in configuration and run manifests, not hidden constants. A shortlist should preserve why a row passed or failed and distinguish transcript-level annotations from candidate-level records.

For allele-fraction analysis, use informative allelic depth consistently:

```math
n_i=AD_{\mathrm{ref},i}+AD_{\mathrm{alt},i},\qquad
\widehat f_i=\frac{AD_{\mathrm{alt},i}}{n_i},\quad n_i>0.
```

Total depth and informative allelic depth need not be identical. Missing annotation or population-frequency information must remain distinguishable from an observed negative result.

Output: a keyed candidate evidence table, transcript provenance, QC flags and an exclusion audit.

### A4: phenotype similarity

Compare the case's observed phenotype terms with gene-associated term sets using a declared ontology and annotation snapshot. Information content gives rarer terms greater specificity:

```math
IC(t)=-\log\frac{n_t}{N},\qquad
s_{\mathrm{Lin}}(a,b)=\frac{2IC(\mathrm{MICA}(a,b))}{IC(a)+IC(b)}.
```

Here, n_t counts annotated entities associated with a term or its descendants in the declared universe. MICA is the most informative common ancestor. Edge cases such as a zero denominator need an explicit implementation convention.

Aggregate term similarities with symmetric best-match average:

```math
\mathrm{BMA}(P,G)=\frac12\left[
\frac1{|P|}\sum_{a\in P}\max_{b\in G}s(a,b)
+\frac1{|G|}\sum_{b\in G}\max_{a\in P}s(a,b)
\right].
```

P is the observed case-term set and G is the gene annotation set. Resnik, depth-based Wu–Palmer and ancestor-set Jaccard are complementary views, not independent confirmations. Missing gene annotations are unscorable coverage, not evidence against the gene.

Output: gene-level similarity views, ranks, annotation coverage and resource provenance.

### A5: network propagation

Use phenotype-derived seeds and, in a separately labeled informed view, a declared prior gene panel. Keep the agnostic view free of that panel. A current top-ranked candidate must not become an undeclared prior seed.

For a column-normalized transition matrix P, seed distribution s, restart probability r and dangling-column indicator d, the build-guide update is:

```math
p_{t+1}=r s+(1-r)\left(Pp_t+s\,d^{\mathsf T}p_t\right).
```

The dangling mass returns to the seed distribution. This equation uses r as restart probability, not damping probability. Record the graph, seed policy and implementation used by each run; corrections must not silently replace historical outputs.

A multiplex view connects network layers while keeping their identities available for comparison. For two layers, one possible declared transition structure is:

```math
T=\begin{pmatrix}
(1-\lambda)P_1&\lambda I\\
\lambda I&(1-\lambda)P_2
\end{pmatrix}.
```

Normalization, common node ordering and dangling handling must be checked. If no valid seeds remain, return an unavailable channel rather than a fabricated zero score.

Output: informed/agnostic network scores, graph coverage and propagation settings.

### A6–A7: multiview comparison and coverage-aware integration

Align genomic, phenotype and network evidence by candidate identity. Agreement statistics describe relationships between views; an integrated rank is not an independent input view.

For channel score s_ik, declared weight w_k and availability M_ik, the available-channel mean is:

```math
\overline S_i=\frac{\sum_k w_k M_{ik}s_{ik}}{\sum_k w_k M_{ik}},
\qquad n_{\mathrm{views},i}=\sum_k M_{ik}.
```

The score is undefined when the denominator is zero. Missing evidence is not replaced with biological negative evidence. Score transformations and percentile universes must be documented before combination.

A declared prior-weighted normalization may then be written:

```math
Q_i=\frac{\alpha_i L_i}{\sum_j\alpha_j L_j}.
```

The meaning of L_i must be stated: a transformed evidence score is not automatically a calibrated likelihood. Report informed and agnostic versions separately. Normalization alone does not make Q_i a clinical posterior probability.

### Dependence and channel attribution

The dependence analysis examines whether phenotype and network evidence overlap. A rank-based Gaussian transform on paired available observations gives:

```math
z_i=\Phi^{-1}\left(\frac{\mathrm{rank}(s_i)}{n+1}\right),
\qquad \widehat\rho=\mathrm{cor}(z^{(H)},z^{(N)}).
```

Declare tie handling and the complete-case sample size. A modest correlation does not prove independence or validate a chosen copula.

Exact Shapley attribution describes each channel's contribution to a declared scoring function v over the available evidence channels K:

```math
\phi_k=\sum_{S\subseteq K\setminus\{k\}}
\frac{|S|!(|K|-|S|-1)!}{|K|!}
\left[v(S\cup\{k\})-v(S)\right].
```

Specify the empty-set and missing-data conventions and check that contributions sum to v(K) minus v(empty). Attribution explains the score function; it does not establish biological causality.

### A8: allele-fraction and mosaicism review

With prior Beta(a_0,b_0) and informative counts:

```math
f_i\mid\mathrm{counts}\sim
\mathrm{Beta}(a_0+AD_{\mathrm{alt},i},b_0+AD_{\mathrm{ref},i}).
```

Report credible intervals and read-level/QC qualifications. A beta-binomial versus fixed heterozygous-binomial comparison is a model comparison, not proof of mosaicism. A low fraction alone does not identify its biological cause.

### A9–A10: robustness summaries and evidence-preserving handoff

Join coverage, allele-fraction review, informed/agnostic rankings and robustness summaries into the retained candidate registry. A gene-level summary must preserve its variant links and aggregation rule; it must not erase distinct alleles or compound-heterozygous hypotheses.

Output: ranked candidates with evidence coverage, assumptions, stability summaries and provenance. Case-derived tables remain inside the authorized analysis boundary unless specifically cleared for the intended export.

## A11 — Ranking sensitivity under declared assumptions

A11 keeps sensitivity results separate from the original evidence scores. Perturb declared weights or assumptions, rerank, and record how often each candidate remains in the top m:

```math
\widehat\pi_i^{(m)}=\frac1B\sum_{b=1}^{B}
\mathbf 1\{\mathrm{rank}_b(i)\le m\}.
```

For two top-m sets with overlap r from n candidates, the Kuncheva index is:

```math
KI=\frac{rn-m^2}{m(n-m)},\qquad 0<m<n.
```

Leave-one-channel-out comparisons measure score and rank changes. Record the perturbation distribution, seed, candidate universe and tie rules. These frequencies describe robustness to the chosen perturbations, not the probability that a candidate is causal.

The build guide also groups stability/influence calculations with the A9 robustness stage; A11 is retained here as the separate sensitivity view in the current repository overview.

## B — External-cohort preparation and model development

B uses external public training data. It does not read or train on the hackathon child's genomic data. The retained model-fold artifacts come from this external-training branch of the analysis.

### B1–B7: eligibility, labels, group audits and deduplication

Screen and parse source records, distinguish observed from explicitly excluded phenotype terms, normalize gene labels, assess resource mapping, audit repeated case groups, select deterministic representatives, and canonicalize positive labels.

The deduplication contract must record the actual grouping fields, conflict policy, representative-selection rule and exclusion reasons. Publication or phenotype resemblance can require review; a packet identifier alone is not evidence that two records are independent.

Output: external cases, canonical labels, duplicate-group identifiers and auditable inclusion/exclusion records. This README does not invent duplicate counts or certify a no-overlap result without the corresponding audit.

### B8 → B10 → B9: candidates, splits, then features

1. B8 constructs candidate sets with a declared negative/background policy and positive-label contract.
2. B10 assigns whole duplicate groups/cases to splits or held-out folds and checks disjointness.
3. B9 fits transformations on the appropriate training partition only, then applies them to held-out rows.

For any held-out partition, require disjoint case and duplicate-group identities from its training partition. Record the actual protocol: a grouped OOF design is not the same as a fixed train/validation/test design.

The guide's gene-frequency feature excludes the training case's own label contribution:

```math
G_{gc}=\log\left(1+\sum_{c'\in\mathcal T,\,c'\ne c}
\mathbf 1\{\mathrm{causal}(c')=g\}\right).
```

Phenotype-overlap, network and interaction features must follow their recorded definitions. A gene-frequency feature is not a raw genomic-evidence channel. Stan accepts X as input; the feature meanings are supplied by the matching preprocessing schema, not inferred from model names.

### B11–B13: four Stan specifications

Each training case has one labeled candidate position within its supplied candidate set I_c. All four models use the same conditional categorical-logit likelihood:

```math
p_{ci}(\theta)=\frac{\exp(\eta_i)}{\sum_{j\in I_c}\exp(\eta_j)},
\qquad
\log L(\theta)=\sum_c\left[\eta_{y_c}-\log\sum_{j\in I_c}e^{\eta_j}\right].
```

Stan's `y_pos` is the local position within a contiguous case segment. These probabilities depend on which candidates were supplied; they are not clinical pathogenicity probabilities.

#### M2a: global feature coefficients

```math
\eta_i=x_i^{\mathsf T}\beta,\qquad
\beta_k\sim\mathcal N(0,s_{\beta,k}^2).
```

M2a has no added copula score. Its declared group and pseudo-observation inputs do not enter this predictor. Omitting a dependence term does not prove empirical independence.

#### M2b: global coefficients with a global copula contribution

Let z_hi and z_ni be inverse-Normal transforms of the supplied HPO/network pseudo-observations. The implemented Gaussian-copula log density is:

```math
\ell_i(\rho)=-\frac12\log(1-\rho^2)
-\frac{\rho^2(z_{hi}^2+z_{ni}^2)-2\rho z_{hi}z_{ni}}{2(1-\rho^2)}.
```

With availability mask m_i and gamma denoting `beta_cop`:

```math
\rho=\tanh(r),\qquad
\eta_i=x_i^{\mathsf T}\beta+\gamma m_i\ell_i(\rho),
```
```math
r\sim\mathcal N(0,s_{\mu\rho}^2),\qquad
\gamma\sim\mathcal N(0,s_\gamma^2),\qquad
\beta_k\sim\mathcal N(0,s_{\beta,k}^2).
```

The raw correlation parameter is `rho_raw`. Correlation and copula utility are global. The utility coefficient is not constrained to be nonnegative.

#### M3a: partially pooled group coefficients

```math
\beta_{gk}=\beta_k+\tau_{\beta,k}z_{\beta,gk},\qquad
\eta_i=\sum_k x_{ik}\beta_{g(i),k},
```
```math
z_{\beta,gk}\sim\mathcal N(0,1),\quad
\beta_k\sim\mathcal N(0,s_{\beta,k}^2),\quad
\tau_{\beta,k}\sim\mathrm{HalfNormal}(s_{\tau\beta,k}).
```

The grouping factor must be declared in the run metadata. M3a adds no copula contribution. For the unseen-group slot `G + 1`, generated quantities sample a population-derived coefficient vector; rows assigned to that slot share it within a posterior draw.

#### M3b: partially pooled coefficients, dependence and copula utility

M3b retains the M3a coefficient hierarchy and adds:

```math
\rho_g=\tanh(\mu_\rho+\tau_\rho z_{\rho,g}),\qquad
\gamma_g=\gamma+\tau_\gamma z_{\gamma,g},
```
```math
\eta_i=\sum_k x_{ik}\beta_{g(i),k}
+\gamma_{g(i)}m_i\ell_i(\rho_{g(i)}).
```

The added priors are:

```math
z_{\rho,g},z_{\gamma,g}\overset{\mathrm{ind}}\sim\mathcal N(0,1),\quad
\mu_\rho\sim\mathcal N(0,s_{\mu\rho}^2),\quad
\gamma\sim\mathcal N(0,s_\gamma^2),
```
```math
\tau_\rho\sim\mathrm{HalfNormal}(s_{\tau\rho}),\qquad
\tau_\gamma\sim\mathrm{HalfNormal}(s_{\tau\gamma}).
```

Here, `mu_rho_raw` is the unconstrained location before the tanh transform; gamma and tau_gamma correspond to `beta_cop` and `tau_cop`. The Stan definitions supply the prior scales as data and include the half-Normal normalization constants. New-group coefficients, correlation and utility are sampled from the fitted populations rather than silently assigned to an existing group.

### Evaluation and optional calibration

Keep held-out candidate recovery, ranking/log-score evaluation and posterior sampling diagnostics separate. Model selection and transformation fitting must not use the held-out evaluation labels. If conformal calibration is enabled, document its calibration set, nonconformity definition and threshold; do not claim coverage without checking the required calibration/selection assumptions and actual evaluation.

R-hat, effective sample size, divergences and treedepth answer different sampling questions. The retained-fit documentation reports divergence qualifications separately and leaves treedepth unverified. Passing retained-parameter thresholds is not evidence of predictive accuracy or complete sampler clearance.

### B14: freezing and the retained predictor

Freeze the feature schema, transformations, group mapping, model definitions, priors/configuration, artifact identifiers and checksums. The documented M3b deployment adapter averages three retained fold-specific predictors equally:

```math
\widehat p_i^{\mathrm{frozen}}=
\frac13\sum_{f=1}^{3}\mathbb E_{\theta_f\mid D_{-f}}[p_i(\theta_f)].
```

Average probabilities after within-case softmax, not coefficient vectors before softmax. This is neither one full-data posterior nor evidence-weighted Bayesian model averaging across the four specifications. An OOF case uses its excluding component, not the new-case mixture of all components.

## C — Frozen inference and candidate comparison

### C1: reproduce the external feature contract

Apply the same feature function and frozen training state used by B9:

```math
X_{B9}=f(\mathrm{external\ candidates};\widehat\psi_{\mathrm{train}}),\qquad
X_{C1}=f(\mathrm{case\ candidates};\widehat\psi_{\mathrm{train}}).
```

The input identities differ; the feature meanings, column order and fitted state must not. Similar column names do not make a phenotype BMA score interchangeable with an overlap feature, or an A integrated score interchangeable with a training gene-frequency prior.

### C2–C3: keep model views distinct

Compare private-case baseline views with frozen external predictions without refitting B. Preserve candidate-set definitions, missingness, posterior uncertainty and unseen-group handling. Candidate-set sizes and construction affect softmax interpretation.

### C4: retain the hypothesis registry

Align views by stable candidate/hypothesis identifiers and retain the declared weighted set, not just rank one. Record the rule producing nonnegative hypothesis weights q_j that sum to one; do not invent combination weights in documentation.

For Track 2, combine the retained A and C registries by their union, preserving source membership and evidence. Duplicate identities need explicit reconciliation, not repeated counting. Missing C predictions are not negative evidence against an A-only hypothesis. The q_j assignment across the union must be declared and checked before downstream scenario weighting.

## D — Mechanism-aware Track 2 exploration

### D0–D2: uncertainty handoff, mechanisms and targets

Expand each retained hypothesis into its applicable mechanism scenarios. Join mechanisms to the matching hypothesis, not an unrestricted cross-join. With declared mechanism prior p0(m|j):

```math
w_{jm}=\frac{q_j\,p_0(m\mid j)}
{\sum_{j',m'}q_{j'}p_0(m'\mid j')}.
```

Keep the declared prior separate from evidence-updated mechanism assessments. Competing loss-of-function, gain-of-function, dominant-negative and unknown states may imply different intervention directions. Every removal or renormalization needs an audit reason.

Map genes to protein/target identifiers using versioned resources, preserving ambiguous mappings and action direction rather than guessing.

### D3.1–D3.8: complementary research methods

1. Target/action queries and direction compatibility.
2. Signature reversal when a suitable expression resource is available.
3. Pathway-based hypotheses using actual pathway evidence.
4. Phenotype/indication similarity when suitable annotations exist.
5. Literature-supported hypotheses with traceable sources.
6. Signed network propagation with uncertainty analysis.
7. Controllability analyses only with a declared model and interpretation.
8. Rescue optimization using its declared upstream method outputs.

Methods return evidence with provenance, or an explicit unavailable/disabled/not-applicable status. A missing data source is not replaced with toy evidence in production.

For pathway target sets T_d and hypothesis-associated pathway genes P_j, a declared overlap score is:

```math
J_{dj}=\frac{|T_d\cap P_j|}{|T_d\cup P_j|}.
```

A signed-propagation contract may use the source-to-target adjacency W:

```math
x_{t+1}=\tanh(\alpha W^{\mathsf T}x_t+x_0).
```

Edge orientation, signs, scale and perturbation definitions must match the implementation. Signed methods require their validation/enablement checks; a simulated rescue score is not evidence of treatment efficacy. Controllability is a separately qualified analysis, not an automatic clinical gate.

### D4–D5: evidence integration and noncompensable gates

Integrate declared evidence families with availability recorded:

```math
S_d=\frac{\sum_f v_f M_{df}e_{df}}{\sum_f v_f M_{df}}.
```

The denominator must be positive. Preserve the number and identity of available families. Do not count mechanism evidence twice by copying it into the pathway family.

Direction compatibility, safety and feasibility remain separate gates. A large aggregate score cannot compensate for a failed gate. Unknown or unconfigured pediatric safety is not a pass. Eligibility labels describe research review, not clinical suitability.

### D6–D7: scenario sensitivity and research reporting

Perturb declared scenario and evidence-family weights, recompute scores, and rank all drugs within each draw. For k selected in configuration:

```math
\widehat\pi_d^{(k)}=\frac1B\sum_{b=1}^{B}
\mathbf1\{\mathrm{rank}_b(d)\le k\}.
```

Report rank ranges, coverage, gating status and assumptions that drive the ranking. These frequencies measure stability under the chosen scenarios; they are not probabilities of benefit. Final outputs remain research hypotheses, not treatment recommendations.

## E — Case-independent reference resources

Prepare versioned ontology, gene/protein mappings, interaction networks, pathway sets, drug-action mappings and safety resources where licensed and configured. Record source, version/date, license, transformations and checksums.

Validate identifier joins, species, edge direction/sign, duplicate/conflicting records and resource readiness. Optional methods stay disabled if their required resource or validation sentinel is absent. A resource checksum establishes identity, not scientific appropriateness or permission to redistribute.

## F — Reporting, provenance and release boundaries

F converts saved results into tables, figures and reports without silently recomputing different scientific results. Each run should record inputs, configuration, software versions, seeds, schema, checksums, exclusions and method status.

Keep these quantities separately labeled:

| Quantity | Interpretation |
|---|---|
| Coverage-aware integrated score | Evidence summary under a declared scaling/weighting rule |
| Prior-weighted normalized score | Conditional on the declared score/likelihood map and prior rule |
| A11 top-m frequency | Robustness under the specified ranking perturbations |
| Stan p_new | Within-case probability over the supplied candidate set |
| Track 2 scenario weight | Declared hypothesis/mechanism weighting |
| Drug top-k frequency | Ranking sensitivity, not efficacy or safety |

Reproduction requires more than matching filenames: check identifiers, schemas, row counts, split/group separation, B9/C1 feature parity, frozen-model behavior and the declared validation outputs.

## Data boundaries and privacy status

This portfolio's public-facing export is restricted to reviewed code, documentation, schemas and synthetic/public examples with appropriate permissions. Keep controlled case inputs, private phenotype lists/configuration, genotype-scale derivatives, case prediction objects, candidate registries, generated case reports, credentials and private consolidation material outside that export unless specifically cleared for the intended destination.

External provenance is established for the retained B training fits; file-content and release reviews still remain separate checks. Repository privacy does not itself authorize uploading controlled data. Ignore rules do not remove previously tracked files or history.

The current local preflight inventories files/history and release assets but does not certify that opaque artifacts, images, credentials and every publication surface have been cleared. A clean working tree, successful push or matching checksum is not privacy approval.

Three export-boundary tests were previously reported as passing in the local development workspace. The tests and reviewed implementation are not yet imported into `main`; that result is not end-to-end scientific validation or complete privacy clearance.

## License

The repository license file is pending consolidation. Participant-submission licensing must follow the hackathon requirements; third-party materials retain their own licenses. Do not infer redistribution rights from public availability.

## AI assistance

I used Perplexity Pro primarily for documentation, debugging and code/workflow integration. I configured the relevant privacy and data-usage settings from the beginning to prevent private data from being used for model training, and checked those settings throughout the project. 


See the [AI assistance disclosure](docs/AI_ASSISTANCE.md).

Scientific literature, data-resource attribution and software citations are collected in the [citation register](docs/CITATIONS.md). The [observed software-version table](environment/software_versions.observed.tsv) distinguishes participant-observed no-training/runtime versions from unverified historical training versions.

## Intended use

Computational research and human review only. This is not a clinically validated diagnostic or treatment-selection system. Candidate ranks, allele-fraction models and repurposing hypotheses do not by themselves establish diagnosis, pathogenicity, mosaicism, causality, safety or treatment efficacy.
