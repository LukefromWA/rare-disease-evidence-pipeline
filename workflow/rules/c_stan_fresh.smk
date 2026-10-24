
# Fresh frozen-Stan scoring and handoff; no fitting rules.
import glob
import hashlib
import json

_CV_SOURCE = SCRIPTS
_CV_RUNNER = os.path.join(SCRIPTS, "..", "lib", "stan_runner.R")
_CV_STAN_DIR = os.path.join(SCRIPTS, "..", "stan")
_CV_BUNDLE = config["frozen_stan"]["bundle"]
_UNION_CODE = os.path.join(SCRIPTS, "c4_d0_union_stan_v1_2.R")
_UNION_LIBS = sorted(glob.glob(os.path.join(SCRIPTS, "..", "lib", "*.R")))
_UNION_OUT = OUT + "/{sample}/c4_stan_union"
_REVIEW_WF = os.path.dirname(SCRIPTS)
_A11_EVIDENCE = OUT + "/{sample}/a11_review/bayesian_candidate_evidence.csv"
_A11_CALC = OUT + "/{sample}/a11_review/calculation"

def _fresh_union_history(wc):
    directory = config.get("track2", {}).get(
        "frozen_union_dirs", {}
    ).get(wc.sample)
    if not directory:
        raise ValueError("Missing configured historical union for this sample")
    if not os.path.isabs(directory):
        directory = os.path.join(ROOT, directory)

    manifest = os.path.join(directory, "manifest.json")
    with open(manifest, encoding="utf-8") as stream:
        meta = json.load(stream)

    if meta.get("version") != "A10_C4_union_Stan_v1_2":
        raise ValueError("Unsupported historical union manifest")
    if (
        meta.get("current_build") != "GRCh38"
        or meta.get("legacy_build") != "GRCh38"
    ):
        raise ValueError("Historical union does not declare expected GRCh38 builds")

    result = {"manifest": manifest}
    for role in ("legacy_c4", "legacy_a10"):
        entries = [
            entry for entry in meta.get("inputs", [])
            if entry.get("role") == role
        ]
        if len(entries) != 1:
            raise ValueError("Historical manifest needs one input role: " + role)
        entry = entries[0]
        path = entry["path"]
        if not os.path.isabs(path):
            path = os.path.join(ROOT, path)
        with open(path, "rb") as stream:
            actual = hashlib.sha256(stream.read()).hexdigest()
        if actual != entry.get("sha256"):
            raise ValueError("Historical input checksum mismatch: " + role)
        result[role] = path
    return result

rule c3_stan_cv_score:
    input:
        features=OUT + "/{sample}/c1/mva_feature_matrix.tsv.gz",
        bundle=_CV_BUNDLE,
        bundle_sha=_CV_BUNDLE + ".sha256",
        scorer=os.path.join(_CV_SOURCE, "b14_stan_cv_pipeline.R"),
        helper=os.path.join(_CV_SOURCE, "b14_stan_cv_helpers.R"),
        runner=_CV_RUNNER,
        stan=os.path.join(_CV_STAN_DIR, "m3b.stan"),
        config=CFG
    output:
        predictions=OUT + "/{sample}/c3_stan_cv/stan_cv_predictions.rds",
        checksum=OUT + "/{sample}/c3_stan_cv/stan_cv_predictions.rds.sha256",
        summary=OUT + "/{sample}/c3_stan_cv/stan_cv_predictions_summary.tsv"
    params:
        protocol=config["frozen_stan"]["candidate_protocol"],
        stan_dir=_CV_STAN_DIR,
        rscript=config["frozen_stan"]["rscript"]
    threads: 1
    log:
        OUT + "/{sample}/logs/c3_stan_cv.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"

        env -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE \
            -u LD_LIBRARY_PATH -u LIBRARY_PATH \
            {params.rscript:q} --vanilla {input.scorer:q} \
              --mode score \
              --bundle {input.bundle:q} \
              --features {input.features:q} \
              --candidate-protocol {params.protocol:q} \
              --runner {input.runner:q} \
              --helpers {input.helper:q} \
              --stan-dir {params.stan_dir:q} \
              --confirm-raw-features yes \
              --out {output.predictions:q} \
              > {log:q} 2>&1
        """

rule c4_fresh_stan_union:
    input:
        a10=OUT + "/{sample}/a10/ranked_candidates.tsv.gz",
        a2=OUT + "/{sample}/a2/shortlist.tsv.gz",
        features=OUT + "/{sample}/c1/mva_feature_matrix.tsv.gz",
        predictions=OUT + "/{sample}/c3_stan_cv/stan_cv_predictions.rds",
        prediction_sha=OUT + "/{sample}/c3_stan_cv/stan_cv_predictions.rds.sha256",
        history=lambda wc: _fresh_union_history(wc)["manifest"],
        legacy_c4=lambda wc: _fresh_union_history(wc)["legacy_c4"],
        legacy_a10=lambda wc: _fresh_union_history(wc)["legacy_a10"],
        code=_UNION_CODE,
        helpers=_UNION_LIBS,
        config=CFG
    output:
        registry=_UNION_OUT + "/candidate_registry.tsv",
        variants=_UNION_OUT + "/shortlisted_variant_evidence.tsv",
        weights=_UNION_OUT + "/current_Stan_hypothesis_weights.tsv",
        hypotheses=_UNION_OUT + "/current_Stan_final_hypothesis_set.tsv",
        scenarios=_UNION_OUT + "/current_Stan_scenarios.tsv",
        protected=_UNION_OUT + "/protected_unscored_candidates.tsv",
        historical=_UNION_OUT + "/historical_C4_preserved.tsv",
        unresolved=_UNION_OUT + "/historical_unresolved_rows.tsv",
        contract=_UNION_OUT + "/HANDOFF_CONTRACT.txt",
        manifest=_UNION_OUT + "/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/c4_stan_union"
    log:
        OUT + "/{sample}/logs/c4_stan_union.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        stage=$(mktemp -d {params.out:q}/.stage_XXXXXXXX)
        trap 'rm -rf "$stage"' EXIT
        Rscript --vanilla {input.code:q} \
          --config {input.config:q} \
          --sample {wildcards.sample:q} \
          --a10 {input.a10:q} \
          --a2 {input.a2:q} \
          --legacy_c4 {input.legacy_c4:q} \
          --legacy_a10 {input.legacy_a10:q} \
          --features {input.features:q} \
          --predictions {input.predictions:q} \
          --current_build GRCh38 \
          --legacy_build GRCh38 \
          --outdir "$stage/run" > {log:q} 2>&1
        test -s "$stage/run/manifest.json"
        for artifact in "$stage/run/"*; do
          if [ "$(basename "$artifact")" != "manifest.json" ]; then
            mv -- "$artifact" {params.out:q}/
          fi
        done
        mv -- "$stage/run/manifest.json" {params.out:q}/
        """

rule review_a11_evidence:
    input:
        a10=OUT + "/{sample}/a10/ranked_candidates.tsv.gz",
        code=os.path.join(_REVIEW_WF, "scripts", "make_a11_input.R")
    output:
        evidence=_A11_EVIDENCE
    params:
        directory=lambda wc: OUT + "/" + wc.sample + "/a11_review"
    log:
        OUT + "/{sample}/logs/review_a11_evidence.log"
    shell:
        r"""
        mkdir -p {params.directory:q} "$(dirname {log:q})"
        stage=$(mktemp {params.directory:q}/.a10_input.XXXXXXXX)
        trap 'rm -f "$stage"' EXIT
        gzip -cd {input.a10:q} > "$stage"
        Rscript --vanilla {input.code:q} "$stage" \
          {output.evidence:q} > {log:q} 2>&1
        """

rule review_a11_calculation:
    input:
        evidence=_A11_EVIDENCE,
        code=os.path.join(
            _REVIEW_WF, "layer_a", "scripts",
            "a11_dependence_sensitivity.R"
        ),
        runtime=os.path.join(_REVIEW_WF, "lib", "a11_runtime.R")
    output:
        consensus=_A11_CALC + "/bayesian_dependence_sensitivity_ensemble_consensus.csv",
        coverage=_A11_CALC + "/candidate_coverage.csv",
        copula=_A11_CALC + "/copula_dependence_summary.csv",
        provenance=_A11_CALC + "/run_provenance.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/a11_review/calculation",
        seed=int(config.get("a11", {}).get("seed", 20260912)),
        draws=int(config.get("a11", {}).get("draws", 10000))
    log:
        OUT + "/{sample}/logs/review_a11_calculation.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        stage=$(mktemp -d {params.out:q}/.stage_XXXXXXXX)
        trap 'rm -rf "$stage"' EXIT
        Rscript --vanilla {input.code:q} \
          --input {input.evidence:q} \
          --outdir "$stage/run" \
          --seed {params.seed} \
          --draws {params.draws} > {log:q} 2>&1
        test -s "$stage/run/run_provenance.json"
        for artifact in "$stage/run/"*; do
          if [ "$(basename "$artifact")" != "run_provenance.json" ]; then
            mv -- "$artifact" {params.out:q}/
          fi
        done
        mv -- "$stage/run/run_provenance.json" {params.out:q}/
        """


# D0 consumes the NEW union and NEW model views.
# Historical support is resolved only by _fresh_union_history above.
def _t2_union_dir(wc):
    return os.path.join(OUT, wc.sample, "c4_stan_union")

def _t2_model_sources(wc):
    return {
        "a10": os.path.join(OUT, wc.sample, "a10", "ranked_candidates.tsv.gz"),
        "dependence": os.path.join(
            OUT, wc.sample, "a11_review", "calculation",
            "bayesian_dependence_sensitivity_ensemble_consensus.csv"
        ),
    }

def _t2_union_inputs(wc):
    directory = _t2_union_dir(wc)
    return [
        os.path.join(directory, name)
        for name in (
            "manifest.json", "current_Stan_scenarios.tsv",
            "shortlisted_variant_evidence.tsv",
            "candidate_registry.tsv", "protected_unscored_candidates.tsv"
        )
    ] + [
        os.path.join(RES, "manifests", "resource_manifest.json"),
        os.path.join(RES, "manifests", "checksums.sha256"),
        os.path.join(
            OUT, wc.sample, "a11_review", "calculation",
            "run_provenance.json"
        ),
        os.path.join(
            OUT, wc.sample, "a11_review", "calculation",
            "candidate_coverage.csv"
        ),
    ]


# Include explicit symbol-key coverage and retained evidence in D0 dependencies.
_t2_union_inputs_before_gene_key_audit = _t2_union_inputs

def _t2_union_inputs(wc):
    return _t2_union_inputs_before_gene_key_audit(wc) + [
        os.path.join(OUT, wc.sample, "d7", "gene_key_coverage.json"),
        os.path.join(
            OUT, wc.sample, "d7",
            "symbol_unresolved_variant_evidence.tsv.gz"
        ),
    ]
