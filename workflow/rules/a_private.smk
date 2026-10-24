# a_private.smk -- Layer A (private). Each rule: explicit inputs/outputs/logs; no scientific logic here.
PV = config["private"]
A = OUT + "/{sample}"

# Candidate scope is explicit; no AF or technical-QC threshold is changed here.
import glob
import json
A1_OPTIONS = config.get("a1", {})
if not isinstance(A1_OPTIONS, dict):
    raise ValueError("a1 must be a mapping")
A1_POLICY = A1_OPTIONS.get("candidate_policy", "legacy_any_csq_token_position_v1")
A1_CHUNK = A1_OPTIONS.get("chunk_records", 1000)
A1_PROGRESS = A1_OPTIONS.get("progress_every", 100000)
if A1_POLICY != "legacy_any_csq_token_position_v1":
    raise ValueError("Unsupported a1.candidate_policy")
if isinstance(A1_CHUNK, bool) or not isinstance(A1_CHUNK, int) or not 1 <= A1_CHUNK <= 20000:
    raise ValueError("a1.chunk_records must be an integer from 1 to 20000")
if isinstance(A1_PROGRESS, bool) or not isinstance(A1_PROGRESS, int) or A1_PROGRESS < 1:
    raise ValueError("a1.progress_every must be a positive integer")

def _candidate_a1_input_vcf(wc):
    # Evaluate before job execution, so Snakemake cannot remove an old full A1.
    od = os.path.join(OUT, wc.sample, "a1")
    products = [os.path.join(od, name) for name in
                ("evidence_matrix.tsv.gz", "qc_summary.tsv", "manifest.json")]
    marker = os.path.join(od, "candidate_scope.json")
    if any(os.path.exists(p) for p in products):
        if not os.path.isfile(marker):
            raise ValueError("Existing A1 has no candidate-scope marker. Preserve it and use a fresh run.output_root: " + od)
        with open(marker) as handle:
            scope = json.load(handle)
        if scope.get("status") != "complete" or scope.get("policy") != A1_POLICY:
            raise ValueError("Existing A1 scope is incompatible; use a fresh run.output_root: " + od)
    return os.path.join(OUT, wc.sample, "a0_candidates", "candidates.vcf.gz")


if config["run"]["mode"] == "existing_vep":
    rule a0_existing_vep:
        input: vcf=PV["existing_biallelic_vep_vcf"], tbi=PV["existing_biallelic_vep_vcf_index"], metadata=PV["vep_run_metadata"]
        output: vcf=A + "/a0/annotated.biallelic.vcf.gz", tbi=A + "/a0/annotated.biallelic.vcf.gz.tbi", manifest=A + "/a0/a0_manifest.json", schema=A + "/a0/csq_schema.txt"
        log: A + "/logs/a0_existing_vep.log"
        shell: "bash {SCRIPTS}/a0_existing_vep.sh --config {CFG} --sample {wildcards.sample} --vcf {input.vcf} --tbi {input.tbi} --metadata {input.metadata} --outdir {OUT}/{wildcards.sample}/a0 > {log} 2>&1"
else:
    rule a0_normalize_and_vep:
        input: vcf=PV["raw_vcf"], tbi=PV["raw_vcf_index"], fasta=PV["reference_fasta"] or "MISSING_reference_fasta", fai=PV["reference_fasta_index"] or "MISSING_reference_fasta_index"
        output: normalized=A + "/a0/normalized.biallelic.vcf.gz", vcf=A + "/a0/annotated.biallelic.vcf.gz", tbi=A + "/a0/annotated.biallelic.vcf.gz.tbi", manifest=A + "/a0/a0_manifest.json", schema=A + "/a0/csq_schema.txt"
        log: A + "/logs/a0_rerun_vep.log"
        threads: 4
        shell: "bash {SCRIPTS}/a0_normalize_and_vep.sh --config {CFG} --sample {wildcards.sample} --vcf {input.vcf} --tbi {input.tbi} --fasta {input.fasta} --fai {input.fai} --cache_dir {PV[vep_cache_dir]} --outdir {OUT}/{wildcards.sample}/a0 --threads {threads} > {log} 2>&1"

rule a0_candidate_subset:
    input:
        vcf=A + "/a0/annotated.biallelic.vcf.gz",
        provenance=A + "/a0/a0_manifest.json",
        code=SCRIPTS + "/a0_prepare_candidates.sh"
    output:
        vcf=A + "/a0_candidates/candidates.vcf.gz",
        tbi=A + "/a0_candidates/candidates.vcf.gz.tbi",
        manifest=A + "/a0_candidates/candidate_manifest.json"
    params:
        policy=A1_POLICY,
        progress=A1_PROGRESS,
        outdir=lambda wc: os.path.join(OUT, wc.sample, "a0_candidates")
    log: A + "/logs/a0_candidates.log"
    shell:
        "bash {input.code:q} --input {input.vcf:q} --outdir {params.outdir:q} "
        "--policy {params.policy:q} --progress-every {params.progress} > {log:q} 2>&1"

rule a1_streaming_selftest:
    input:
        vcf=A + "/a0_candidates/candidates.vcf.gz",
        scope=A + "/a0_candidates/candidate_manifest.json",
        code=SCRIPTS + "/a1_build_evidence_matrix.R",
        config=CFG,
        helpers=sorted(glob.glob(os.path.join(SCRIPTS, "..", "lib", "*.R")))
    output:
        parity=A + "/a1_selftest/parity_report.json"
    params:
        parent=lambda wc: os.path.join(OUT, wc.sample, "a1_selftest")
    log: A + "/logs/a1_selftest.log"
    shell:
        r"""
        test_parent=$(mktemp -d {params.parent:q}/attempt_XXXXXXXX)
        test_dir="$test_parent/run"
        {RSCRIPT} {input.code:q} --config {input.config:q} \
            --sample {wildcards.sample:q} --vcf {input.vcf:q} \
            --outdir "$test_dir" --mode selftest > {log:q} 2>&1
        test -s "$test_dir/parity_report.json"
        cp -- "$test_dir/parity_report.json" {output.parity:q}
        """

rule a1_evidence_matrix:
    input:
        vcf=_candidate_a1_input_vcf,
        scope=A + "/a0_candidates/candidate_manifest.json",
        parity=A + "/a1_selftest/parity_report.json",
        code=SCRIPTS + "/a1_build_evidence_matrix.R",
        config=CFG,
        helpers=sorted(glob.glob(os.path.join(SCRIPTS, "..", "lib", "*.R")))
    output:
        m=A + "/a1/evidence_matrix.tsv.gz",
        q=A + "/a1/qc_summary.tsv",
        manifest=A + "/a1/manifest.json",
        scope=A + "/a1/candidate_scope.json"
    params:
        chunk=A1_CHUNK,
        outdir=lambda wc: os.path.join(OUT, wc.sample, "a1")
    log: A + "/logs/a1.log"
    shell:
        r"""
        {RSCRIPT} {input.code:q} --config {input.config:q} \
            --sample {wildcards.sample:q} --vcf {input.vcf:q} \
            --outdir {params.outdir:q} --mode build \
            --parity {input.parity:q} --chunk_records {params.chunk} > {log:q} 2>&1
        cp -- {input.scope:q} {output.scope:q}
        """

rule a2_shortlist:
    input: m=A + "/a1/evidence_matrix.tsv.gz"
    output: sl=A + "/a2/shortlist.tsv.gz", ex=A + "/a2/exclusion_summary.tsv", fa=A + "/a2/filter_audit.tsv.gz", pr=A + "/a2/pair_candidates.tsv", manifest=A + "/a2/manifest.json"
    log: A + "/logs/a2.log"
    shell: "{RSCRIPT} {SCRIPTS}/a2_build_candidate_shortlist.R --config {CFG} --sample {wildcards.sample} --matrix {input.m} --outdir {OUT}/{wildcards.sample}/a2 > {log} 2>&1"

rule a3_hgvs:
    input: vcf=A + "/a0_candidates/candidates.vcf.gz", sl=A + "/a2/shortlist.tsv.gz"
    output: sl=A + "/a3/shortlist_hgvs.tsv.gz", audit=A + "/a3/hgvs_merge_audit.tsv", manifest=A + "/a3/manifest.json"
    log: A + "/logs/a3.log"
    shell: "python3 {SCRIPTS}/a3_merge_hgvs.py --vcf {input.vcf} --shortlist {input.sl} --out {output.sl} --audit {output.audit} > {log} 2>&1 && "
           "{RSCRIPT} {SCRIPTS}/write_manifest.R --config {CFG} --node a3 --outdir {OUT}/{wildcards.sample}/a3 --inputs {input.vcf},{input.sl} --outputs {output.sl},{output.audit} >> {log} 2>&1"

rule a4_hpo:
    input: sl=A + "/a3/shortlist_hgvs.tsv.gz"
    output: hpo=A + "/a4/hpo_gene_scores.tsv", manifest=A + "/a4/manifest.json"
    log: A + "/logs/a4.log"
    shell: "{RSCRIPT} {SCRIPTS}/a4_hpo_multimetric.R --config {CFG} --sample {wildcards.sample} --shortlist {input.sl} --outdir {OUT}/{wildcards.sample}/a4 > {log} 2>&1"

rule a5_build_network:
    input: sl=A + "/a3/shortlist_hgvs.tsv.gz"
    output: g=A + "/a5/graphs.rds", manifest=A + "/a5/manifest_graphs.json"
    log: A + "/logs/a5_graphs.log"
    shell: "{RSCRIPT} {SCRIPTS}/a5_build_network.R --config {CFG} --sample {wildcards.sample} --outdir {OUT}/{wildcards.sample}/a5 > {log} 2>&1"

rule a5_rwr:
    input: sl=A + "/a3/shortlist_hgvs.tsv.gz", hpo=A + "/a4/hpo_gene_scores.tsv", g=A + "/a5/graphs.rds"
    output: n=A + "/a5/network_gene_scores.tsv", grid=A + "/a5/rwr_grid.tsv", manifest=A + "/a5/manifest.json"
    log: A + "/logs/a5_rwr.log"
    shell: "{RSCRIPT} {SCRIPTS}/a5_rwr.R --config {CFG} --sample {wildcards.sample} --shortlist {input.sl} --hpo {input.hpo} --graphs {input.g} --outdir {OUT}/{wildcards.sample}/a5 > {log} 2>&1"

rule a6_multiview:
    input: sl=A + "/a3/shortlist_hgvs.tsv.gz", hpo=A + "/a4/hpo_gene_scores.tsv", n=A + "/a5/network_gene_scores.tsv"
    output: a=A + "/a6/view_agreement.tsv", c=A + "/a6/coverage_counts.tsv", manifest=A + "/a6/manifest.json"
    log: A + "/logs/a6.log"
    shell: "{RSCRIPT} {SCRIPTS}/a6_multiview.R --config {CFG} --sample {wildcards.sample} --shortlist {input.sl} --hpo {input.hpo} --network {input.n} --outdir {OUT}/{wildcards.sample}/a6 > {log} 2>&1"

rule a7_integrate:
    input: sl=A + "/a3/shortlist_hgvs.tsv.gz", hpo=A + "/a4/hpo_gene_scores.tsv", n=A + "/a5/network_gene_scores.tsv", v=A + "/a6/view_agreement.tsv"
    output: i=A + "/a7/integrated_scores.tsv.gz", manifest=A + "/a7/manifest.json"
    log: A + "/logs/a7.log"
    shell: "{RSCRIPT} {SCRIPTS}/a7_integrate.R --config {CFG} --sample {wildcards.sample} --shortlist {input.sl} --hpo {input.hpo} --network {input.n} --outdir {OUT}/{wildcards.sample}/a7 > {log} 2>&1"

rule a7b_dependence:
    input: i=A + "/a7/integrated_scores.tsv.gz"
    output: d=A + "/a7/dependence.tsv", manifest=A + "/a7/manifest_a7b.json"
    log: A + "/logs/a7b.log"
    shell: "{RSCRIPT} {SCRIPTS}/a7b_dependence.R --config {CFG} --sample {wildcards.sample} --integrated {input.i} --outdir {OUT}/{wildcards.sample}/a7 > {log} 2>&1"

rule a7c_shapley:
    input: i=A + "/a7/integrated_scores.tsv.gz"
    output: s=A + "/a7/shapley_top.tsv", manifest=A + "/a7/manifest_a7c.json"
    log: A + "/logs/a7c.log"
    shell: "{RSCRIPT} {SCRIPTS}/a7c_shapley.R --config {CFG} --sample {wildcards.sample} --integrated {input.i} --outdir {OUT}/{wildcards.sample}/a7 > {log} 2>&1"

rule a8_mosaicism:
    input: sl=A + "/a2/shortlist.tsv.gz"
    output: m=A + "/a8/mosaicism.tsv", manifest=A + "/a8/manifest.json"
    log: A + "/logs/a8.log"
    shell: "{RSCRIPT} {SCRIPTS}/a8_mosaicism.R --config {CFG} --sample {wildcards.sample} --shortlist {input.sl} --outdir {OUT}/{wildcards.sample}/a8 > {log} 2>&1"

rule a9_stability:
    input: i=A + "/a7/integrated_scores.tsv.gz", m=A + "/a8/mosaicism.tsv"          # A7 + A8 must both exist before A9
    output: f=A + "/a9/top_m_frequency.tsv", l=A + "/a9/leave_one_channel_out.tsv.gz", s=A + "/a9/stability_summary.tsv", manifest=A + "/a9/manifest.json"
    log: A + "/logs/a9.log"
    shell: "{RSCRIPT} {SCRIPTS}/a9_stability_influence.R --config {CFG} --sample {wildcards.sample} --integrated {input.i} --outdir {OUT}/{wildcards.sample}/a9 > {log} 2>&1"


# Preserve symbol-unresolved variants outside the symbol-keyed gene ranking.
rule a10_gene_key_partition:
    input:
        integrated=A + "/a7/integrated_scores.tsv.gz",
        vcf=A + "/a0_candidates/candidates.vcf.gz",
        code=SCRIPTS + "/a10_partition_gene_keys.py"
    output:
        named=A + "/a10_gene_keys/named_integrated_scores.tsv.gz",
        unresolved=A + "/d7/symbol_unresolved_variant_evidence.tsv.gz",
        audit=A + "/d7/gene_key_coverage.json"
    log:
        A + "/logs/a10_gene_key_partition.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        python3 {input.code:q} \
          --input {input.integrated:q} \
          --vcf {input.vcf:q} \
          --named {output.named:q} \
          --unresolved {output.unresolved:q} \
          --audit {output.audit:q} > {log:q} 2>&1
        """

rule a10_aggregate:
    input: i=A + "/a10_gene_keys/named_integrated_scores.tsv.gz", gene_key_audit=A + "/d7/gene_key_coverage.json", symbol_unresolved=A + "/d7/symbol_unresolved_variant_evidence.tsv.gz", m=A + "/a8/mosaicism.tsv", f=A + "/a9/top_m_frequency.tsv", v=A + "/a6/view_agreement.tsv", d=A + "/a7/dependence.tsv", s=A + "/a7/shapley_top.tsv"
    output: r=A + "/a10/ranked_candidates.tsv.gz", manifest=A + "/a10/manifest.json"
    log: A + "/logs/a10.log"
    shell: "{RSCRIPT} {SCRIPTS}/a10_aggregate.R --config {CFG} --sample {wildcards.sample} --integrated {input.i} --mosaicism {input.m} --stability {input.f} --outdir {OUT}/{wildcards.sample}/a10 > {log} 2>&1"
