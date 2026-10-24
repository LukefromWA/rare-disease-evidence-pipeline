
def _t2_hierarchy_inputs(wc):
    directory = config.get("track2", {}).get("reactome_hierarchy_dir")
    if not isinstance(directory, str) or not directory:
        raise ValueError("Configure track2.reactome_hierarchy_dir")
    return [
        directory.rstrip("/") + "/" + name
        for name in (
            "checksums.sha256",
            "manifest.json",
            "human_pathway_hierarchy.tsv",
            "ReactomePathways.txt",
            "ReactomePathwaysRelation.txt",
        )
    ]


def _t2_model_sources(wc):
    value = config.get("track2", {}).get(
        "model_view_sources", {}
    ).get(wc.sample)
    if not isinstance(value, dict) or set(value) != {"a10", "dependence"}:
        raise ValueError(
            "Configure track2.model_view_sources for this sample: "
            "a10 and dependence"
        )
    return value

def _t2_model_inputs(wc):
    value = _t2_model_sources(wc)
    return [value["a10"], value["dependence"]]

def _t2_model_args(wc):
    value = _t2_model_sources(wc)
    quote = __import__("shlex").quote
    return (
        "--a10_views " + quote(value["a10"]) +
        " --dependence_views " + quote(value["dependence"])
    )


def _t2_union_dir(wc):
    return config.get("track2", {}).get(
        "frozen_union_dirs", {}
    ).get(wc.sample, "")

def _t2_union_inputs(wc):
    directory = _t2_union_dir(wc)
    if not directory:
        return []
    return [
        os.path.join(directory, name)
        for name in (
            "manifest.json",
            "current_Stan_scenarios.tsv",
            "shortlisted_variant_evidence.tsv",
            "candidate_registry.tsv",
            "protected_unscored_candidates.tsv",
        )
    ]

def _t2_union_args(wc):
    directory = _t2_union_dir(wc)
    return (
        "--frozen_union_dir " +
        __import__("shlex").quote(directory)
    ) if directory else ""

# Targeted labeling integration: no bulk label database.
from shlex import quote as _t2_quote

def _t2_review_inputs(wc):
    path = config.get("track2", {}).get(
        "label_review_overrides", {}
    ).get(wc.sample)
    return [path] if path else []

# d_track2.smk -- Layer D. Every rule depends on the frozen resource bundle (.ready) and, through D0, on the C4 handoff.
READY = RES + "/.ready"; D = OUT + "/{sample}"
rule d0_scenarios:
    input:
        c4=lambda wc: (
            os.path.join(
                _t2_union_dir(wc),
                "current_Stan_final_hypothesis_set.tsv"
            ) if _t2_union_dir(wc) else
            (OUT + "/{sample}/c4/final_hypothesis_set.tsv.gz").format(
                sample=wc.sample
            )
        ),
        frozen=lambda wc: _t2_union_inputs(wc),
        model_sources=_t2_model_inputs,
        ready=READY,
        code=SCRIPTS + "/d0_scenarios.R",
        cfg=CFG,
        private_cfg=(
            [os.path.join(os.path.dirname(CFG), "config.private.yaml")]
            if os.path.isfile(os.path.join(
                os.path.dirname(CFG), "config.private.yaml"
            )) else []
        )
    output:
        s=D + "/d0/scenarios.tsv.gz",
        a=D + "/d0/scenario_audit.tsv",
        manifest=D + "/d0/manifest.json",
        views=D + "/d0/candidate_model_views.tsv.gz",
        view_coverage=D + "/d0/model_view_coverage.tsv"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d0",
        frozen_args=_t2_union_args,
        model_args=_t2_model_args
    log: D + "/logs/d0.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --config {CFG:q} "
        "--sample {wildcards.sample:q} --c4 {input.c4:q} "
        "--outdir {params.out:q} {params.frozen_args} {params.model_args} > {log:q} 2>&1"

rule d1_mechanism:
    input: s=D + "/d0/scenarios.tsv.gz"
    output: m=D + "/d1/mechanism.tsv.gz", manifest=D + "/d1/manifest.json"
    log: D + "/logs/d1.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d1_mechanism.R --config {CFG:q} --sample {wildcards.sample} --scenarios {input.s:q} --outdir {OUT:q}/{wildcards.sample}/d1 > {log:q} 2>&1"
rule d2_targets:
    input:
        m=D + "/d1/mechanism.tsv.gz",
        views=D + "/d0/candidate_model_views.tsv.gz",
        ready=READY,
        drug_map=RES + "/derived/drug_target_action_map.tsv.gz",
        pathways=RES + "/derived/reactome_gene_sets.tsv.gz",
        hierarchy=_t2_hierarchy_inputs,
        code=SCRIPTS + "/d2_targets.R"
    output:
        t=D + "/d2/target_universe.tsv.gz",
        audit=D + "/d2/hypothesis_discovery_audit.tsv.gz",
        manifest=D + "/d2/manifest.json",
        pathway_pairs=D + "/d2/pathway_discovery.tsv.gz",
        pathway_witnesses=D + "/d2/pathway_discovery_witnesses.tsv.gz",
        pathway_audit=D + "/d2/pathway_discovery_audit.tsv.gz"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d2",
        hierarchy_dir=lambda wc: config["track2"]["reactome_hierarchy_dir"]
    log: D + "/logs/d2.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {RSCRIPT} {input.code:q} \
          --config {CFG:q} \
          --sample {wildcards.sample:q} \
          --mechanism {input.m:q} \
          --model_views {input.views:q} \
          --hierarchy_dir {params.hierarchy_dir:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule d3_1:
    input: gates=SCRIPTS + "/../lib/gates.R", t=D + "/d2/target_universe.tsv.gz"
    output: o=D + "/d3/method1_target_query.tsv.gz", manifest=D + "/d3/manifest_d3_1.json"
    log: D + "/logs/d3_1.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d3_1_target_query.R --config {CFG:q} --sample {wildcards.sample} --targets {input.t:q} --outdir {OUT:q}/{wildcards.sample}/d3 > {log:q} 2>&1"
rule d3_2:
    output: D + "/d3/d3_2_signature_reversal_status.tsv"
    log: D + "/logs/d3_2.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d3_2_signature_reversal.R --config {CFG:q} --sample {wildcards.sample} --outdir {OUT:q}/{wildcards.sample}/d3 > {log:q} 2>&1"
rule d3_3:
    input: t=D + "/d2/target_universe.tsv.gz", ready=READY
    output: o=D + "/d3/method3_pathway.tsv.gz", manifest=D + "/d3/manifest_d3_3.json"
    log: D + "/logs/d3_3.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d3_3_pathway.R --config {CFG:q} --sample {wildcards.sample} --targets {input.t:q} --outdir {OUT:q}/{wildcards.sample}/d3 > {log:q} 2>&1"
rule d3_4:
    output: D + "/d3/d3_4_phenotype_similarity_status.tsv"
    log: D + "/logs/d3_4.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d3_4_phenotype_similarity.R --config {CFG:q} --sample {wildcards.sample} --outdir {OUT:q}/{wildcards.sample}/d3 > {log:q} 2>&1"
rule d3_5:
    input:
        scenarios=D + "/d0/scenarios.tsv.gz",
        code=SCRIPTS + "/d3_5_literature.R",
        cfg=CFG,
        registry=lambda wc: config["track2"]["mechanism_search"]["compound_seed_registry"],
        registry_manifest=lambda wc: config["track2"]["mechanism_search"]["compound_seed_registry"] + ".manifest.json"
    output:
        status=D + "/d3/d3_5_literature_status.tsv",
        queries=D + "/d3/mechanism_search_queries.tsv",
        publications=D + "/d3/mechanism_publications.tsv.gz",
        links=D + "/d3/hypothesis_mechanism_search_links.tsv.gz",
        audit=D + "/d3/mechanism_search_hypothesis_audit.tsv.gz",
        worksheet=D + "/d3/mechanism_compound_review_worksheet.tsv",
        manifest=D + "/d3/manifest_d3_5_literature.json",
        compound_evidence=D + "/d3/mechanism_compound_evidence.tsv.gz",
        chemical_mentions=D + "/d3/mechanism_chemical_index_mentions.tsv.gz",
        extraction_audit=D + "/d3/mechanism_compound_extraction_audit.tsv"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3",
        search_settings=lambda wc: repr(
            config.get("track2", {}).get("mechanism_search", {})
        )
    log: D + "/logs/d3_5.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {RSCRIPT} {input.code:q} \
          --config {CFG:q} \
          --sample {wildcards.sample:q} \
          --scenarios {input.scenarios:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule d3_6:
    input:
        code=SCRIPTS + "/d3_6_signed_propagation.R",
        gates=SCRIPTS + "/../lib/gates.R",
        t=D + "/d2/target_universe.tsv.gz",
        ready=READY
    output:
        o=D + "/d3/method6_propagation.tsv.gz",
        s=D + "/d3/d3_6_signed_propagation_status.tsv",
        manifest=D + "/d3/manifest_d3_6_signed_propagation.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3"
    log: D + "/logs/d3_6.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {RSCRIPT} {input.code:q} \
          --config {CFG:q} \
          --sample {wildcards.sample:q} \
          --targets {input.t:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """
# END disabled_rescue_status_v1

rule d3_7:
    output: D + "/d3/d3_7_controllability_status.tsv"
    log: D + "/logs/d3_7.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d3_7_controllability.R --config {CFG:q} --sample {wildcards.sample} --outdir {OUT:q}/{wildcards.sample}/d3 > {log:q} 2>&1"
rule d3_8:
    input:
        s7=D + "/d3/d3_7_controllability_status.tsv",
        s6=D + "/d3/d3_6_signed_propagation_status.tsv",
        code=SCRIPTS + "/d3_8_rescue_optimization.R"
    output: D + "/d3/d3_8_rescue_optimization_status.tsv"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3"
    log: D + "/logs/d3_8.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {RSCRIPT} {input.code:q} \
          --config {CFG:q} \
          --sample {wildcards.sample:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule d4_integrate:
    input:
        t=D + "/d2/target_universe.tsv.gz", m1=D + "/d3/method1_target_query.tsv.gz", m3=D + "/d3/method3_pathway.tsv.gz", m6=D + "/d3/method6_propagation.tsv.gz",
        s2=D + "/d3/d3_2_signature_reversal_status.tsv", s4=D + "/d3/d3_4_phenotype_similarity_status.tsv", s5=D + "/d3/d3_5_literature_status.tsv", s8=D + "/d3/d3_8_rescue_optimization_status.tsv"   # all D3 outputs or explicit statuses
    output: e=D + "/d4/evidence.tsv.gz", manifest=D + "/d4/manifest.json"
    log: D + "/logs/d4.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d4_integrate.R --config {CFG:q} --sample {wildcards.sample} --targets {input.t:q} --method1 {input.m1:q} --method3 {input.m3:q} --method6 {input.m6:q} --outdir {OUT:q}/{wildcards.sample}/d4 > {log:q} 2>&1"
rule d5_safety_gate:
    input:
        e=D + "/d4/evidence.tsv.gz",
        ready=READY,
        safety=D + "/label_evidence/pediatric_safety_map.tsv.gz",
        lm=D + "/label_evidence/run_manifest.json",
        ls=D + "/label_evidence/checksums.sha256",
        targets=D + "/d2/target_universe.tsv.gz",
        code=SCRIPTS + "/d5_safety_gate.R",
        gates=SCRIPTS + "/../lib/gates.R"
    output:
        g=D + "/d5/gated.tsv.gz",
        x=D + "/d5/excluded.tsv",
        held=D + "/d5/held_or_excluded.tsv",
        manifest=D + "/d5/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d5"
    log: D + "/logs/d5.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --config {CFG:q} "
        "--sample {wildcards.sample:q} --evidence {input.e:q} "
        "--safety {input.safety:q} --targets {input.targets:q} "
        "--outdir {params.out:q} > {log:q} 2>&1"

rule d6_sensitivity:
    input:
        code=SCRIPTS + "/d6_sensitivity.R",
        gates=SCRIPTS + "/../lib/gates.R",
        g=D + "/d5/gated.tsv.gz",
        scenarios=D + "/d0/scenarios.tsv.gz"
    output:
        s=D + "/d6/sensitivity.tsv",
        coverage=D + "/d6/scenario_coverage.tsv",
        status=D + "/d6/status.tsv",
        manifest=D + "/d6/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d6"
    log: D + "/logs/d6.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        {RSCRIPT} {input.code:q} \
          --config {CFG:q} \
          --sample {wildcards.sample:q} \
          --gated {input.g:q} \
          --scenarios {input.scenarios:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule d7_report:
    input: code=SCRIPTS + "/d7_report.R", g=D + "/d5/gated.tsv.gz", s=D + "/d6/sensitivity.tsv"
    output: r=D + "/d7/final_report.tsv", manifest=D + "/d7/manifest.json"
    log: D + "/logs/d7.log"
    shell: "{RSCRIPT} {SCRIPTS:q}/d7_report.R --config {CFG:q} --sample {wildcards.sample} --gated {input.g:q} --sensitivity {input.s:q} --outdir {OUT:q}/{wildcards.sample}/d7 > {log:q} 2>&1"


rule d2_label_evidence:
    input:
        drugs=D + "/d2/target_universe.tsv.gz",
        code=SCRIPTS + "/fetch_label_evidence.py",
        review=_t2_review_inputs
    output:
        safety=D + "/label_evidence/pediatric_safety_map.tsv.gz",
        status=D + "/label_evidence/drug_lookup_status.tsv",
        evidence=D + "/label_evidence/label_evidence.tsv.gz",
        sections=D + "/label_evidence/label_sections.tsv.gz",
        worksheet=D + "/label_evidence/review_worksheet.tsv",
        manifest=D + "/label_evidence/run_manifest.json",
        sums=D + "/label_evidence/checksums.sha256"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/label_evidence",
        cache=RES + "/cache/openfda",
        review=lambda wc: (
            "--review-overrides " +
            _t2_quote(_t2_review_inputs(wc)[0])
        ) if _t2_review_inputs(wc) else "",
        maxd=int(config.get("track2", {}).get("label_max_drugs", 200))
    log: D + "/logs/d2_label_evidence.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q}) {params.out:q} {params.cache:q}
        stage=$(mktemp -d {params.out:q}/.staging.XXXXXXXX)
        trap 'rm -rf "$stage"' EXIT

        python3 {input.code:q} \
          --drugs {input.drugs:q} \
          --out-dir "$stage" \
          --cache-dir {params.cache:q} \
          --max-drugs {params.maxd} \
          {params.review} > {log:q} 2>&1

        (cd "$stage" && sha256sum -c checksums.sha256) >> {log:q} 2>&1

        for name in drug_lookup_status.tsv label_evidence.tsv.gz label_sections.tsv.gz review_worksheet.tsv pediatric_safety_map.tsv.gz checksums.sha256; do
          mv -f "$stage/$name" {params.out:q}/
        done
        mv -f "$stage/run_manifest.json" {params.out:q}/
        """


# BEGIN d3_codon_context_v1
def _t2_codon_variant_evidence(wc):
    return OUT + "/" + wc.sample + "/c4_stan_union/shortlisted_variant_evidence.tsv"

def _t2_codon_a10(wc):
    return OUT + "/" + wc.sample + "/a10/ranked_candidates.tsv.gz"

def _t2_codon_gene_args(wc):
    genes = config["track2"]["mechanism_search"]["context_genes"][wc.sample]
    if not isinstance(genes, list) or not genes:
        raise ValueError("Configure a nonempty context_genes list")
    quote = __import__("shlex").quote
    return " ".join("--gene " + quote(str(gene)) for gene in genes)

rule d3_codon_context:
    input:
        evidence=_t2_codon_variant_evidence,
        a10=_t2_codon_a10,
        code=SCRIPTS + "/d3_codon_context.py"
    output:
        table=D + "/d3/codon_context/codon_context.tsv",
        fasta=D + "/d3/codon_context/flanking_context.fasta",
        manifest=D + "/d3/codon_context/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3/codon_context",
        cache=RES + "/supplemental/ensembl_cds",
        genes=_t2_codon_gene_args,
        server=lambda wc: config["track2"]["mechanism_search"]["context_rest_server"],
        flank=lambda wc: int(config["track2"]["mechanism_search"]["context_flank"]),
        live=lambda wc: (
            "--allow-live"
            if config["track2"]["mechanism_search"].get("context_allow_live") is True
            else ""
        )
    log: D + "/logs/d3_codon_context.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 {input.code:q} \
          --evidence {input.evidence:q} \
          --a10 {input.a10:q} \
          {params.genes} \
          --server {params.server:q} \
          --flank {params.flank} \
          --cache {params.cache:q} \
          --outdir {params.out:q} \
          {params.live} > {log:q} 2>&1
        """
# END d3_codon_context_v1

# BEGIN d3_readthrough_route_v1
rule d3_readthrough_route:
    input:
        scenarios=D + "/d0/scenarios.tsv.gz",
        code=SCRIPTS + "/d3_readthrough_route.py"
    output:
        routing=D + "/d3/readthrough_route/routing.tsv",
        manifest=D + "/d3/readthrough_route/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_route",
        mode=lambda wc: config.get("track2", {}).get(
            "rtdetective", {}
        ).get("mode", "auto"),
        scope=lambda wc: config.get("track2", {}).get(
            "rtdetective", {}
        ).get("hypothesis_scope", "primary_q_j")
    log:
        D + "/logs/d3_readthrough_route.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --scenarios {input.scenarios:q} \
          --outdir {params.out:q} \
          --mode {params.mode:q} \
          --scope {params.scope:q} > {log:q} 2>&1
        """

rule layer_d_readthrough_route:
    input:
        expand(
            OUT + "/{sample}/d3/readthrough_route/routing.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_route/manifest.json",
            sample=SAMPLES
        )
# END d3_readthrough_route_v1

# BEGIN d3_readthrough_multiview_v1
rule d3_readthrough_multiview:
    input:
        scenarios=D + "/d0/scenarios.tsv.gz",
        views=D + "/d0/candidate_model_views.tsv.gz",
        code=SCRIPTS + "/d3_readthrough_multiview.py",
        helper=SCRIPTS + "/d3_readthrough_route.py"
    output:
        routing=D + "/d3/readthrough_route_multiview/routing.tsv",
        witnesses=D + "/d3/readthrough_route_multiview/selection_witnesses.tsv",
        coverage=D + "/d3/readthrough_route_multiview/view_coverage.tsv",
        manifest=D + "/d3/readthrough_route_multiview/manifest.json"
    params:
        out=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_route_multiview",
        settings=lambda wc: __import__("json").dumps({
            "mode": config.get("track2", {}).get(
                "readthrough_selection", {}
            ).get("mode", "auto"),
            "views": config.get("track2", {}).get(
                "readthrough_selection", {}
            ).get("views", [
                "stan", "a10_informed",
                "a10_agnostic", "dependence_ensemble"
            ])
        }, sort_keys=True)
    log:
        D + "/logs/d3_readthrough_multiview.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --scenarios {input.scenarios:q} \
          --model-views {input.views:q} \
          --settings-json {params.settings:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule layer_d_readthrough_multiview:
    input:
        expand(
            OUT + "/{sample}/d3/readthrough_route_multiview/routing.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_route_multiview/selection_witnesses.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_route_multiview/view_coverage.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_route_multiview/manifest.json",
            sample=SAMPLES
        )
# END d3_readthrough_multiview_v1

# BEGIN fresh_codon_genomic_mapping_v2
rule d3_codon_genomic_mapping:
    input:
        context=D + "/d3/codon_context/codon_context.tsv",
        fasta=D + "/d3/codon_context/flanking_context.fasta",
        manifest=D + "/d3/codon_context/manifest.json",
        code=SCRIPTS + "/d3_codon_genomic_mapping.py"
    output:
        mapping=D + "/d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv",
        bed=D + "/d3/codon_context/genomic_mapping/stop_codon_lookup.bed",
        manifest=D + "/d3/codon_context/genomic_mapping/manifest.json"
    params:
        context=lambda wc: OUT + "/" + wc.sample + "/d3/codon_context",
        out=lambda wc: OUT + "/" + wc.sample + "/d3/codon_context/genomic_mapping",
        cache=lambda wc: OUT + "/" + wc.sample + "/d3/codon_context/.genomic_mapping_request_cache_v2",
        server=lambda wc: config.get("track2", {}).get("mechanism_search", {}).get("context_rest_server", "https://rest.ensembl.org"),
        live=lambda wc: "--allow-live" if config.get("track2", {}).get("mechanism_search", {}).get("mapping_allow_live") is True else ""
    log: D + "/logs/d3_codon_genomic_mapping.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} --context-dir {params.context:q} \
          --outdir {params.out:q} --request-cache {params.cache:q} \
          --server {params.server:q} {params.live} > {log:q} 2>&1
        """

rule layer_d_codon_genomic_mapping:
    input:
        expand(OUT + "/{sample}/d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/codon_context/genomic_mapping/stop_codon_lookup.bed", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/codon_context/genomic_mapping/manifest.json", sample=SAMPLES)
# END fresh_codon_genomic_mapping_v2

# BEGIN d3_readthrough_context_bridge_v1
def _t2_readthrough_cached_context_inputs(wc):
    directory = OUT + "/" + wc.sample + "/d3/codon_context"
    paths = [
        directory + "/codon_context.tsv",
        directory + "/manifest.json",
        directory + "/genomic_mapping/codon_genomic_mapping.tsv",
        directory + "/genomic_mapping/manifest.json"
    ]
    return paths

rule d3_readthrough_context_bridge:
    input:
        routing=D + "/d3/readthrough_route_multiview/routing.tsv",
        witnesses=D + "/d3/readthrough_route_multiview/selection_witnesses.tsv",
        coverage=D + "/d3/readthrough_route_multiview/view_coverage.tsv",
        route_manifest=D + "/d3/readthrough_route_multiview/manifest.json",
        cache=_t2_readthrough_cached_context_inputs,
        code=SCRIPTS + "/d3_readthrough_context_bridge.py"
    output:
        verification=D + "/d3/readthrough_context_bridge/verification.tsv",
        requests=D + "/d3/readthrough_context_bridge/context_requests.tsv",
        manifest=D + "/d3/readthrough_context_bridge/manifest.json"
    params:
        route=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_route_multiview",
        context=lambda wc: OUT + "/" + wc.sample + "/d3/codon_context",
        out=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_context_bridge",
        cache_json=lambda wc, input: __import__("json").dumps(list(input.cache)),
        assembly="GRCh38"
    log:
        D + "/logs/d3_readthrough_context_bridge.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --route-dir {params.route:q} \
          --context-dir {params.context:q} \
          --cache-inputs-json {params.cache_json:q} \
          --assembly {params.assembly:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule layer_d_readthrough_context_bridge:
    input:
        expand(
            OUT + "/{sample}/d3/readthrough_context_bridge/verification.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_context_bridge/context_requests.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/readthrough_context_bridge/manifest.json",
            sample=SAMPLES
        )
# END d3_readthrough_context_bridge_v1

# BEGIN d3_rtdetective_adopt_v1
def _t2_rtdetective_adoption_settings():
    value = config.get("track2", {}).get("rtdetective", {})
    return {
        "mode": value.get("mode", "auto"),
        "tracks": value.get("tracks", [
            "DAP_tga.txt", "Clitocine_tga.txt", "SRI_tga.txt"
        ]),
        "article_id": value.get("article_id", 23708901),
        "article_version": value.get("article_version", 2)
    }

def _t2_rtdetective_inventory():
    return config.get("track2", {}).get("rtdetective", {}).get(
        "inventory_dir",
        RES + "/supplemental/rtdetective_inventory_20261021T043409788319Z"
    )

def _t2_rtdetective_track_cache():
    return config.get("track2", {}).get("rtdetective", {}).get(
        "track_cache_dir",
        RES + "/supplemental/rtdetective_tracks_v2"
    )

def _t2_rtdetective_adoption_cache_inputs(wc):
    # Required producer outputs, not existing-file discovery.
    inventory = _t2_rtdetective_inventory()
    cache = _t2_rtdetective_track_cache()
    lookups = OUT + "/" + wc.sample + "/d3/rtdetective_lookup"
    paths = [inventory + "/article_metadata.json", inventory + "/manifest.json"]
    for name in _t2_rtdetective_adoption_settings()["tracks"]:
        if name not in {"DAP_tga.txt", "Clitocine_tga.txt", "SRI_tga.txt"}:
            raise ValueError("Unsupported RTDetective adoption track")
        paths.append(cache + "/" + name)
        directory = os.path.join(lookups, os.path.splitext(name)[0])
        paths.extend(os.path.join(directory, item) for item in
                     ("lookup_status.tsv", "exact_matches.tsv", "manifest.json"))
    return sorted(set(paths))

# BEGIN rtdetective_fresh_coordinate_lookup_v1
def _t2_rtdetective_lookup_track(wc):
    if wc.rtd_track not in {"DAP_tga", "Clitocine_tga", "SRI_tga"}:
        raise ValueError("Unsupported RTDetective track")
    return _t2_rtdetective_track_cache() + "/" + wc.rtd_track + ".txt"

rule d3_rtdetective_lookup:
    input:
        mapping=D + "/d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv",
        mapping_manifest=D + "/d3/codon_context/genomic_mapping/manifest.json",
        metadata=lambda wc: _t2_rtdetective_inventory() + "/article_metadata.json",
        inventory_manifest=lambda wc: _t2_rtdetective_inventory() + "/manifest.json",
        track=_t2_rtdetective_lookup_track,
        code=SCRIPTS + "/d3_rtdetective_lookup.py"
    output:
        status=D + "/d3/rtdetective_lookup/{rtd_track}/lookup_status.tsv",
        exact=D + "/d3/rtdetective_lookup/{rtd_track}/exact_matches.tsv",
        manifest=D + "/d3/rtdetective_lookup/{rtd_track}/manifest.json"
    params:
        inventory=lambda wc: _t2_rtdetective_inventory(),
        out=lambda wc: OUT + "/" + wc.sample + "/d3/rtdetective_lookup/" + wc.rtd_track
    wildcard_constraints:
        rtd_track="DAP_tga|Clitocine_tga|SRI_tga"
    log: D + "/logs/d3_rtdetective_lookup_{rtd_track}.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} --mapping-table {input.mapping:q} \
          --mapping-manifest {input.mapping_manifest:q} --track {input.track:q} \
          --inventory-dir {params.inventory:q} --outdir {params.out:q} > {log:q} 2>&1
        """
# END rtdetective_fresh_coordinate_lookup_v1

rule d3_rtdetective_adopt:
    input:
        verification=D + "/d3/readthrough_context_bridge/verification.tsv",
        requests=D + "/d3/readthrough_context_bridge/context_requests.tsv",
        bridge_manifest=D + "/d3/readthrough_context_bridge/manifest.json",
        mapping=D + "/d3/codon_context/genomic_mapping/codon_genomic_mapping.tsv",
        cache=_t2_rtdetective_adoption_cache_inputs,
        code=SCRIPTS + "/d3_rtdetective_adopt.py"
    output:
        predictions=D + "/d3/rtdetective_adoption/predictions.tsv",
        status=D + "/d3/rtdetective_adoption/adoption_status.tsv",
        manifest=D + "/d3/rtdetective_adoption/manifest.json"
    params:
        bridge=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_context_bridge",
        lookups=lambda wc: OUT + "/" + wc.sample + "/d3/rtdetective_lookup",
        inventory=lambda wc: _t2_rtdetective_inventory(),
        cache=lambda wc: _t2_rtdetective_track_cache(),
        out=lambda wc: OUT + "/" + wc.sample + "/d3/rtdetective_adoption",
        declared=lambda wc, input: __import__("json").dumps(list(input)),
        settings=lambda wc: __import__("json").dumps(
            _t2_rtdetective_adoption_settings(), sort_keys=True
        )
    log:
        D + "/logs/d3_rtdetective_adopt.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --bridge-dir {params.bridge:q} \
          --mapping-table {input.mapping:q} \
          --lookup-dir {params.lookups:q} \
          --inventory-dir {params.inventory:q} \
          --track-cache {params.cache:q} \
          --declared-inputs-json {params.declared:q} \
          --settings-json {params.settings:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule layer_d_rtdetective_adopt:
    input:
        expand(
            OUT + "/{sample}/d3/rtdetective_adoption/predictions.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/rtdetective_adoption/adoption_status.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/rtdetective_adoption/manifest.json",
            sample=SAMPLES
        )
# END d3_rtdetective_adopt_v1

# BEGIN d3_claim_ledger_v1
def _t2_claim_fulltext_dir(wc):
    configured = config.get("track2", {}).get(
        "evidence_ledger", {}
    ).get("fulltext_dirs", {}).get(wc.sample)
    if configured:
        return configured
    return (
        OUT + "/" + wc.sample
        + "/d3/approved_drug_review/ataluren/20261021T063427869589Z"
        + "/priority_fulltext_review_20261021T064023702342Z"
    )

def _t2_claim_bibliography_dir(wc):
    configured = config.get("track2", {}).get(
        "evidence_ledger", {}
    ).get("bibliography_dirs", {}).get(wc.sample)
    if configured:
        return configured
    return (
        OUT + "/" + wc.sample
        + "/d3/approved_drug_review/ataluren/20261021T063427869589Z"
        + "/bibliography_20261021T065312328652Z"
    )

def _t2_claim_seeds():
    return config.get("track2", {}).get(
        "evidence_ledger", {}
    ).get(
        "claim_seeds",
        RES + "/curation/ataluren_claim_seeds.v1.json"
    )

def _t2_claim_snapshot_inputs(wc):
    paths = []
    for directory in [
        _t2_claim_fulltext_dir(wc),
        _t2_claim_bibliography_dir(wc)
    ]:
        manifest_path = os.path.join(directory, "manifest.json")
        with open(manifest_path, encoding="utf-8") as stream:
            manifest = __import__("json").load(stream)
        paths.append(manifest_path)
        for name in manifest["output_sha256"]:
            if os.path.isabs(name) or ".." in name.replace("\\", "/").split("/"):
                raise ValueError("Unsafe snapshot manifest path")
            paths.append(os.path.join(directory, name))
    return sorted(set(paths))

rule d3_claim_ledger:
    input:
        predictions=D + "/d3/rtdetective_adoption/predictions.tsv",
        adoption_status=D + "/d3/rtdetective_adoption/adoption_status.tsv",
        adoption_manifest=D + "/d3/rtdetective_adoption/manifest.json",
        bridge_verification=D + "/d3/readthrough_context_bridge/verification.tsv",
        bridge_requests=D + "/d3/readthrough_context_bridge/context_requests.tsv",
        bridge_manifest=D + "/d3/readthrough_context_bridge/manifest.json",
        snapshots=_t2_claim_snapshot_inputs,
        inventory=lambda wc: _t2_rtdetective_inventory() + "/article_metadata.json",
        seeds=lambda wc: _t2_claim_seeds(),
        code=SCRIPTS + "/d3_claim_ledger.py"
    output:
        sources=D + "/d3/claim_ledger/sources.tsv",
        claims=D + "/d3/claim_ledger/claims.tsv",
        links=D + "/d3/claim_ledger/claim_hypothesis_links.tsv",
        review=D + "/d3/claim_ledger/review_queue.tsv",
        references=D + "/d3/claim_ledger/used_references.ris",
        manifest=D + "/d3/claim_ledger/manifest.json"
    params:
        adoption=lambda wc: OUT + "/" + wc.sample + "/d3/rtdetective_adoption",
        bridge=lambda wc: OUT + "/" + wc.sample + "/d3/readthrough_context_bridge",
        fulltext=_t2_claim_fulltext_dir,
        bibliography=_t2_claim_bibliography_dir,
        out=lambda wc: OUT + "/" + wc.sample + "/d3/claim_ledger"
    log:
        D + "/logs/d3_claim_ledger.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --bridge-dir {params.bridge:q} \
          --adoption-dir {params.adoption:q} \
          --fulltext-dir {params.fulltext:q} \
          --bibliography-dir {params.bibliography:q} \
          --inventory-metadata {input.inventory:q} \
          --seeds {input.seeds:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule layer_d_claim_ledger:
    input:
        expand(OUT + "/{sample}/d3/claim_ledger/sources.tsv", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/claim_ledger/claims.tsv", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/claim_ledger/claim_hypothesis_links.tsv", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/claim_ledger/review_queue.tsv", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/claim_ledger/used_references.ris", sample=SAMPLES),
        expand(OUT + "/{sample}/d3/claim_ledger/manifest.json", sample=SAMPLES)
# END d3_claim_ledger_v1

# BEGIN d3_mechanistic_ledger_v1
def _t2_primary_mechanistic_bundle():
    return config.get("track2", {}).get(
        "evidence_ledger", {}
    ).get(
        "mechanistic_bundle",
        RES + "/supplemental/primary_mechanisms/20261021T100841842445206Z"
    )

def _t2_primary_mechanistic_inputs(wc):
    directory = _t2_primary_mechanistic_bundle()
    manifest_path = os.path.join(directory, "manifest.json")
    with open(manifest_path, encoding="utf-8") as stream:
        manifest = __import__("json").load(stream)
    paths = [manifest_path]
    for name in manifest["output_sha256"]:
        if os.path.isabs(name) or ".." in name.replace("\\", "/").split("/"):
            raise ValueError("Unsafe mechanistic bundle path")
        paths.append(os.path.join(directory, name))
    return paths

rule d3_mechanistic_ledger:
    input:
        sources=D + "/d3/claim_ledger/sources.tsv",
        claims=D + "/d3/claim_ledger/claims.tsv",
        links=D + "/d3/claim_ledger/claim_hypothesis_links.tsv",
        review=D + "/d3/claim_ledger/review_queue.tsv",
        references=D + "/d3/claim_ledger/used_references.ris",
        manifest=D + "/d3/claim_ledger/manifest.json",
        augmentation=_t2_primary_mechanistic_inputs,
        code=SCRIPTS + "/d3_merge_mechanistic_ledger.py"
    output:
        sources=D + "/d3/claim_ledger_mechanistic/sources.tsv",
        snapshots=D + "/d3/claim_ledger_mechanistic/source_snapshots.tsv",
        claims=D + "/d3/claim_ledger_mechanistic/claims.tsv",
        links=D + "/d3/claim_ledger_mechanistic/claim_hypothesis_links.tsv",
        source_links=D + "/d3/claim_ledger_mechanistic/claim_source_links.tsv",
        review=D + "/d3/claim_ledger_mechanistic/review_queue.tsv",
        source_review=D + "/d3/claim_ledger_mechanistic/source_review_queue.tsv",
        references=D + "/d3/claim_ledger_mechanistic/used_references.ris",
        manifest=D + "/d3/claim_ledger_mechanistic/manifest.json"
    params:
        base=lambda wc: OUT + "/" + wc.sample + "/d3/claim_ledger",
        augmentation=lambda wc: _t2_primary_mechanistic_bundle(),
        out=lambda wc: OUT + "/" + wc.sample + "/d3/claim_ledger_mechanistic"
    log:
        D + "/logs/d3_mechanistic_ledger.log"
    shell:
        r"""
        mkdir -p $(dirname {log:q})
        python3 -B {input.code:q} \
          --base-ledger {params.base:q} \
          --augmentation {params.augmentation:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule layer_d_mechanistic_ledger:
    input:
        expand(
            OUT + "/{sample}/d3/claim_ledger_mechanistic/manifest.json",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/claim_ledger_mechanistic/claims.tsv",
            sample=SAMPLES
        ),
        expand(
            OUT + "/{sample}/d3/claim_ledger_mechanistic/used_references.ris",
            sample=SAMPLES
        )
# END d3_mechanistic_ledger_v1
