# e_track2_resources.smk -- independent resource DAG (no private input, no C4). Ends at RES/.ready, written ONLY by the validator.
RAW = RES + "/raw"; DER = RES + "/derived"; GR = RES + "/graph"; MAN = RES + "/manifests"
SRC = {"signor": "signor/interactions.tsv.gz", "uniprot": "uniprot/reviewed_gene_map.tsv.gz", "reactome": "reactome/UniProt2Reactome.txt",
       "drug_targets": "drug_targets/drug_target_actions.tsv.gz", "pediatric_safety": "safety/pediatric_safety.tsv.gz"}


# Native metadata is configuration; payloads are read only by their jobs.
with open(REG) as _e_registry_handle:
    _E_SOURCES = yaml.safe_load(_e_registry_handle)["track2_resources"]["sources"]

def _e_value(source, key):
    value = _E_SOURCES[source].get(key)
    if not isinstance(value, str) or not value:
        raise ValueError("Resource build requires sources." + source + "." + key)
    return value

def _e_native(source):
    path = _e_value(source, "native_path")
    digest = _e_value(source, "native_sha256")
    if not os.path.isabs(path):
        raise ValueError("Resolve native_path to an absolute path: " + source)
    if not __import__("re").fullmatch("[0-9a-fA-F]{64}", digest):
        raise ValueError("Invalid native_sha256: " + source)
    return path

_e_reactome_rel = _e_value("reactome", "raw_path")
if os.path.isabs(_e_reactome_rel) or ".." in _e_reactome_rel.split("/"):
    raise ValueError("Reactome raw_path must stay inside the resource root")
_E_REACTOME_LOCAL = os.path.join(RES, _e_reactome_rel)


rule e02_prepare_uniprot:
    input:
        native=lambda wc: _e_native("uniprot"),
        registry=REG,
        code=SCRIPTS + "/e02_prepare_uniprot.r",
        helper=SCRIPTS + "/../lib/track2_resources.R"
    output:
        table=RAW + "/uniprot_reviewed_human.normalized.tsv.gz",
        mf=RAW + "/uniprot_reviewed_human.normalized.tsv.gz.resource_manifest.json",
        audit=RAW + "/uniprot.normalization_audit.tsv.gz"
    params:
        sha=lambda wc: _e_value("uniprot", "native_sha256"),
        url=lambda wc: _e_value("uniprot", "url"),
        version=lambda wc: _e_value("uniprot", "version"),
        license=lambda wc: _e_value("uniprot", "license")
    log: RES + "/logs/e02_prepare_uniprot.log"
    shell:
        "{RSCRIPT} {input.code:q} --native {input.native:q} "
        "--native_sha256 {params.sha:q} --out {output.table:q} "
        "--audit {output.audit:q} "
        "--source_url {params.url:q} --source_version {params.version:q} "
        "--license {params.license:q} > {log:q} 2>&1"

rule e03_prepare_signor:
    input:
        native=lambda wc: _e_native("signor"),
        registry=REG,
        code=SCRIPTS + "/e03_prepare_signor.r",
        helper=SCRIPTS + "/../lib/track2_resources.R"
    output:
        table=RAW + "/signor_human.normalized.tsv.gz",
        mf=RAW + "/signor_human.normalized.tsv.gz.resource_manifest.json",
        audit=RAW + "/signor.normalization_audit.tsv.gz", summary=RAW + "/signor.normalization_summary.json"
    params:
        sha=lambda wc: _e_value("signor", "native_sha256"),
        url=lambda wc: _e_value("signor", "url"),
        version=lambda wc: _e_value("signor", "version"),
        license=lambda wc: _e_value("signor", "license")
    log: RES + "/logs/e03_prepare_signor.log"
    shell:
        "{RSCRIPT} {input.code:q} --native {input.native:q} "
        "--native_sha256 {params.sha:q} --out {output.table:q} "
        "--audit {output.audit:q} --summary {output.summary:q} "
        "--source_url {params.url:q} --source_version {params.version:q} "
        "--license {params.license:q} > {log:q} 2>&1"

rule e04_stage_reactome:
    input:
        native=lambda wc: _e_native("reactome"),
        registry=REG,
        code=SCRIPTS + "/e04_stage_resource.R",
        helper=SCRIPTS + "/../lib/track2_resources.R"
    output:
        table=_E_REACTOME_LOCAL,
        mf=_E_REACTOME_LOCAL + ".resource_manifest.json"
    params:
        sha=lambda wc: _e_value("reactome", "native_sha256"),
        url=lambda wc: _e_value("reactome", "url"),
        version=lambda wc: _e_value("reactome", "version"),
        license=lambda wc: _e_value("reactome", "license")
    log: RES + "/logs/e04_stage_reactome.log"
    shell:
        "{RSCRIPT} {input.code:q} --native {input.native:q} "
        "--native_sha256 {params.sha:q} --out {output.table:q} "
        "--source_url {params.url:q} --source_version {params.version:q} "
        "--license {params.license:q} > {log:q} 2>&1"

# Offline preparation of the already acquired ChEMBL snapshot.
_T2_DRUG_SOURCE = yaml.safe_load(open(REG))["track2_resources"]["sources"]["drug_targets"]
_T2_CHEMBL_NATIVE = os.path.join(RES, _T2_DRUG_SOURCE["native_dir"])
_T2_CHEMBL_SUMS = os.path.join(_T2_CHEMBL_NATIVE, "native.sha256")
def _t2_chembl_inputs(wildcards):
    inputs = [
        os.path.join(_T2_CHEMBL_NATIVE, "download_manifest.json"),
        _T2_CHEMBL_SUMS,
    ]

    # Keep the missing inventory as a required input. Do not read it
    # while loading unrelated parts of the workflow.
    if not os.path.isfile(_T2_CHEMBL_SUMS):
        return inputs

    with open(_T2_CHEMBL_SUMS) as stream:
        for line in stream:
            if line.strip():
                inputs.append(
                    os.path.join(
                        _T2_CHEMBL_NATIVE,
                        line.rstrip("\n").split("  ", 1)[1]
                    )
                )
    return inputs

ruleorder: e06_prepare_drugs > e01_pull

rule e06_prepare_drugs:
    input:
        native=_t2_chembl_inputs,
        code=SCRIPTS + "/e13_build_drug_target_action_map.R",
        helper=SCRIPTS + "/../lib/track2_resources.R"
    output:
        raw=RAW + "/drug_targets/drug_target_actions.tsv.gz",
        mf=RAW + "/drug_targets/drug_target_actions.tsv.gz.resource_manifest.json",
        audit=RAW + "/drug_targets/drug_target_actions.tsv.gz.preparation_audit.tsv.gz",
        summary=RAW + "/drug_targets/drug_target_actions.tsv.gz.preparation_summary.json"
    params:
        native=_T2_CHEMBL_NATIVE,
        inventory_sha=lambda wc: _e_value("drug_targets", "native_inventory_sha256"),
        acquisition_sha=lambda wc: _e_value("drug_targets", "download_manifest_sha256")
    log: RES + "/logs/e06_prepare_drugs.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "printf '%s  %s\\n' {params.inventory_sha:q} {params.native:q}/native.sha256 "
        "| sha256sum -c - > {log:q} 2>&1; "
        "printf '%s  %s\\n' {params.acquisition_sha:q} {params.native:q}/download_manifest.json "
        "| sha256sum -c - >> {log:q} 2>&1; "
        "(cd {params.native:q} && sha256sum -c native.sha256) >> {log:q} 2>&1; "
        "{RSCRIPT} {input.code:q} --registry {REG:q} "
        "--raw {output.raw:q} --native_dir {params.native:q} "
        "--prepare_only true >> {log:q} 2>&1"


rule e00_registry:
    input:
        reg=ancient(REG),
        code=SCRIPTS + "/e00_validate_registry.R",
        uniprot=rules.e02_prepare_uniprot.output.table,
        uniprot_mf=rules.e02_prepare_uniprot.output.mf,
        signor=rules.e03_prepare_signor.output.table,
        signor_mf=rules.e03_prepare_signor.output.mf,
        reactome=rules.e04_stage_reactome.output.table,
        reactome_mf=rules.e04_stage_reactome.output.mf,
        drug_raw=rules.e06_prepare_drugs.output.raw,
        drug_mf=rules.e06_prepare_drugs.output.mf
    output: MAN + "/.registry_ok"
    log: RES + "/logs/e00.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --registry {input.reg:q} "
        "--root {RES:q} --out {output:q} --scope core > {log:q} 2>&1"

rule e01_pull:
    input: reg=ancient(REG)
    output: RAW + "/{rel}"
    wildcard_constraints:
        rel="(?:" + "|".join(__import__("re").escape(value) for value in SRC.values()) + ")"
    log: RES + "/logs/e01_{rel}.log"
    params: src=lambda wc: [k for k, v in SRC.items() if v == wc.rel][0]
    shell: "{RSCRIPT} {SCRIPTS}/e01_pull_resource.R --registry {input.reg} --root {RES} --source {params.src} > {log} 2>&1"

rule e10_gene_map:
    input: raw=rules.e02_prepare_uniprot.output.table, ok=rules.e00_registry.output
    output: m=DER + "/gene_protein_map.tsv.gz", lk=DER + "/uniprot_gene_lookup.tsv.gz", m_mf=DER + "/gene_protein_map.tsv.gz.resource_manifest.json", lk_mf=DER + "/uniprot_gene_lookup.tsv.gz.resource_manifest.json"
    log: RES + "/logs/e10.log"
    shell: "{RSCRIPT} {SCRIPTS}/e10_build_gene_protein_map.R --registry {REG} --raw {input.raw} --out_map {output.m} --out_lookup {output.lk} > {log} 2>&1"
rule e11_signed_edges:
    input: raw=rules.e03_prepare_signor.output.table, lk=rules.e10_gene_map.output.lk
    output: e=DER + "/signed_directed_edges.tsv.gz", c=DER + "/signed_network_conflicts.tsv.gz", a=DER + "/signed_network_drop_audit.tsv", e_mf=DER + "/signed_directed_edges.tsv.gz.resource_manifest.json", c_mf=DER + "/signed_network_conflicts.tsv.gz.resource_manifest.json"
    log: RES + "/logs/e11.log"
    shell: "{RSCRIPT} {SCRIPTS}/e11_build_signed_edges.R --registry {REG} --raw {input.raw} --lookup {input.lk} --out_edges {output.e} --out_conflicts {output.c} --out_audit {output.a} > {log} 2>&1"
rule e12_reactome:
    input: raw=rules.e04_stage_reactome.output.table, lk=rules.e10_gene_map.output.lk
    output: table=DER + "/reactome_gene_sets.tsv.gz", mf=DER + "/reactome_gene_sets.tsv.gz.resource_manifest.json"
    log: RES + "/logs/e12.log"
    shell: "{RSCRIPT} {SCRIPTS}/e12_build_reactome_gene_sets.R --registry {REG} --raw {input.raw} --lookup {input.lk} --out {output.table:q} > {log} 2>&1"
rule e13_drugs:
    input:
        raw=rules.e06_prepare_drugs.output.raw,
        raw_mf=rules.e06_prepare_drugs.output.mf,
        lk=rules.e10_gene_map.output.lk,
        code=SCRIPTS + "/e13_build_drug_target_action_map.R",
        helper=SCRIPTS + "/../lib/track2_resources.R"
    output:
        table=DER + "/drug_target_action_map.tsv.gz",
        mf=DER + "/drug_target_action_map.tsv.gz.resource_manifest.json"
    log: RES + "/logs/e13.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --registry {REG:q} "
        "--raw {input.raw:q} --lookup {input.lk:q} "
        "--out {output.table:q} > {log:q} 2>&1"

rule e14_safety:
    input: raw=RAW + "/safety/pediatric_safety.tsv.gz"
    output: DER + "/pediatric_safety_map.tsv.gz"
    log: RES + "/logs/e14.log"
    shell: "{RSCRIPT} {SCRIPTS}/e14_build_pediatric_safety_map.R --registry {REG} --raw {input.raw} --out {output} > {log} 2>&1"
rule e20_audit_graph:
    input: e=rules.e11_signed_edges.output.e, d=rules.e13_drugs.output.table
    output: g=GR + "/signed_directed_graph.rds", a=DER + "/signed_network_audit.tsv", m=GR + "/signed_directed_graph_metadata.json"
    log: RES + "/logs/e20.log"
    shell: "{RSCRIPT} {SCRIPTS}/e20_audit_signed_network.R --edges {input.e} --drugs {input.d} --out_graph {output.g} --out_audit {output.a} --out_meta {output.m} > {log} 2>&1"
rule e21_enablement:
    input: code=SCRIPTS + "/e21_signed_propagation_checks.R", g=GR + "/signed_directed_graph.rds"
    output: report=MAN + "/signed_propagation_checks.tsv"
    log: RES + "/logs/e21.log"
    params: sentinel=RES + "/.signed_propagation_enabled"
    shell: "{RSCRIPT} {SCRIPTS}/e21_signed_propagation_checks.R --graph {input.g} --sentinel {params.sentinel} --report {output.report} > {log} 2>&1"
rule e22_validate_bundle:
    input:
        g=rules.e20_audit_graph.output.g,
        e=rules.e11_signed_edges.output.e,
        c=rules.e11_signed_edges.output.c,
        a=rules.e20_audit_graph.output.a,
        genes=rules.e10_gene_map.output.m,
        lk=rules.e10_gene_map.output.lk,
        p=rules.e12_reactome.output.table,
        d=rules.e13_drugs.output.table,
        code=SCRIPTS + "/e22_validate_resource_bundle.R",
        registry_ok=rules.e00_registry.output
    output:
        ready=RES + "/.ready",
        manifest=MAN + "/resource_manifest.json",
        sums=MAN + "/checksums.sha256"
    log: RES + "/logs/e22.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --ready {output.ready:q} "
        "--manifest {output.manifest:q} --checksums {output.sums:q} "
        "--root {RES:q} --graph {input.g:q} --edges {input.e:q} "
        "--conflicts {input.c:q} --audit {input.a:q} "
        "--genes {input.genes:q} --lookup {input.lk:q} "
        "--pathways {input.p:q} --drugs {input.d:q} "
        "--scope core > {log:q} 2>&1"




# Complete core resource build, independent of any private case analysis.
rule resources:
    input:
        rules.e02_prepare_uniprot.output,
        rules.e03_prepare_signor.output,
        rules.e04_stage_reactome.output,
        rules.e06_prepare_drugs.output,
        rules.e00_registry.output,
        rules.e10_gene_map.output,
        rules.e11_signed_edges.output,
        rules.e12_reactome.output,
        rules.e13_drugs.output,
        rules.e20_audit_graph.output,
        rules.e22_validate_bundle.output
