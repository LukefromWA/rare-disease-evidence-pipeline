# e_track2_resources.smk -- independent resource DAG (no private input, no C4). Ends at RES/.ready, written ONLY by the validator.
RAW = RES + "/raw"; DER = RES + "/derived"; GR = RES + "/graph"; MAN = RES + "/manifests"
SRC = {"signor": "signor/interactions.tsv.gz", "uniprot": "uniprot/reviewed_gene_map.tsv.gz", "reactome": "reactome/UniProt2Reactome.txt",
       "drug_targets": "drug_targets/drug_target_actions.tsv.gz", "pediatric_safety": "safety/pediatric_safety.tsv.gz"}

rule e00_registry:
    input:
        reg=ancient(REG),
        code=SCRIPTS + "/e00_validate_registry.R",
        drug_raw=RAW + "/drug_targets/drug_target_actions.tsv.gz",
        drug_mf=RAW + "/drug_targets/drug_target_actions.tsv.gz.resource_manifest.json"
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
    input: raw=RES + "/raw/uniprot_reviewed_human.normalized.tsv.gz", ok=MAN + "/.registry_ok"
    output: m=DER + "/gene_protein_map.tsv.gz", lk=DER + "/uniprot_gene_lookup.tsv.gz"
    log: RES + "/logs/e10.log"
    shell: "{RSCRIPT} {SCRIPTS}/e10_build_gene_protein_map.R --registry {REG} --raw {input.raw} --out_map {output.m} --out_lookup {output.lk} > {log} 2>&1"
rule e11_signed_edges:
    input: raw=RES + "/raw/signor_human.normalized.tsv.gz", lk=DER + "/uniprot_gene_lookup.tsv.gz"
    output: e=DER + "/signed_directed_edges.tsv.gz", c=DER + "/signed_network_conflicts.tsv.gz", a=DER + "/signed_network_drop_audit.tsv"
    log: RES + "/logs/e11.log"
    shell: "{RSCRIPT} {SCRIPTS}/e11_build_signed_edges.R --registry {REG} --raw {input.raw} --lookup {input.lk} --out_edges {output.e} --out_conflicts {output.c} --out_audit {output.a} > {log} 2>&1"
rule e12_reactome:
    input: raw=RES + "/native/reactome_20261020T100433Z/UniProt2Reactome_All_Levels.txt", lk=DER + "/uniprot_gene_lookup.tsv.gz"
    output: DER + "/reactome_gene_sets.tsv.gz"
    log: RES + "/logs/e12.log"
    shell: "{RSCRIPT} {SCRIPTS}/e12_build_reactome_gene_sets.R --registry {REG} --raw {input.raw} --lookup {input.lk} --out {output} > {log} 2>&1"
rule e13_drugs:
    input:
        raw=RAW + "/drug_targets/drug_target_actions.tsv.gz",
        raw_mf=RAW + "/drug_targets/drug_target_actions.tsv.gz.resource_manifest.json",
        lk=DER + "/uniprot_gene_lookup.tsv.gz",
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
    input: e=DER + "/signed_directed_edges.tsv.gz", d=DER + "/drug_target_action_map.tsv.gz"
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
        g=GR + "/signed_directed_graph.rds",
        e=DER + "/signed_directed_edges.tsv.gz",
        c=DER + "/signed_network_conflicts.tsv.gz",
        a=DER + "/signed_network_audit.tsv",
        genes=DER + "/gene_protein_map.tsv.gz",
        lk=DER + "/uniprot_gene_lookup.tsv.gz",
        p=DER + "/reactome_gene_sets.tsv.gz",
        d=DER + "/drug_target_action_map.tsv.gz",
        code=SCRIPTS + "/e22_validate_resource_bundle.R",
        registry_ok=MAN + "/.registry_ok"
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
    params: native=_T2_CHEMBL_NATIVE
    log: RES + "/logs/e06_prepare_drugs.log"
    shell:
        "mkdir -p $(dirname {log:q}); "
        "{RSCRIPT} {input.code:q} --registry {REG:q} "
        "--raw {output.raw:q} --native_dir {params.native:q} "
        "--prepare_only true > {log:q} 2>&1"
