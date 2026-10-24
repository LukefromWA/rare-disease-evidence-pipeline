# Gene-only, unweighted Track 2 discovery; no case or model modules.
import os

CFG = os.path.abspath(os.environ["MVA_CONFIG"])
configfile: CFG

SCRIPTS = os.path.join(workflow.basedir, "scripts")
OUT = config["run"]["output_root"]
RES = config["track2"]["root"]
HIERARCHY = config["track2"]["reactome_hierarchy_dir"]
PYTHON = config["cli"]["python"]

rule all:
    default_target: True
    input:
        OUT + "/report/gene_summary.tsv",
        OUT + "/report/report.json"

rule prepare_gene_hypotheses:
    input:
        cfg=CFG,
        ready=RES + "/.ready",
        manifest=RES + "/manifests/resource_manifest.json",
        mapping=RES + "/derived/gene_protein_map.tsv.gz",
        code=SCRIPTS + "/mva_cli.py"
    output:
        hypotheses=OUT + "/inputs/gene_hypotheses.tsv",
        manifest=OUT + "/inputs/gene_input_manifest.json"
    shell:
        "{PYTHON:q} {input.code:q} _prepare --config {input.cfg:q}"

rule discover_gene_associations:
    input:
        hypotheses=rules.prepare_gene_hypotheses.output.hypotheses,
        gene_manifest=rules.prepare_gene_hypotheses.output.manifest,
        ready=RES + "/.ready",
        resource_manifest=RES + "/manifests/resource_manifest.json",
        drugs=RES + "/derived/drug_target_action_map.tsv.gz",
        pathways=RES + "/derived/reactome_gene_sets.tsv.gz",
        hierarchy=[
            HIERARCHY + "/" + name for name in (
                "checksums.sha256", "manifest.json",
                "ReactomePathways.txt", "ReactomePathwaysRelation.txt",
                "human_pathway_hierarchy.tsv"
            )
        ],
        code=SCRIPTS + "/d2_targets.R",
        helpers=sorted(__import__("glob").glob(
            os.path.join(SCRIPTS, "..", "lib", "*.R")
        )),
        cfg=CFG
    output:
        direct=OUT + "/discovery/target_universe.tsv.gz",
        audit=OUT + "/discovery/hypothesis_discovery_audit.tsv.gz",
        pathways=OUT + "/discovery/pathway_discovery.tsv.gz",
        witnesses=OUT + "/discovery/pathway_discovery_witnesses.tsv.gz",
        pathway_audit=OUT + "/discovery/pathway_discovery_audit.tsv.gz",
        manifest=OUT + "/discovery/manifest.json"
    params:
        out=OUT + "/discovery"
    log:
        OUT + "/logs/discovery.log"
    shell:
        r"""
        mkdir -p "$(dirname {log:q})"
        Rscript --vanilla {input.code:q} \
          --config {input.cfg:q} --sample gene_list \
          --mechanism {input.hypotheses:q} \
          --input_mode gene_only \
          --hierarchy_dir {HIERARCHY:q} \
          --outdir {params.out:q} > {log:q} 2>&1
        """

rule gene_exploration_report:
    input:
        hypotheses=rules.prepare_gene_hypotheses.output.hypotheses,
        direct=rules.discover_gene_associations.output.direct,
        pathways=rules.discover_gene_associations.output.pathways,
        manifest=rules.discover_gene_associations.output.manifest,
        code=SCRIPTS + "/mva_cli.py",
        cfg=CFG
    output:
        summary=OUT + "/report/gene_summary.tsv",
        report=OUT + "/report/report.json"
    shell:
        "{PYTHON:q} {input.code:q} _report --config {input.cfg:q}"
