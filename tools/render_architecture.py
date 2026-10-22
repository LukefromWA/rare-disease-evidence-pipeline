import argparse
from pathlib import Path
from diagrams import Diagram, Cluster, Edge, Node

parser = argparse.ArgumentParser()
parser.add_argument("--outdir", required=True)
args = parser.parse_args()
out = Path(args.outdir)
out.mkdir(parents=True, exist_ok=True)

def box(title, subtitle, color):
    return Node(
        title + "\n\n" + subtitle,
        shape="box",
        style="rounded,filled",
        fillcolor="#ffffff",
        color=color,
        fontcolor="#243247",
        penwidth="1.8",
        fontname="Helvetica",
        fontsize="12",
        labelloc="c",
        labeljust="c",
        fixedsize="false",
        margin="0.24,0.20",
        width="3.1",
        height="1.1",
    )

with Diagram(
    "From evidence to prioritization",
    filename=str(out / "architecture"),
    outformat=["svg", "png"],
    show=False,
    direction="TB",
    graph_attr={
        "bgcolor": "#ffffff",
        "fontname": "Helvetica",
        "fontsize": "23",
        "fontcolor": "#172b4d",
        "labelloc": "t",
        "labeljust": "c",
        "pad": "0.5",
        "nodesep": "0.65",
        "ranksep": "0.85",
        "splines": "spline",
    },
    edge_attr={
        "color": "#8290a5",
        "penwidth": "1.4",
        "arrowsize": "0.8",
    },
):
    with Cluster(
        "Public-cohort training",
        graph_attr={
            "bgcolor": "#f5f1ff",
            "color": "#d7ccf0",
            "style": "rounded",
            "margin": "24",
            "fontsize": "15",
        },
    ):
        cohort = box(
            "Public Phenopackets",
            "Version 0.1.27\nEligibility and cohort preparation",
            "#8250c4",
        )
        training = box(
            "B · External training",
            "Features, cross-validation\nand model fitting",
            "#8250c4",
        )
        retained = box(
            "Retained M3b predictor",
            "Three fold components\ncombined with equal weights",
            "#8250c4",
        )
        cohort >> training >> retained

    with Cluster(
        "Case analysis",
        graph_attr={
            "bgcolor": "#effbfc",
            "color": "#bddfe4",
            "style": "rounded",
            "margin": "24",
            "fontsize": "15",
        },
    ):
        inputs = box(
            "Case inputs",
            "Variant evidence\nand phenotype information",
            "#168fa5",
        )
        layer_a = box(
            "A · Evidence preparation",
            "Quality control, phenotype\nand network-based ranking",
            "#168fa5",
        )
        dependence = box(
            "Dependence analysis",
            "Evidence relationships\nand ranking sensitivity",
            "#2980b9",
        )
        inference = box(
            "C · Frozen inference",
            "Apply the retained model\nto upstream case features",
            "#159570",
        )

        inputs >> layer_a
        layer_a >> dependence
        layer_a >> inference

    retained >> inference

    resources = box(
        "E · Reference resources",
        "Versioned mappings,\nnetworks and supporting evidence",
        "#c58b29",
    )
    downstream = box(
        "D · Downstream analysis",
        "Candidate scenarios,\ndiscovery and evidence integration",
        "#159570",
    )
    reporting = box(
        "F · Reporting",
        "Generate tables, figures\nand analysis reports",
        "#617087",
    )

    inference >> downstream
    dependence >> downstream
    resources >> downstream
    downstream >> reporting

for extension in ("svg", "png"):
    path = out / f"architecture.{extension}"
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("Missing or empty diagram: " + extension)

print("Generated SVG and PNG previews.")
