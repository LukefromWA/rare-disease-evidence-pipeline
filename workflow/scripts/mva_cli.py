#!/usr/bin/env python3
"""CLI dispatcher; scientific calculations remain in the workflow scripts."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import shutil
import subprocess
import sys
import uuid
from datetime import datetime, timezone
import yaml

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"

def stop(message):
    raise SystemExit("STOP: " + message)

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for part in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()

def read_table(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_table(path, rows, columns):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)

def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")

def merge(base, overlay):
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            merge(base[key], value)
        else:
            base[key] = value
    return base


def normalize_path(value, base=None):
    """Accept Linux paths or Windows paths for a WSL/Linux process."""
    if value is None:
        stop("a path is required")
    text = os.fspath(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        text = text[1:-1]
    if not text or any(c in text for c in ("\x00", "\n", "\r", "\t")):
        stop("empty path or unsupported control character")
    if "," in text:
        stop("paths containing commas are unsupported by legacy manifest arguments")
    text = os.path.expandvars(os.path.expanduser(text))

    if re.match(r"^[A-Za-z]:[\\/]", text):
        if os.name == "nt":
            return Path(text).resolve()
        converter = shutil.which("wslpath")
        if converter:
            try:
                text = subprocess.check_output(
                    [converter, "-a", "-u", text],
                    text=True, stderr=subprocess.PIPE,
                ).strip()
            except subprocess.CalledProcessError as error:
                stop("Windows-path conversion failed: " + error.stderr.strip())
        else:
            windows = PureWindowsPath(text)
            mount = Path("/mnt") / windows.drive[0].lower()
            if not mount.is_dir():
                stop("Windows drive is not mounted; supply its Linux path")
            text = str(mount.joinpath(*windows.parts[1:]))
    elif re.match(r"^[A-Za-z]:", text):
        stop("drive-relative Windows paths are ambiguous; use C:\\folder\\file")
    elif text.startswith("\\\\"):
        parts = text[2:].split("\\")
        if len(parts) >= 3 and parts[0].lower() in {"wsl$", "wsl.localhost"}:
            current = os.environ.get("WSL_DISTRO_NAME")
            if not current or parts[1].casefold() != current.casefold():
                stop("use a Linux path for this distribution, not another WSL distribution")
            text = "/" + "/".join(parts[2:])
        else:
            stop("mount network shares first and supply their Linux mount path")

    path = Path(text)
    if not path.is_absolute():
        path = Path(base or Path.cwd()) / path
    return path.resolve()

def ask_yes_no(prompt, default=False):
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        answer = input(prompt + suffix).strip().casefold()
        if not answer:
            return default
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please enter y or n.")

def ask_path(prompt, default=None, kind="file"):
    while True:
        raw = input(prompt + (f" [{default}]" if default else "") + ": ").strip()
        try:
            path = normalize_path(raw or default)
            if kind == "file" and not path.is_file():
                print("File not found:", path)
                continue
            if kind == "dir" and not path.is_dir():
                print("Directory not found:", path)
                continue
            if kind == "new" and path.exists():
                print("Choose a new output directory; existing outputs will not be reused.")
                continue
            return path
        except SystemExit as error:
            print(error)

def vcf_samples(path, require_csq=False):
    path = normalize_path(path)
    if not path.name.endswith(".vcf.gz"):
        stop("current input routes require an indexed .vcf.gz file")
    has_csq = False
    try:
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("##INFO=<ID=CSQ,"):
                    has_csq = True
                if line.startswith("#CHROM"):
                    fields = line.rstrip("\r\n").split("\t")
                    samples = fields[9:]
                    if not samples or any(not s for s in samples):
                        stop("VCF has no valid sample columns")
                    if len(samples) != len(set(samples)):
                        stop("VCF sample names are duplicated")
                    if require_csq and not has_csq:
                        stop("selected VEP output has no CSQ header")
                    return samples
    except (OSError, UnicodeError) as error:
        stop("cannot read compressed VCF header: " + str(error))
    stop("VCF #CHROM header not found")

def hpo_terms(path):
    terms = []
    for line in normalize_path(path).read_text().splitlines():
        term = line.split("#", 1)[0].strip()
        if not term:
            continue
        if not re.fullmatch(r"HP:\d{7}", term):
            stop("HPO file must contain one HP:nnnnnnn term per line")
        if term not in terms:
            terms.append(term)
    if not terms:
        stop("at least one observed HPO term is required for this case route")
    return terms

def convert_windows_values(value):
    if isinstance(value, dict):
        return {k: convert_windows_values(v) for k, v in value.items()}
    if isinstance(value, list):
        return [convert_windows_values(v) for v in value]
    if isinstance(value, str) and (
        re.match(r"^[A-Za-z]:[\\/]", value) or value.startswith("\\\\")
    ):
        return str(normalize_path(value))
    return value

def create_case_run(
    profile_path, project_root, output, track, sample_label,
    selected_sample, hpo_path, annotated_vcf=None, annotated_index=None,
    metadata=None, raw_vcf=None, raw_index=None, fasta=None, fai=None,
    cache=None, cache_version=None, family_path=None, history=None,
    assembly=None,
):
    import copy
    profile_path, original = profile(profile_path)
    cfg = convert_windows_values(copy.deepcopy(original))
    root = normalize_path(project_root)
    out = normalize_path(output)
    if not root.is_dir():
        stop("project/data root does not exist")
    if out.exists():
        stop("new run directory already exists")
    if out == REPO or REPO in out.parents:
        stop("private runs must stay outside the source repository")
    if assembly != "GRCh38":
        stop("current annotation implementation supports declared human GRCh38 only")
    if track not in {"track1", "run"}:
        stop("unsupported new-case track")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample_label):
        stop("run sample label must use letters, numbers, dots, underscores or hyphens")
    if cfg.get("frozen_stan", {}).get("enabled") is not True:
        stop("profile must explicitly select the fresh frozen-Stan route")
    for field in ("external_root",):
        if not cfg.get("run", {}).get(field):
            stop("profile lacks run." + field)
    for field in ("root", "registry"):
        if not cfg.get("track2", {}).get(field):
            stop("profile lacks track2." + field)
    if not cfg.get("resources") or not cfg.get("qc"):
        stop("use a complete resource/model/analysis profile, including resources and QC")

    terms = hpo_terms(hpo_path)
    pv = {}
    assets = {}

    def required_file(value, label):
        if value is None:
            stop(label + " is required")
        path = normalize_path(value)
        if not path.is_file() or path.stat().st_size == 0:
            stop(label + " is absent or empty: " + str(path))
        assets[label] = str(path)
        return path

    if annotated_vcf is not None:
        vcf = required_file(annotated_vcf, "annotated_vcf")
        index = required_file(annotated_index, "annotated_index")
        meta = required_file(metadata, "vep_metadata")
        samples = vcf_samples(vcf, require_csq=True)
        text = meta.read_text()
        vep = cfg.get("vep", {})
        if vep.get("require_flags_in_metadata", True):
            absent = [
                flag for flag in vep.get("required_flags", [])
                if isinstance(flag, str) and flag.startswith("--") and flag not in text
            ]
            if absent:
                stop("VEP provenance lacks required flags: " + " ".join(absent))
        expected = vep.get("version_expected")
        version = re.search(r"(?i)vep[^0-9]{0,20}([0-9]+(?:\.[0-9]+)*)", text)
        if expected is not None and version and str(expected) not in version.group(1):
            stop("VEP metadata version differs from the profile")
        pv.update({
            "existing_biallelic_vep_vcf": str(vcf),
            "existing_biallelic_vep_vcf_index": str(index),
            "vep_run_metadata": str(meta),
        })
        mode = "existing_vep"
    else:
        vcf = required_file(raw_vcf, "raw_vcf")
        index = required_file(raw_index, "raw_index")
        reference = required_file(fasta, "reference_fasta")
        reference_index = required_file(fai, "reference_fasta_index")
        samples = vcf_samples(vcf)
        if cache is None or not normalize_path(cache).is_dir():
            stop("VEP cache directory is required")
        try:
            version = int(cache_version)
        except (TypeError, ValueError):
            stop("VEP cache version must be a positive integer")
        if version < 1:
            stop("VEP cache version must be positive")
        cfg.setdefault("vep", {})["cache_version"] = version
        cfg["vep"]["version_expected"] = version
        pv.update({
            "raw_vcf": str(vcf), "raw_vcf_index": str(index),
            "reference_fasta": str(reference),
            "reference_fasta_index": str(reference_index),
            "vep_cache_dir": str(normalize_path(cache)),
        })
        mode = "rerun_vep"

    if selected_sample not in samples:
        stop("selected sample is not present in the input VCF")
    if family_path is not None:
        family = normalize_path(family_path)
        if not family.is_file():
            stop("family-history phenotype file is missing")

    if track == "run" and history is not None:
        history = normalize_path(history)
        manifest_path = history / "manifest.json"
        if not manifest_path.is_file():
            stop("historical-support manifest is missing")
        historical = json.loads(manifest_path.read_text())
        if historical.get("version") != "A10_C4_union_Stan_v1_2":
            stop("unsupported historical-support manifest")
        if historical.get("current_build") != "GRCh38" or historical.get("legacy_build") != "GRCh38":
            stop("historical support does not declare GRCh38 compatibility")

    out.mkdir(parents=True, mode=0o700)
    inputs = out / "inputs"
    inputs.mkdir()
    observed = inputs / "proband_hpo.txt"
    observed.write_text("\n".join(terms) + "\n")
    family_copy = inputs / "family_history_hpo.txt"
    family_copy.write_text(
        normalize_path(family_path).read_text() if family_path is not None else ""
    )
    sample_file = inputs / "samples.tsv"
    write_table(sample_file, [{"sample_id": sample_label}], ["sample_id"])

    pv["proband_hpo"] = str(observed)
    pv["family_history_hpo"] = str(family_copy)
    cfg["private"] = pv
    cfg["samples_file"] = str(sample_file)
    cfg["run"]["run_id"] = "interactive_" + uuid.uuid4().hex
    cfg["run"]["mode"] = mode
    cfg["run"]["output_root"] = str(out / "results")
    cfg.setdefault("vcf", {})["sample_id"] = selected_sample
    cfg["vcf"]["sample_index"] = None
    cfg["vcf"].setdefault("required_format_fields", ["GT", "DP", "GQ", "AD"])
    t2 = cfg["track2"]
    t2["model_view_sources"] = {}
    t2["label_review_overrides"] = {}
    t2["frozen_union_dirs"] = (
        {sample_label: str(history)} if track == "run" and history is not None else {}
    )
    if isinstance(t2.get("mechanism_search"), dict):
        t2["mechanism_search"]["context_genes"] = {}
    if isinstance(t2.get("evidence_ledger"), dict):
        t2["evidence_ledger"]["fulltext_dirs"] = {}
        t2["evidence_ledger"]["bibliography_dirs"] = {}
    cfg.setdefault("cli", {}).update({
        "python": sys.executable,
        "project_root": str(root),
        "selected_track": track,
        "assembly_declared": "GRCh38",
        "source_profile": str(profile_path),
        "source_profile_sha256": sha(profile_path),
        "hpo_source": str(normalize_path(hpo_path)),
        "hpo_source_sha256": sha(normalize_path(hpo_path)),
        "annotation_assets": assets,
        "raw_annotation_execution_validated": False,
    })
    overlay = profile_path.parent / "config.private.yaml"
    if profile_path.name == "config.yaml" and overlay.is_file():
        cfg["cli"]["source_overlay_sha256"] = sha(overlay)
    config = out / "run.yaml"
    json_write(config, cfg)
    print("\nSaved private configuration:", config)
    print("Analysis outputs:", cfg["run"]["output_root"])
    print("Profile assumptions were retained; no model was trained.")
    return config

def interactive_case(track, cores):
    print("\nUse a complete resource/model/analysis profile.")
    print("An existing complete config is accepted; old private case inputs will be replaced.")
    settings = ask_path("Profile YAML")
    _, settings_cfg = profile(settings)
    root = ask_path("Project/data root for relative resource paths", kind="dir")
    if not ask_yes_no(
        "Use this profile's QC thresholds, gene priors and mechanism assumptions",
        default=False,
    ):
        stop("cancelled; review the profile first")

    have_raw = ask_yes_no("Do you have an indexed raw VCF", default=False)
    raw = raw_index = None
    if have_raw:
        raw = ask_path("Raw .vcf.gz")
        raw_index = ask_path("Raw VCF .tbi index", str(raw) + ".tbi")

    have_vep = ask_yes_no(
        "Do you already have VEP-annotated VCF output",
        default=False,
    )
    annotated = annotated_index = metadata = None
    if have_vep:
        annotated = ask_path("Existing VEP-annotated .vcf.gz")
        annotated_index = ask_path("Annotated VCF .tbi index", str(annotated) + ".tbi")
        metadata = ask_path("Original VEP command/provenance file")
        samples = vcf_samples(annotated, require_csq=True)
    else:
        if raw is None:
            stop("supply raw VCF or existing annotated output; gene-only Track 2 needs neither")
        samples = vcf_samples(raw)

    if not ask_yes_no("Are these human GRCh38 inputs", default=False):
        stop("this annotation implementation does not support another assembly")

    if len(samples) == 1:
        selected = samples[0]
        print("VCF sample selected:", selected)
    else:
        for index, sample in enumerate(samples, 1):
            print(f"{index}. {sample}")
        try:
            index = int(input("Select the VCF sample number: ").strip())
            if not 1 <= index <= len(samples):
                raise ValueError()
            selected = samples[index - 1]
        except ValueError:
            stop("invalid sample selection")

    sample_label = input(f"Run sample label [{selected}]: ").strip() or selected
    hpo = ask_path("Observed HPO text file, one HP:nnnnnnn term per line")
    family = (
        ask_path("Separate family-history phenotype file")
        if ask_yes_no("Supply separate family-history terms", False) else None
    )

    fasta = fai = cache = version = None
    if not have_vep:
        fasta = ask_path("GRCh38 reference FASTA")
        fai = ask_path("FASTA .fai index", str(fasta) + ".fai")
        cache = ask_path("Installed VEP cache directory", kind="dir")
        default_version = settings_cfg.get("vep", {}).get("cache_version", 116)
        version = input(f"VEP cache version [{default_version}]: ").strip() or default_version
        print("Raw annotation is an existing backend route, not yet execution-validated here.")

    history = None
    if track == "run":
        print("\nFresh connected analysis does not require historical case outputs.")
        if ask_yes_no("Include verified historical support for this same case", False):
            history = ask_path("Historical union directory", kind="dir")

    default_out = Path.home() / "mva_runs" / (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    )
    output = ask_path("New run directory", str(default_out), kind="new")
    config = create_case_run(
        settings, root, output, track, sample_label, selected, hpo,
        annotated_vcf=annotated, annotated_index=annotated_index,
        metadata=metadata, raw_vcf=raw, raw_index=raw_index,
        fasta=fasta, fai=fai, cache=cache, cache_version=version,
        family_path=family, history=history, assembly="GRCh38",
    )
    plan_only = ask_yes_no("Plan only; do not execute analysis", True)
    if not plan_only:
        if not have_vep and not ask_yes_no(
            "Execute the not-yet-validated raw annotation route", False
        ):
            stop("configuration saved; annotation execution cancelled")
        if not ask_yes_no(
            "Execute the configured analysis, including any configured public evidence queries",
            False,
        ):
            stop("configuration saved; execution cancelled")
    argv = [track, "--config", str(config), "--root", str(root), "--cores", str(cores)]
    if plan_only:
        argv.append("--dry-run")
    return argv


def profile(path):
    path = normalize_path(path)
    if not path.is_file():
        stop("configuration file does not exist: " + str(path))
    cfg = yaml.safe_load(path.read_text())
    overlay = path.parent / "config.private.yaml"
    if path.name == "config.yaml" and overlay.is_file():
        cfg = merge(cfg, yaml.safe_load(overlay.read_text()) or {})
    return path, cfg

def resolve(path, root):
    return normalize_path(path, base=root)

def gene_tokens(values):
    tokens = [
        token for value in values
        for token in re.split(r"[,;\s]+", value.strip()) if token
    ]
    if not tokens:
        stop("provide at least one gene")
    if any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", t)
           for t in tokens):
        stop("gene entries must be identifiers, not descriptions or shell commands")
    return tokens

def verify_resources(root):
    root = Path(root)
    if not (root / ".ready").is_file():
        stop("resource bundle is not ready")
    manifest = json.loads(
        (root / "manifests/resource_manifest.json").read_text()
    )
    if manifest.get("scope") != "core" or not manifest.get("artifacts"):
        stop("gene-only mode requires a validated core bundle")
    artifacts = manifest["artifacts"]
    if isinstance(artifacts, dict):
        entries = list(artifacts.values())
    elif isinstance(artifacts, list):
        entries = artifacts
    else:
        stop("resource artifacts must be an object or list")
    if not entries:
        stop("resource manifest has no artifact entries")

    recorded_paths = set()
    for item in entries:
        if not isinstance(item, dict):
            stop("resource artifact entry must be an object")
        filename = item.get("path")
        expected = item.get("sha256")
        if not isinstance(filename, str) or not filename:
            stop("resource artifact has no valid path")
        if not isinstance(expected, str) or not re.fullmatch(
            r"[0-9a-fA-F]{64}", expected
        ):
            stop("resource artifact has no valid SHA-256")
        path = Path(filename)
        if not path.is_absolute():
            stop("validated resource artifact paths must be absolute")
        if not path.is_file() or sha(path) != expected.lower():
            stop("validated resource changed: " + str(path))
        recorded_paths.add(path.resolve())

    required = {
        (root / "derived/gene_protein_map.tsv.gz").resolve(),
        (root / "derived/drug_target_action_map.tsv.gz").resolve(),
        (root / "derived/reactome_gene_sets.tsv.gz").resolve(),
    }
    if not required.issubset(recorded_paths):
        stop("required gene-discovery artifacts are absent from the manifest")
    return manifest

def prepare_genes(config):
    cfg_path, cfg = profile(config)
    out = Path(cfg["run"]["output_root"])
    resources = Path(cfg["track2"]["root"])
    verify_resources(resources)

    mapping_path = resources / "derived/gene_protein_map.tsv.gz"
    symbols = {
        row["gene_symbol"] for row in read_table(mapping_path)
        if row.get("gene_symbol", "").strip()
    }
    folded = {}
    for symbol in symbols:
        folded.setdefault(symbol.casefold(), set()).add(symbol)

    grouped = {}
    for token in cfg["gene_input"]["genes"]:
        matches = folded.get(token.casefold(), set())
        if token in symbols:
            gene, status = token, "exact_bundle_symbol"
        elif len(matches) == 1:
            gene, status = next(iter(matches)), "case_resolved_bundle_symbol"
        elif len(matches) > 1:
            stop("ambiguous gene spelling: " + token)
        else:
            gene, status = token, "not_in_bundle_gene_map"
        entry = grouped.setdefault(
            gene, {"status": status, "tokens": []}
        )
        entry["tokens"].append(token)

    rows = []
    for gene, entry in sorted(grouped.items()):
        hid = "manual_gene_" + hashlib.sha256(gene.encode()).hexdigest()[:24]
        rows.append({
            "hypothesis_id": hid,
            "scenario_id": hid + "__mechanism_unspecified",
            "gene": gene,
            "variant": "",
            "mechanism": "unknown",
            "direction": "uncertain_direction",
            "input_source": "manual_gene_list",
            "identity_status": entry["status"],
            "input_tokens_json": json.dumps(entry["tokens"]),
        })
    write_table(out / "inputs/gene_hypotheses.tsv", rows, list(rows[0]))
    json_write(out / "inputs/gene_input_manifest.json", {
        "mode": "gene_only_unweighted",
        "config_sha256": sha(cfg_path),
        "gene_map_sha256": sha(mapping_path),
        "hypotheses": len(rows),
        "identity_status_counts": {
            status: sum(r["identity_status"] == status for r in rows)
            for status in sorted({r["identity_status"] for r in rows})
        },
        "weights_supplied": False,
        "variant_coordinates_supplied": False,
        "mechanism_inferred": False,
        "input_gene_tokens_preserved": True,
    })

def gene_report(config):
    cfg_path, cfg = profile(config)
    out = Path(cfg["run"]["output_root"])
    hypotheses = read_table(out / "inputs/gene_hypotheses.tsv")
    direct = read_table(out / "discovery/target_universe.tsv.gz")
    pathways = read_table(out / "discovery/pathway_discovery.tsv.gz")
    discovery_manifest = json.loads(
        (out / "discovery/manifest.json").read_text()
    )
    for entry in discovery_manifest["outputs"].values():
        if sha(entry["path"]) != entry["sha256"]:
            stop("discovery output checksum mismatch")

    rows = []
    for hypothesis in hypotheses:
        hid = hypothesis["hypothesis_id"]
        rows.append({
            "hypothesis_id": hid,
            "gene": hypothesis["gene"],
            "identity_status": hypothesis["identity_status"],
            "direct_drugs": len({
                r["drug_id"] for r in direct if r["hypothesis_id"] == hid
            }),
            "exploratory_pathway_drugs": len({
                r["drug_id"] for r in pathways if r["hypothesis_id"] == hid
            }),
            "mechanism_status": "unknown",
            "direction_status": "not_assessed",
            "safety_status": "not_assessed",
            "weight_status": "not_supplied",
            "tiering_status": "not_performed_gene_only_exploration",
        })

    report_dir = out / "report"
    write_table(report_dir / "gene_summary.tsv", rows, list(rows[0]))
    json_write(report_dir / "report.json", {
        "mode": "gene_only_unweighted",
        "config_sha256": sha(cfg_path),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "hypotheses": len(rows),
        "scope": "Direct-target and shared-pathway discovery only",
        "no_gene_dropped_for_lack_of_mapping": True,
        "hypothesis_probabilities_created": False,
        "weighted_sensitivity_performed": False,
        "label_collection_performed": False,
        "mechanism_direction_and_safety_require_review": True,
        "disclaimer": (
            "Exploratory associations only. Not evidence of rescue, "
            "clinical efficacy, safety, or a treatment recommendation."
        ),
    })

def launch(snakefile, config, root, targets, cores, dry_run, allowed=None):
    executable = shutil.which("snakemake")
    if not executable:
        stop("snakemake is not on PATH; activate the tested runtime")
    env = os.environ.copy()
    env["MVA_CONFIG"] = str(config)
    env["MVA_ROOT"] = str(root)
    env.pop("MVA_LIB_DIR", None)
    cmd = [
        executable, "--snakefile", str(snakefile),
        "--directory", str(root), "--cores", str(cores),
        "--printshellcmds", "--rerun-incomplete",
    ]
    if dry_run:
        cmd.append("--dry-run")
    cmd.extend(targets)
    if allowed:
        cmd.extend(["--allowed-rules", *allowed])
    print("Workflow:", snakefile, flush=True)
    print("Configuration:", config, flush=True)
    return subprocess.run(cmd, env=env).returncode

def run(args):
    if args.cores < 1:
        stop("--cores must be positive")

    if args.command == "track2" and args.genes:
        if args.config or args.from_run:
            stop("--genes cannot be combined with --config or --from-run")
        if not all((args.out, args.resources, args.hierarchy)):
            stop("gene-only mode requires --out, --resources and --hierarchy")
        genes = gene_tokens(args.genes)
        out = normalize_path(args.out)
        if out == REPO or REPO in out.parents:
            stop("keep generated runs outside the source repository")
        if out.exists():
            stop("output path exists; use resume --config <run.yaml>")
        resources = normalize_path(args.resources)
        hierarchy = normalize_path(args.hierarchy)
        for path in (
            resources / ".ready",
            resources / "manifests/resource_manifest.json",
            hierarchy / "manifest.json",
            hierarchy / "checksums.sha256",
        ):
            if not path.is_file():
                stop("required resource input missing: " + str(path))
        out.mkdir(parents=True, mode=0o700)
        config = out / "run.yaml"
        json_write(config, {
            "run": {
                "run_id": "gene_only_" + uuid.uuid4().hex,
                "mode": "gene_only_unweighted",
                "seed": 20260912,
                "output_root": str(out),
            },
            "track2": {
                "root": str(resources),
                "resource_scope": "core",
                "reactome_hierarchy_dir": str(hierarchy),
            },
            "gene_input": {"genes": genes},
            "cli": {"python": sys.executable},
        })
        return launch(
            WF / "track2_genes.Snakefile", config, out,
            ["all"], args.cores, args.dry_run,
        )

    config_arg = args.config
    if args.command == "track2" and args.from_run:
        if config_arg:
            stop("choose --config or --from-run, not both")
        workspace = normalize_path(args.from_run)
        candidates = [workspace / "run.yaml", workspace / "config.yaml"]
        selected = next((p for p in candidates if p.is_file()), None)
        if selected is None:
            stop("previous case workspace contains neither run.yaml nor config.yaml: "
                 + str(workspace))
        config_arg = str(selected)
    if not config_arg:
        stop("provide --config, or use track2 --genes with resource paths")

    config, cfg = profile(config_arg)
    if cfg["run"].get("mode") == "gene_only_unweighted":
        return launch(
            WF / "track2_genes.Snakefile", config,
            Path(cfg["run"]["output_root"]),
            ["all"], args.cores, args.dry_run,
        )
    if args.command == "resume":
        stop("for a master configuration choose track1, track2, or run explicitly")
    if cfg.get("frozen_stan", {}).get("enabled") is not True:
        stop("this CLI version supports the explicit fresh frozen-Stan master route")

    root = normalize_path(
        args.root or cfg.get("cli", {}).get("project_root")
        or os.environ.get("MVA_ROOT", ".")
    )
    if not root.is_dir():
        stop("project root does not exist")

    targets = ["track1"] if args.command == "track1" else ["all"]
    allowed = None
    if args.command == "track2":
        samples_file = resolve(cfg["samples_file"], root)
        samples = [
            r["sample_id"] for r in read_table(samples_file)
        ]
        output = resolve(cfg["run"]["output_root"], root)
        for sample in samples:
            for relative in (
                "c4_stan_union/manifest.json",
                "c4_stan_union/current_Stan_final_hypothesis_set.tsv",
                "a10/ranked_candidates.tsv.gz",
                "a11_review/calculation/run_provenance.json",
                "a11_review/calculation/bayesian_dependence_sensitivity_ensemble_consensus.csv",
                "d7/gene_key_coverage.json",
                "d7/symbol_unresolved_variant_evidence.tsv.gz",
            ):
                if not (output / sample / relative).is_file():
                    stop("existing Track 1 handoff missing: " + str(output / sample / relative))
        verify_resources(resolve(cfg["track2"]["root"], root))
        targets = ["layer_d"]
        allowed = ["layer_d"] + re.findall(
            r"(?m)^\s*rule\s+(d[A-Za-z0-9_]*)\s*:",
            (WF / "rules/d_track2.smk").read_text(),
        )
    return launch(
        WF / "Snakefile", config, root, targets,
        args.cores, args.dry_run, allowed,
    )

def wizard():
    if os.name == "nt":
        stop("run the workflow under WSL; Windows path syntax is accepted there")
    print("\nMVA research workflow")
    print("1. Track 1\n2. Track 2\n3. Track 1 -> Track 2")
    choice = input("Choose 1, 2 or 3: ").strip()
    if choice not in {"1", "2", "3"}:
        stop("invalid track selection")
    try:
        cores = int(input("Cores [1]: ").strip() or "1")
        if cores < 1:
            raise ValueError()
    except ValueError:
        stop("cores must be a positive integer")
    if choice != "2":
        return interactive_case("track1" if choice == "1" else "run", cores)

    print("\n1. Existing pipeline run\n2. Candidate genes")
    method = input("Choose 1 or 2: ").strip()
    if method == "1":
        previous = ask_path("Previous case workspace", kind="dir")
        root = ask_path("Original project/data root", kind="dir")
        argv = [
            "track2", "--from-run", str(previous),
            "--root", str(root), "--cores", str(cores),
        ]
    elif method == "2":
        genes = input("Genes, separated by commas or spaces: ").strip()
        gene_tokens([genes])
        resources = ask_path("Validated core resource directory", kind="dir")
        hierarchy = ask_path("Matching Reactome hierarchy directory", kind="dir")
        output = ask_path("New output directory", kind="new")
        argv = [
            "track2", "--genes", genes,
            "--resources", str(resources), "--hierarchy", str(hierarchy),
            "--out", str(output), "--cores", str(cores),
        ]
        print("Gene-only mode performs unweighted discovery, not D0-D7 research tiering.")
    else:
        stop("invalid Track 2 input mode")

    if ask_yes_no("Plan only; do not execute analysis", True):
        argv.append("--dry-run")
    elif not ask_yes_no("Execute this configured Track 2 route", False):
        stop("execution cancelled")
    return argv


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] in {"_prepare", "_report"}:
        parser = argparse.ArgumentParser()
        parser.add_argument("operation")
        parser.add_argument("--config", required=True)
        args = parser.parse_args(argv)
        (prepare_genes if args.operation == "_prepare" else gene_report)(args.config)
        return 0
    parser = argparse.ArgumentParser(prog="mva")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("track1", "track2", "run", "resume"):
        p = sub.add_parser(name)
        p.add_argument("--config")
        p.add_argument("--root")
        p.add_argument("--cores", type=int, default=1)
        p.add_argument("--dry-run", action="store_true")
        if name == "track2":
            p.add_argument("--genes", nargs="+")
            p.add_argument("--from-run")
            p.add_argument("--resources")
            p.add_argument("--hierarchy")
            p.add_argument("--out")
    return run(parser.parse_args(argv or wizard()))

if __name__ == "__main__":
    sys.exit(main())
