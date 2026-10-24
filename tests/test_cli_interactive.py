"""Synthetic CLI configuration tests; no VEP execution or patient data."""
from pathlib import Path
import csv
import gzip
import importlib.util
import json
import os
import tempfile
import unittest
from unittest.mock import patch

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "mva_interactive_under_test", REPO / "workflow/scripts/mva_cli.py"
)
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)

class InteractiveTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.vcf = self.root / "input with spaces.vcf.gz"
        with gzip.open(self.vcf, "wt") as handle:
            handle.write(
                "##fileformat=VCFv4.2\n"
                '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Feature">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_1\n"
            )
        self.index = Path(str(self.vcf) + ".tbi")
        self.index.write_bytes(b"synthetic configuration fixture")
        self.meta = self.root / "vep metadata.txt"
        self.meta.write_text("VEP 116 --hgvs --symbol; synthetic configuration test\n")
        self.hpo = self.root / "phenotypes.txt"
        self.hpo.write_text("# synthetic configuration fixture\nHP:0000001\n")
        self.profile = self.root / "profile.yaml"
        self.profile.write_text(json.dumps({
            "run": {"external_root": "external", "output_root": "old"},
            "resources": {"hpo_obo": "resources/hp.obo"},
            "qc": {"min_gq": 20},
            "private": {"old_case_input": "must-not-survive"},
            "track2": {
                "root": "resources/core", "registry": "resources/registry.yaml",
                "label_review_overrides": {"OLD_CASE": "old.tsv"},
                "frozen_union_dirs": {"OLD_CASE": "old_union"},
            },
            "frozen_stan": {"enabled": True},
            "vep": {
                "required_flags": ["--hgvs", "--symbol"],
                "version_expected": 116,
                "cache_version": 116,
            },
        }))

    def create(self, **changes):
        options = dict(
            profile_path=self.profile, project_root=self.root,
            output=self.root / "new run", track="track1",
            sample_label="SYNTHETIC_1", selected_sample="SYNTHETIC_1",
            hpo_path=self.hpo, annotated_vcf=self.vcf,
            annotated_index=self.index, metadata=self.meta,
            assembly="GRCh38",
        )
        options.update(changes)
        return cli.create_case_run(**options)

    def test_linux_path_and_quotes(self):
        self.assertEqual(
            cli.normalize_path('"' + str(self.vcf) + '"'), self.vcf.resolve()
        )

    def test_windows_path_uses_converter(self):
        with patch.object(cli.shutil, "which", return_value="/usr/bin/wslpath"):
            with patch.object(
                cli.subprocess, "check_output",
                return_value="/mnt/c/Data Folder/sample.vcf.gz\n",
            ) as convert:
                result = cli.normalize_path(r'"C:\Data Folder\sample.vcf.gz"')
        self.assertEqual(str(result), "/mnt/c/Data Folder/sample.vcf.gz")
        self.assertIn(r"C:\Data Folder\sample.vcf.gz", convert.call_args.args[0])

    def test_drive_relative_path_rejected(self):
        with self.assertRaises(SystemExit):
            cli.normalize_path(r"C:relative\sample.vcf.gz")

    def test_yes_no_reprompts(self):
        with patch("builtins.input", side_effect=["maybe", "yes"]):
            self.assertTrue(cli.ask_yes_no("Test", False))

    def test_existing_vep_configuration(self):
        config = self.create()
        value = json.loads(config.read_text())
        self.assertEqual(value["run"]["mode"], "existing_vep")
        self.assertNotIn("old_case_input", value["private"])
        self.assertEqual(value["track2"]["frozen_union_dirs"], {})
        self.assertEqual(value["track2"]["label_review_overrides"], {})
        self.assertEqual(value["vcf"]["sample_id"], "SYNTHETIC_1")
        self.assertEqual(value["cli"]["project_root"], str(self.root))
        self.assertEqual(
            Path(value["private"]["proband_hpo"]).read_text(), "HP:0000001\n"
        )

    def test_missing_csq_rejected(self):
        raw = self.root / "raw.vcf.gz"
        with gzip.open(raw, "wt") as handle:
            handle.write(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC_1\n"
            )
        with self.assertRaises(SystemExit):
            cli.vcf_samples(raw, require_csq=True)

    def test_no_fake_metadata(self):
        with self.assertRaises(SystemExit):
            self.create(metadata=None)
        self.assertFalse((self.root / "new run").exists())

    def test_connected_requires_history(self):
        with self.assertRaises(SystemExit):
            self.create(track="run", history=None)
        self.assertFalse((self.root / "new run").exists())

    def test_raw_configuration(self):
        fasta = self.root / "reference.fa"
        fasta.write_text(">synthetic\nA\n")
        fai = self.root / "reference.fa.fai"
        fai.write_text("synthetic fixture\n")
        cache = self.root / "cache"
        cache.mkdir()
        config = self.create(
            annotated_vcf=None, annotated_index=None, metadata=None,
            raw_vcf=self.vcf, raw_index=self.index,
            fasta=fasta, fai=fai, cache=cache, cache_version=116,
        )
        value = json.loads(config.read_text())
        self.assertEqual(value["run"]["mode"], "rerun_vep")
        self.assertEqual(value["vep"]["cache_version"], 116)
        self.assertNotIn("existing_biallelic_vep_vcf", value["private"])
        self.assertFalse(value["cli"]["raw_annotation_execution_validated"])

if __name__ == "__main__":
    unittest.main()
