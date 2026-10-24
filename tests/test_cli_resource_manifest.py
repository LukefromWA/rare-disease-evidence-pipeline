"""Synthetic integrity-reader tests; no case data or biological analysis."""
from pathlib import Path
import hashlib
import importlib.util
import json
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "mva_cli_under_test", REPO / "workflow/scripts/mva_cli.py"
)
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)

class ResourceManifestTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        (self.root / "manifests").mkdir()
        (self.root / "derived").mkdir()
        (self.root / ".ready").write_text("synthetic integrity fixture\n")
        self.entries = []
        for name in (
            "gene_protein_map.tsv.gz",
            "drug_target_action_map.tsv.gz",
            "reactome_gene_sets.tsv.gz",
        ):
            path = self.root / "derived" / name
            payload = ("synthetic integrity fixture: " + name).encode()
            path.write_bytes(payload)
            self.entries.append({
                "path": str(path),
                "sha256": hashlib.sha256(payload).hexdigest(),
            })

    def manifest(self, artifacts):
        (self.root / "manifests/resource_manifest.json").write_text(
            json.dumps({"scope": "core", "artifacts": artifacts})
        )

    def test_named_object(self):
        self.manifest({
            str(index): entry for index, entry in enumerate(self.entries)
        })
        cli.verify_resources(self.root)

    def test_list(self):
        self.manifest(self.entries)
        cli.verify_resources(self.root)

    def test_changed_payload_rejected(self):
        self.manifest(self.entries)
        Path(self.entries[0]["path"]).write_bytes(b"changed")
        with self.assertRaises(SystemExit):
            cli.verify_resources(self.root)

    def test_malformed_entry_rejected(self):
        self.manifest(["not-an-artifact-object"])
        with self.assertRaises(SystemExit):
            cli.verify_resources(self.root)

    def test_invalid_digest_rejected(self):
        self.entries[0]["sha256"] = "not-a-digest"
        self.manifest(self.entries)
        with self.assertRaises(SystemExit):
            cli.verify_resources(self.root)

    def test_missing_required_artifact_rejected(self):
        self.manifest(self.entries[:-1])
        with self.assertRaises(SystemExit):
            cli.verify_resources(self.root)

if __name__ == "__main__":
    unittest.main()
