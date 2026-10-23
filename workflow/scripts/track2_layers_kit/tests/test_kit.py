"""Self-test for the kit. Offline, synthetic data only:  python3 tests/test_kit.py"""
import contextlib, gzip, importlib.util, io, json, tempfile, unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def load(rel, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
bridge = load("layer2_handoff/make_d0_input.py", "bridge")
labels = load("layer3_evidence/fetch_label_evidence.py", "labels")
export = load("tools/public_export_check.py", "export")

def tsv(path, header, rows):
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "wt", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")

class Bridge(unittest.TestCase):
    def setUp(self):
        self.t = tempfile.TemporaryDirectory(); self.d = Path(self.t.name)
        tsv(self.d / "a10.tsv", ["gene_symbol", "best_candidate_id", "posterior_informed", "best_IMPACT"],
            [["G1", "v1", 0.6, "HIGH"], ["G2", "v2", 0.3, "MODERATE"], ["G3", "v3", 0.1, "LOW"]])
        tsv(self.d / "em.tsv", ["candidate_id", "consequence"], [["v1", "stop_gained"], ["v2", "missense_variant"], ["v3", "intron_variant"]])
        tsv(self.d / "prior.tsv", ["consequence_class", "mechanism", "pm"],
            [["default", "loss_of_function", 0.7], ["default", "unknown", 0.3]])
    def run_bridge(self, *extra, out="out/h.tsv.gz"):
        args = ["--a10", str(self.d / "a10.tsv"), "--evidence-matrix", str(self.d / "em.tsv"), "--out", str(self.d / out)] + list(extra)
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()): bridge.main(args)
        return [l.split("\t") for l in gzip.open(self.d / out, "rt").read().splitlines()[1:]]
    def test_full_declared(self):
        rows = self.run_bridge("--prior-mode", "declared", "--mechanism-prior", str(self.d / "prior.tsv"))
        self.assertEqual(len(rows), 6); self.assertAlmostEqual(sum({r[1]: float(r[2]) for r in rows}.values()), 1)
    def test_top_k_renormalizes_and_reports_truncation(self):
        rows = self.run_bridge("--selection", "top_k", "--top-k", "2", "--prior-mode", "uniform_neutral", out="o2/h.tsv.gz")
        self.assertEqual({r[1] for r in rows}, {"v1", "v2"}); self.assertAlmostEqual(sum({r[1]: float(r[2]) for r in rows}.values()), 1)
        m = json.loads((self.d / "o2/h.manifest.json").read_text()); self.assertAlmostEqual(m["selection"]["truncated_mass_fraction"], 0.1)
        self.assertFalse(m["prior"]["is_declared_biological_prior"])
    def test_placeholder_prior_fails_loudly(self):
        tsv(self.d / "bad.tsv", ["consequence_class", "mechanism", "pm"], [["default", "loss_of_function", "TBD"]])
        with self.assertRaises(ValueError): self.run_bridge("--prior-mode", "declared", "--mechanism-prior", str(self.d / "bad.tsv"), out="o3/h.tsv.gz")
    def test_refuses_overwrite(self):
        self.run_bridge("--prior-mode", "uniform_neutral")
        with self.assertRaises(ValueError): self.run_bridge("--prior-mode", "uniform_neutral")

class Labels(unittest.TestCase):
    def fake(self, log):
        def transport(url):
            log.append(url)
            if "ZZZ" in url: return 404, b'{"error":{"code":"NOT_FOUND"}}'
            hits = {"ASPIRIN": [{"id": "L1", "set_id": "S1", "version": "1", "openfda": {"route": ["ORAL"]}, "pediatric_use": ["text"], "boxed_warning": ["warn"]}]}
            for k, v in hits.items():
                if k in url: return 200, json.dumps({"meta": {"last_updated": "2026-10-16", "results": {"total": len(v)}}, "results": v}).encode()
            return 404, b'{"error":{"code":"NOT_FOUND"}}'
        return transport
    def test_lookup_cache_overrides(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); tsv(d / "drugs.tsv", ["drug_id", "drug_name"], [["C1", "Aspirin"], ["C1", "Aspirin"], ["C2", "ZZZ"], ["C3", ""]])
            log = []; common = ["--drugs", str(d / "drugs.tsv"), "--cache-dir", str(d / "cache"), "--delay", "0"]
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                labels.main(common + ["--out-dir", str(d / "o1")], transport=self.fake(log))
                n = len(log); labels.main(common + ["--out-dir", str(d / "o2")], transport=self.fake(log))
            self.assertEqual(len(log), n)  # second run is fully cached
            st = {r["drug_id"]: r["status"] for r in labels.read_tsv(d / "o1/drug_lookup_status.tsv")[1]}
            self.assertEqual(st, {"C1": "matched_exact", "C2": "no_match", "C3": "missing_name"})
            m = {r["drug_id"]: r for r in labels.read_tsv(d / "o1/pediatric_safety_map.tsv.gz")[1]}
            self.assertEqual((m["C1"]["pediatric_evidence_status"], m["C1"]["route_feasible"], m["C1"]["review_state"]), ("unknown", "", "unreviewed"))
            tsv(d / "ov.tsv", ["drug_id", "pediatric_evidence_status", "contraindication_status", "route_feasible", "evidence_label_id", "reviewer_note"],
                [["C1", "limited", "none", "TRUE", "L1", "reviewed label text"]])
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                labels.main(common + ["--out-dir", str(d / "o3"), "--review-overrides", str(d / "ov.tsv")], transport=self.fake(log))
            self.assertEqual({r["drug_id"]: r["review_state"] for r in labels.read_tsv(d / "o3/pediatric_safety_map.tsv.gz")[1]}["C1"], "reviewed")
            tsv(d / "bad.tsv", ["drug_id", "pediatric_evidence_status", "contraindication_status", "route_feasible", "evidence_label_id", "reviewer_note"],
                [["C1", "limited", "none", "TRUE", "NOT_A_LABEL", "x"]])
            with self.assertRaises(ValueError):
                with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                    labels.main(common + ["--out-dir", str(d / "o4"), "--review-overrides", str(d / "bad.tsv")], transport=self.fake(log))
    def test_max_drugs_guard(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); tsv(d / "drugs.tsv", ["drug_id", "drug_name"], [["A", "x"], ["B", "y"]])
            with self.assertRaises(ValueError):
                labels.main(["--drugs", str(d / "drugs.tsv"), "--cache-dir", str(d / "c"), "--out-dir", str(d / "o"), "--max-drugs", "1"], transport=lambda u: (404, b"NOT_FOUND"))

class Export(unittest.TestCase):
    def test_blocks_and_exports(self):
        with tempfile.TemporaryDirectory() as td:
            d = Path(td); (d / "src").mkdir(); (d / "src/ok.md").write_text("clean\n"); (d / "src/bad.md").write_text("see C:\\Users\\someone\\x and PATIENT42\n")
            (d / "terms.txt").write_text("PATIENT42\n")
            with contextlib.redirect_stdout(io.StringIO()) as buf:
                rc = export.main(["--paths", str(d / "src"), "--terms", str(d / "terms.txt"), "--export-dir", str(d / "pub")])
            self.assertEqual(rc, 1); self.assertFalse((d / "pub").exists()); self.assertNotIn("someone", buf.getvalue()); self.assertNotIn("PATIENT42", buf.getvalue())
            (d / "src/bad.md").unlink()
            with contextlib.redirect_stdout(io.StringIO()):
                rc = export.main(["--paths", str(d / "src"), "--export-dir", str(d / "pub")])
            self.assertEqual(rc, 0); self.assertTrue((d / "pub/src/ok.md").exists())

if __name__ == "__main__":
    unittest.main(verbosity=2)
