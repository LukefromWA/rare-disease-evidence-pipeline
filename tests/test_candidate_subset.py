#!/usr/bin/env python3
"""Synthetic, patient-free candidate-selection regression tests."""
import hashlib, json, pathlib, shutil, subprocess, tempfile, unittest
ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/a0_prepare_candidates.sh"
HEADER = """##fileformat=VCFv4.2
##contig=<ID=toy,length=1000>
##FILTER=<ID=LowQual,Description="Synthetic low quality">
##INFO=<ID=CSQ,Number=.,Type=String,Description="Synthetic annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Feature|HGVSc|HGVSp">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Phase genotype">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Phase identifier">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC
"""
ROWS = [
 (10,"C","C|missense_variant|MODERATE|GENE_A|TX1|c.1A>C|p.Ala1Val"),
 (10,"G","G|intron_variant|MODIFIER|GENE_A|TX1||"),
 (20,"T","T|intron_variant|MODIFIER|GENE_B|TX2||"),
 (30,"G","G|intron_variant|MODIFIER|HIGH|TX3||"),
 (40,"T","T|missense_variant|moderate|GENE_C|TX4||"),
 (50,"C","C|other|MODERATE_OTHER|GENE_D|TX5||"),
 (60,"G","G|other|HIGHER|GENE_E|TX6||"),
 (70,"C","C|intron_variant|MODIFIER|GENE_F|TX7||,C|stop_gained|HIGH|GENE_F|TX8|c.1A>C|p.Ala1Ter")]

def record(pos, alt, csq):
    return f"toy\t{pos}\t.\tA\t{alt}\t5\tLowQual\tCSQ={csq}\tGT:DP:GQ:AD:PGT:PID:PS\t0/1:2:1:1,1:0|1:synthetic:10\n"

def body(path):
    return subprocess.check_output(["bcftools","view","--no-header",str(path)], text=True).splitlines()

@unittest.skipUnless(shutil.which("bcftools") and shutil.which("awk"), "bcftools and awk required")
class CandidateTests(unittest.TestCase):
    def test_legacy_scope_and_annotation_retention(self):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td); src = root/'input.vcf'; dest = root/'out'
            src.write_text(HEADER + ''.join(record(*row) for row in ROWS))
            command = ['bash',str(SCRIPT),'--input',str(src),'--outdir',str(dest),'--progress-every','2']
            run = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            positions = {pos for pos, alt, csq in ROWS if any(
                field.upper() in {'HIGH','MODERATE'}
                for annotation in csq.split(',') for field in annotation.split('|'))}
            expected = [line for line in body(src) if int(line.split('\t')[1]) in positions]
            self.assertEqual(body(dest/'candidates.vcf.gz'), expected)
            self.assertEqual(len(expected), 5)
            report = json.loads((dest/'candidate_manifest.json').read_text())
            self.assertEqual(report['input_records'], 8)
            self.assertEqual(report['selected_input_records'], 4)
            self.assertEqual(report['kept_records'], 5)
            self.assertFalse(report['contig_filter_applied'])
            self.assertFalse(report['af_filter_applied'])
            self.assertEqual(report['hard_qc_filters_applied'], [])
            self.assertIn('scanned 8 records', run.stderr)
            original = hashlib.sha256((dest/'candidates.vcf.gz').read_bytes()).hexdigest()
            retry = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(retry.returncode, 0)
            self.assertEqual(hashlib.sha256((dest/'candidates.vcf.gz').read_bytes()).hexdigest(), original)

    def test_no_candidates_publishes_no_products(self):
        with tempfile.TemporaryDirectory() as td:
            root = pathlib.Path(td); src = root/'empty_selection.vcf'; dest = root/'out'
            src.write_text(HEADER + record(*ROWS[2]))
            run = subprocess.run(['bash',str(SCRIPT),'--input',str(src),'--outdir',str(dest)], capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn('No candidates selected', run.stderr)
            self.assertFalse((dest/'candidate_manifest.json').exists())
            self.assertFalse((dest/'candidates.vcf.gz').exists())

if __name__ == '__main__': unittest.main(verbosity=2)
