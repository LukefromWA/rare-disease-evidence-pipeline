import copy
import unittest
from d3_codon_genomic_mapping import map_record, verify_transcript

def fixture():
    ref = "A"*15 + "CAA" + "C"*15
    row = dict(gene="SYNTHETIC",candidate_id="1:1016:C:T",transcript_version="ENST00000000001.1",
        HGVSc="ENST00000000001.1:c.16C>T",HGVSp="p.Gln6Ter",alternate_codon_RNA="UAA",
        alternate_context_bracketed_RNA="A"*15+"[UAA]"+"C"*15,
        sequence_verified="True",alternate_is_stop="True",context_boundary_truncated="False",
        reference_assembly="GRCh38",upstream_bases="15",downstream_bases="15",
        reference_context_DNA=ref,alternate_context_DNA="A"*15+"TAA"+"C"*15)
    response={"mappings":[{"original":{"start":1,"end":33},"mapped":{
        "start":1001,"end":1033,"strand":1,"seq_region_name":"1","assembly_name":"GRCh38"}}]}
    return row,response,lambda chrom,start,end: ref

class MappingTests(unittest.TestCase):
    def test_full_33_base_window_and_bed(self):
        row,response,sequence=fixture()
        result=map_record(row,response,sequence)
        self.assertEqual(result["bed_chromStart"],"1015")
        self.assertEqual(result["bed_chromEnd"],"1018")
        self.assertEqual(result["genomic_reference_context_DNA"],row["reference_context_DNA"])
    def test_wrong_genomic_flank(self):
        row,response,_=fixture()
        with self.assertRaisesRegex(ValueError,"Full genomic"):
            map_record(row,response,lambda *args: "G"*33)
    def test_negative_strand(self):
        row,response,sequence=fixture(); response["mappings"][0]["mapped"]["strand"]=-1
        with self.assertRaisesRegex(ValueError,"Negative strand"):
            map_record(row,response,sequence)
    def test_candidate_position(self):
        row,response,sequence=fixture(); row["candidate_id"]="1:1017:C:T"
        with self.assertRaisesRegex(ValueError,"position mismatch"):
            map_record(row,response,sequence)
    def test_candidate_alleles(self):
        row,response,sequence=fixture(); row["candidate_id"]="1:1016:A:T"
        with self.assertRaisesRegex(ValueError,"alleles"):
            map_record(row,response,sequence)
    def test_truncation(self):
        row,response,sequence=fixture(); row["context_boundary_truncated"]="True"
        with self.assertRaisesRegex(ValueError,"Truncated"):
            map_record(row,response,sequence)
    def test_mapping_gap(self):
        row,response,sequence=fixture(); response["mappings"][0]["original"]["start"]=2
        with self.assertRaises(ValueError): map_record(row,response,sequence)
    def test_duplicate_mapping(self):
        row,response,sequence=fixture(); response["mappings"]*=2
        with self.assertRaisesRegex(ValueError,"Overlapping"):
            map_record(row,response,sequence)
    def test_nonstop_rejected(self):
        row,response,sequence=fixture(); row["alternate_is_stop"]="False"
        with self.assertRaises(ValueError): map_record(row,response,sequence)
    def test_wrong_transcript_version(self):
        row,_,_=fixture()
        meta=dict(id="ENST00000000001",version=2,species="homo_sapiens",object_type="Transcript",assembly_name="GRCh38")
        with self.assertRaisesRegex(ValueError,"version mismatch"): verify_transcript(meta,row)
    def test_exon_split_flank_supported(self):
        row,_,_=fixture(); ref=row["reference_context_DNA"]
        response={"mappings":[
            {"original":{"start":1,"end":9},"mapped":{"start":501,"end":509,"strand":1,"seq_region_name":"1","assembly_name":"GRCh38"}},
            {"original":{"start":10,"end":33},"mapped":{"start":1010,"end":1033,"strand":1,"seq_region_name":"1","assembly_name":"GRCh38"}}]}
        def sequence(chrom,start,end):
            return ref[:9] if start==500 else ref[9:]
        self.assertEqual(map_record(row,response,sequence)["bed_chromStart"],"1015")
    def test_exon_split_codon_rejected(self):
        row,_,_=fixture(); ref=row["reference_context_DNA"]
        response={"mappings":[
            {"original":{"start":1,"end":16},"mapped":{"start":1001,"end":1016,"strand":1,"seq_region_name":"1","assembly_name":"GRCh38"}},
            {"original":{"start":17,"end":33},"mapped":{"start":2017,"end":2033,"strand":1,"seq_region_name":"1","assembly_name":"GRCh38"}}]}
        with self.assertRaisesRegex(ValueError,"discontinuity"):
            map_record(row,response,lambda chrom,start,end: ref[:16] if start==1000 else ref[16:])

if __name__=="__main__": unittest.main()
