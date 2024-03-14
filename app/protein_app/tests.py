import unittest
from django.test import TestCase
from protein_app.modules.alignment import is_exact_match, get_aligner, search_for_match
from .models import AlignmentRequest, AlignmentResult

# just adding tests here for now.

DNA_QUERY_1 = "ATGGTTCACAATTGCAAAAGTATATTCACTGCAGCAGAAAATGGCCACGACGTTTGTTTGAAAACGCTCATTGAAGCAGGTGCCCCCTTTGACAACGTCGGTGATTCAGGGTGGACCGCGTTGCATTACGCGATTCATTATGATCATACTGCGTGCGTAAAGATGCTCATTGATGCGGGTGCAAATATTGACATCACAGATAATTCGGGATGCACACCACTTCATCGTGCGGTTTTTAATGGCCATGATGCATGTGTGAAGCTGCTCGTCGAAGCAGGTGCAACTCTTGACGTCATTGATGATACTGGAACCATGCCACTGCATCACGCAGTTTATTATGGTTATGATGCATGCGTAAAGATGCTCATAGAGGCAGGTGCCGGTCTTAACATCGACGGTGATGGGTATGCACCGTTACATTACGCGGTTTATAAAGGTCACGATGTGTGTGTGAAGCTGCTCGTCGAAGCCGGTGCACCCCTTGACATCACAGATATTTCGGGATGCACACCACTTCATCGTGCGGTTTTTAATGGCCACGATGCATGTGCGAGCATGTTAGTCAACAAGATCGTTTCGGAGCGGCCGTTGCGTCCGAGTGAGTTGTGTGTCATACCACAAACATCTGCTGTATTAGGTGATGTGTTGCGAACGACGATGCAGCTTCATGGGCGATCGGAAGCTGCAAAGATCACAGCGCATCTTCCTGTGGGCGCAAGGGATACTCTGCGGACTACTATGCTGTGTTTGAACAGGACCATGGTCCCGAGAGACCTCATTGACAGCATAGTACTCCAATGTGTGTA"
ALIGNMENT_RESULT =  {
    "protein_id": ['YP_004678874.1'],
    "alignment_detail": "target           67 ATGGTTCACAATTGCAAAAGTATATTCACTGCAGCAGAAAATGGCCACGACGTTTGTTTG\n                  0 ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nquery             0 ATGGTTCACAATTGCAAAAGTATATTCACTGCAGCAGAAAATGGCCACGACGTTTGTTTG\n\ntarget          127 AAAACGCTCATTGA 141\n                 60 ||||||||||||||  74\nquery            60 AAAACGCTCATTGA  74\n",
    "protein_dna_seq": "ATGAGAAAATACTATAAAACCATACAAACGTCATACAAACCACCATACAAACACCTTAAAGTATACAATGGTTCACAATTGCAAAAGTATATTCACTGCAGCAGAAAATGGCCACGACGTTTGTTTGAAAACGCTCATTGA",
    "organism": "Paramecium bursaria Chlorella virus 1",
    "filename": "NC_000852.gb",
}

class AlignmentTestCase(unittest.TestCase):
    def test_is_exact_match(self):
        query = DNA_QUERY_1
        aligner = get_aligner("local", -5000.0)
        path = "./genome_genbank_files"
        file_list = ["NC_000852.gb"]
        result = search_for_match(query, 8, file_list, path, aligner)
        self.assertEqual(
            result["organism"], ALIGNMENT_RESULT["organism"]
        )

        self.assertEqual(
            result["protein_id"], ALIGNMENT_RESULT["protein_id"]
        )

        self.assertEqual(
            result["protein_dna_seq"], ALIGNMENT_RESULT["protein_dna_seq"]
        )



class ModelTestCase(TestCase):

    def test_alignment_request_creation(self):
        alignment_request = AlignmentRequest(dna_string=DNA_QUERY_1, status="NEW")
        self.assertEqual(alignment_request.dna_string, DNA_QUERY_1 ) 

    def test_alignment_result_creation(self):
        alignment_request = AlignmentRequest(dna_string=DNA_QUERY_1, status="NEW")
        alignment_result = AlignmentResult(
            alignment_request = alignment_request,
            protein_id = ALIGNMENT_RESULT["protein_id"],
            alignment_detail = ALIGNMENT_RESULT["alignment_detail"],
            protein_dna_seq = ALIGNMENT_RESULT["protein_dna_seq"],
            organism =ALIGNMENT_RESULT["organism"],
            filename = ALIGNMENT_RESULT["filename"],
        )
        self.assertEqual(alignment_result.alignment_request, alignment_request)
        self.assertEqual(alignment_result.protein_id, ALIGNMENT_RESULT["protein_id"])




    




