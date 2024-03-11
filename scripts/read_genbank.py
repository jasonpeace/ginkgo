from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
import sys
from Bio.Align import AlignInfo


def parse_CDS_data(gb_record):
    print("Name: {}".format(gb_record.name))
    for index, feature in enumerate(gb_record.features):
        if feature.type == "CDS":
            if "protein_id" in feature.qualifiers:
                location = feature.location
                product = feature.qualifiers["product"]
                protein_id = feature.qualifiers["protein_id"]
                print(
                    "Product: {} ID: {} Location Start:{} End: {} Direction: {}".format(
                        product,
                        protein_id,
                        location.start,
                        location.end,
                        location.strand,
                    )
                )


def align_on_protein():
    aligner = Align.PairwiseAligner(match_score=1.0, mode="global")
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        for feature in gb_record.features:
            protein_dna_seq = feature.extract(gb_record.seq)
            test_dna = "ggtgggattgacagtggtaaatgtgttgac".upper()

            # print("Protein DNA: {}".format(str(protein_dna_seq)))
            # print("Test DNA: {}".format(test_dna))

            score = aligner.score(str(protein_dna_seq), test_dna)
            print("Score: {}".format(score))
            alignments = aligner.align(str(protein_dna_seq), test_dna)
            # import pdb;
            # pdb.set_trace()





def score_test():
    
    dna1 = "ATATATATATATATATATGTATAT"
    dna2 = "AT"
    aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
    aligner.open_gap_score = -10.0
    score = aligner.score(dna1, dna2)
    print(score)
    alignments = aligner.align(dna1, dna2)
    print("Number of alignments: {}".format(len(alignments)))
    for alignment in alignments:
        #a = alignment[1, :]
        #if "-" not in a: #dumb but weeds out alignments with gaps
        print(alignment)


def align_test():
    aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        parse_CDS_data(gb_record)

        # print the origin/full sequence out
        # print("Origin: {}".format(str(gb_record.seq)))

        test_seq = "ggtgggattgacagtggtaaatgtgttgac".upper()
        origin_seq = str(gb_record.seq)

        # test_seq = "AGAACTC"
        # origin_seq = "ggtgggattgacagtggtaaatgtgttgac".upper()
        try:
            score = aligner.score(origin_seq, test_seq)

            print("Score: {}".format(score))

            alignments = aligner.align(origin_seq, test_seq)

            for a in alignments:
                print(a)
            # indices = alignment.indices[0]
            # print("Location: [{}..{}]".format(indices[0], indices[-1]))

            print(alignment)
        except:
            print("ORIGIN: {}".format(origin_seq))
            print("Test: {}".format(test_seq))

        # for alignment in alignments:
        #    print(alignment)


def snp_test():
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        # parse_CDS_data(gb_record)

        # print the origin/full sequence out
        # print("Origin: {}".format(str(gb_record.seq)))
        origin_seq = str(gb_record.seq)
        test_seq = origin_seq[
            1408 - 1 : 1539 - 1
        ]  # note genbank's location indexes appear to not be zero based

        test_seq = test_seq.upper()
        origin_seq = str(gb_record.seq)

        # idx = gb_record.seq.find(test_seq)
        # idxs = [idx]

        # find all indexes where the test sequence match, not just the first
        idxs = [i for i in range(len(origin_seq)) if origin_seq.startswith(test_seq, i)]
        for idx in idxs:
            for feature in gb_record.features:

                if idx in feature:
                    if "protein_id" in feature.qualifiers:
                        print("Starts in protein region:")
                        print(
                            "%s %s"
                            % (feature.type, feature.qualifiers.get("protein_id"))
                        )

                start = feature.location.start - 1
                end = feature.location.end - 1
                if idx >= start and idx + len(test_seq) <= end:
                    if "protein_id" in feature.qualifiers:
                        print("Starts AND ends in protein region:")
                        print(
                            "%s %s"
                            % (feature.type, feature.qualifiers.get("protein_id"))
                        )

def score_hit_test():
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        origin_seq = str(gb_record.seq)
        test_seq = origin_seq[
            1408 - 1 : 1539 - 1
        ]
        test_seq = test_seq.upper()
        origin_seq = str(gb_record.seq)
        max_score = len(test_seq)
        aligner = Align.PairwiseAligner(match_score=max_score, mode="local")

            
    dna1 = "ATATATATATATATATATGTATAT"
    dna2 = "TGT"
    score = aligner.score(dna1, dna2)
    print(score)



def location_test():
    #   CDS             1408..1539
    #                  /gene="a002dR"
    #                  /locus_tag="PBCV1_a002dR"
    #                  /codon_start=1
    #                  /product="hypothetical protein"
    #                  /protein_id="YP_004678873.1"
    #                  /db_xref="GeneID:10971237"
    #                  /translation="MGSAPLRPKSLHSATLRGLLRAPLRSATLRSASELRLPVVIVI"
    gb_file = "../genome_genbank_files/NC_000852.gb"
    translated = "MGSAPLRPKSLHSATLRGLLRAPLRSATLRSASELRLPVVIVI"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        origin_seq = str(gb_record.seq)
        test_seq = origin_seq[
            1408 - 1 : 1539 - 1
        ]  # note genbank's location indexes appear to not be zero based
        test_rna = Seq(test_seq).transcribe()
        test_p = test_rna.translate()
        print("REAL: {}".format(translated))
        print("slicing:")
        print(" GOT: {}".format(str(test_p)))
        (print("MATCH!") if test_p == translated else print("NO match"))

        print("Using SimpleLocation")
        feature = SeqFeature(
            SimpleLocation(1408 - 1, 1539 - 1, strand=1), type="CDS"
        )  # note still need to -1 the indexes, shouldn't be necessary if you use location.start location.end
        test_seq2 = feature.extract(gb_record.seq)
        test_rna2 = test_seq2.transcribe()
        test_p2 = test_rna2.translate()
        print(" GOT: {}".format(str(test_p2)))
        (print("MATCH!") if test_p2 == translated else print("NO match"))


def aligntest():
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    #aligner.extend_gap_score = -0.5
    # aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        for feature in gb_record.features:
            if "protein_id" in feature.qualifiers:
                protein_dna_seq = feature.extract(gb_record.seq)
                test_dna = str(protein_dna_seq)[
                    len(str(protein_dna_seq)) // 2 : len(str(protein_dna_seq)) 
                ]  # half the dna seq
                score = aligner.score(str(protein_dna_seq), test_dna)
                print("Score: {}".format(score))
                alignments = aligner.align(str(protein_dna_seq), test_dna)
                print("Num alignments: {}".format(len(alignments)))
            
                # no_gap_alignment = find_first_no_gap(alignments)
                # if no_gap_alignment:

                # print("--------------------------------------------------------------------------------")
                # print(protein_dna_seq)
                # print(test_dna)
                # print(alignments[0])
                # print("--------------------------------------------------------------------------------")

                
  
def find_first_no_gap(alignments):
    for idx, alignment in enumerate(alignments):
        a = alignment[1, :]
        if "-" not in a: #dumb but weeds out alignments with gaps
            print("INDEX:{}".format(idx))
            return alignment


if __name__ == "__main__":
    score_test()

