from Bio import SeqIO, Align
from Bio.Seq import Seq

aligner = Align.PairwiseAligner(match_score=1.0)
aligner.mode = "local"


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


def align_test():
    gb_file = "../genome_genbank_files/NC_000852.gb"
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        parse_CDS_data(gb_record)

        # print the origin/full sequence out
        # print("Origin: {}".format(str(gb_record.seq)))

        test_seq = "ggtgggattgacagtggtaaatgtgttgac".upper()
        origin_seq = str(gb_record.seq)

        # test_seq = "AGAACTC"
        # origin_seq = "GAACT"
        try:
            score = aligner.score(test_seq, origin_seq)

            print("Score: {}".format(score))

            alignments = aligner.align(test_seq, origin_seq)

            alignment = alignments[0]
            indices = alignment.indices[0]
            print("Location: [{}..{}]".format(indices[0], indices[-1]))

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
        ]  # note Biopython's location indexes appear to not be zero based
        test_seq = test_seq.upper()
        # origin_seq = str(gb_record.seq)

        idx = gb_record.seq.find(test_seq)
        for feature in gb_record.features:
            
            if idx in feature:
                if "protein_id" in feature.qualifiers:
                    print("Starts in protein region:")
                    print(
                        "%s %s" % (feature.type, feature.qualifiers.get("protein_id"))
                    )
            
            start = feature.location.start -1
            end = feature.location.end -1
            if idx >= start and idx+len(test_seq) <= end:
                if "protein_id" in feature.qualifiers:
                        print("Starts and ends in protein region:")
                        print(
                            "%s %s" % (feature.type, feature.qualifiers.get("protein_id"))
                        )    


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
    for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
        origin_seq = str(gb_record.seq)
        test_seq = origin_seq[
            1408 - 1 : 1539 - 1
        ]  # note Biopython's location indexes appear to not be zero based
        test_rna = Seq(test_seq).transcribe()
        test_p = test_rna.translate()
        print("REAL: MGSAPLRPKSLHSATLRGLLRAPLRSATLRSASELRLPVVIVI")
        print(" GOT: {}".format(str(test_p)))


if __name__ == "__main__":
    snp_test()
