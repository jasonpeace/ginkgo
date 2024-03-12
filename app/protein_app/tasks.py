from celery import shared_task
from Bio import SeqIO, Align
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Align import AlignInfo
import os



@shared_task()
def find_filename(path):
    dir_list = os.listdir(path)
    print("Files and directories in '", path, "' :")
    print(dir_list)

@shared_task()
def find_first_match(query, path):
    query = query.upper()
    aligner = Align.PairwiseAligner(match_score=1.0, mode="local")
    aligner.open_gap_score = -5000.0
    file_list = os.listdir(path)
    for filename in file_list:
        print("Filename: {}".format(filename))
        gb_file = path + "/" + filename
        try:
            records = SeqIO.parse(open(gb_file, "r"), "genbank")
            for gb_record in records:
                for feature in gb_record.features:
                    if "protein_id" in feature.qualifiers:
                        target = feature.extract(gb_record.seq)
                        score = aligner.score(target, query)
                        alignments = aligner.align(target, query)
                        num_alignments = len(alignments)
                        if num_alignments == 1 and is_exact_match(alignments[0]):
                            print("Found a match!")
                            print("Protein ID: {} from Genome File: {}".format(feature.qualifiers['protein_id'], filename))
                            print("Alignment:{}".format(str(alignments[0])))
                            return(feature.qualifiers['protein_id'])
        except UndefinedSequenceError as e:
            print("File: {} is missing its origin sequence. Skipping.".format(filename))
            continue
        except:
            print("File: {} had an unexpected error. Skipping.".format(filename))
    

def is_exact_match(alignment):
    a = alignment[1, :]
    b = alignment[0, :]
    return a==b #dumb but weeds out alignments with mismatches
        



