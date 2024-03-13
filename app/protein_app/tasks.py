from celery import shared_task
from Bio import SeqIO, Align
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Align import AlignInfo
import os
import requests


@shared_task()
def find_filename(path):
    dir_list = os.listdir(path)
    print("Files and directories in '", path, "' :")
    print(dir_list)


@shared_task()
def find_first_match(query, path, alignment_request_id):
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
                        target = str(feature.extract(gb_record.seq))
                        score = aligner.score(target, query)
                        alignments = aligner.align(target, query)
                        num_alignments = len(alignments)
                        if num_alignments == 1 and is_exact_match(alignments[0]):

                            alignment_result = {
                                "protein_id": feature.qualifiers["protein_id"],
                                "alignment_detail": str(alignments[0]),
                                "alignment_request_id": alignment_request_id,
                                "protein_dna_seq": target,
                                "organism": gb_record.annotations["organism"],
                                "filename": filename,
                            }
                            print(alignment_result)
                            save_result(alignment_result)
                            return alignment_result

        except UndefinedSequenceError as e:
            print(e)
            print("File: {} is missing its origin sequence. Skipping.".format(filename))
            continue
        except Exception as e:
            print(e)
            print("File: {} had an unexpected error. Skipping.".format(filename))


def is_exact_match(alignment):
    counts = alignment.counts()
    return counts.gaps == 0 and counts.mismatches == 0 #should only let an exact match past


def save_result(data):
    url = "http://app:8000/api/detail/" + str(data["alignment_request_id"])
    x = requests.post(url, json=data)
