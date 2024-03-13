import os
import requests
from Bio import SeqIO, Align
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.Align import AlignInfo
from typing import Type


Type_PairwiseAligner = Type[Align.PairwiseAligner]
Type_Alignment = Type[Align.Alignment]


def find_first_match(query: str, path: str, alignment_request_id: int) -> dict:
    """
    This method uses the PairwiseAligner from Biopython to find a match on an
    alignment between a target and the query.

    query: A String that consists of the DNA we want to find in a protein
    path: A String that points to where the genbank genome files are located
    alignment_request_id: An Int that lets us look up the associated AlignmentRequest

    The aligner is set up to run in mode='local' with an open_gap_score set to -5000.0.
    The open_gap_score is set so low as to make it harder for the aligner to return alignments
    with gaps. Otherwise the number of alignments can be be quite large.

    """

    aligner = get_aligner("local", -5000.0)
    file_list = get_filenames(path)
    alignment_result = search_for_match(
        query, alignment_request_id, file_list, path, aligner
    )
    save_result(alignment_result)
    return alignment_result


def search_for_match(
    query: str,
    alignment_request_id: int,
    file_list: list[str],
    path: str,
    aligner: Type_PairwiseAligner,
) -> dict | None:
    for filename in file_list:
        gb_file = path + "/" + filename
        try:
            records = SeqIO.parse(open(gb_file, "r"), "genbank")
            for gb_record in records:
                for feature in gb_record.features:
                    target = str(feature.extract(gb_record.seq))
                    match = get_alignment_match(aligner, query, target, feature)
                    if match:
                        alignment_result = {
                            "protein_id": feature.qualifiers["protein_id"],
                            "alignment_detail": str(match),
                            "alignment_request_id": alignment_request_id,
                            "protein_dna_seq": target,
                            "organism": gb_record.annotations["organism"],
                            "filename": filename,
                        }
                        return alignment_result
        except UndefinedSequenceError as e:
            print("File: {} is missing its origin sequence. Skipping.".format(filename))
            continue
        except Exception as e:
            print(
                "File: {} had an unexpected error. Skipping. Error: {}".format(
                    filename, e
                )
            )


def get_filenames(path: str):
    """
    Grabs the filenames from the specified dir
    """
    file_list = os.listdir(path)
    return file_list


def get_aligner(mode: str, open_gap_score: float) -> Type_PairwiseAligner:
    """
    Sets up an aligner allowing for mode and open_gap_scote to be passed in.
    """
    aligner = Align.PairwiseAligner(match_score=1.0, mode=mode)
    aligner.open_gap_score = open_gap_score
    return aligner


def get_alignment_match(
    aligner: Type_PairwiseAligner, query: str, target: str, feature
) -> Type_Alignment | None:
    """
    Returns the match or None from an alignment between target and query if there is only only alignment
    and if there are no gaps or mismatches
    """
    if "protein_id" in feature.qualifiers:
        score = aligner.score(target, query)
        alignments = aligner.align(target, query)
        num_alignments = len(alignments)
        if num_alignments == 1 and is_exact_match(alignments[0]):
            return alignments[0]
    return None


def is_exact_match(alignment: Type_Alignment) -> bool:
    """
    Uses the counts of gaps and mismatches to determine if the alignment was a match.
    """
    counts = alignment.counts()
    return (
        counts.gaps == 0 and counts.mismatches == 0
    )  # should only let an exact match past


def save_result(data: dict) -> None:
    """
    Calls the detail api to save the match found in an alignment.
    alignment_request_id is required on the data passed in to build the correct url.
    """
    url = "http://app:80/api/detail/" + str(data["alignment_request_id"])
    x = requests.post(url, json=data)
