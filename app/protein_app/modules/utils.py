def alignment_request_dto(alignment_request:dict):
    """
    Data Transfer Object for an AlignmentRequest
    """
    return {
        "id": alignment_request.id,
        "dna_string": alignment_request.dna_string,
        "status": alignment_request.status,
        "date_submitted": alignment_request.date_submitted.strftime(
            "%m/%d/%Y, %H:%M:%S"
        ),
        "date_updated": alignment_request.date_updated.strftime("%m/%d/%Y, %H:%M:%S"),
    }


def alignment_result_dto(alignment_result:dict):
    """
    Data Transfer Object for an AlignmentResult
    """
    return {
        "id": alignment_result.id,
        "protein_id": alignment_result.protein_id,
        "alignment_detail": alignment_result.alignment_detail,
        "protein_dna_seq": alignment_result.protein_dna_seq,
        "organism": alignment_result.organism,
        "filename": alignment_result.filename,
    }


def alignment_requests_dto(alignment_requests:dict):
    """
    Data Transfer Object for a list of AlignmentRequests
    """
    r = []
    for ar in alignment_requests:
        r.append(alignment_request_dto(ar))
    return r
