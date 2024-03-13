import json
from django.core import serializers
from django.shortcuts import render, HttpResponse
from django.http import HttpRequest
from django.views.decorators.csrf import csrf_exempt
from .utils import (
    alignment_request_dto,
    alignment_requests_dto,
    alignment_result_dto,
)
from protein_app.models import AlignmentRequest, AlignmentResult
from .tasks import alignment_match


def home(request, id=0)->HttpResponse:
    return render(request, "home.html")


# todo work on these csrf stuff
@csrf_exempt
def get_requests(request: HttpRequest)->HttpResponse:
    """
    Returns all of the AlignmentRequests sorted so the newest are first.
    No pagination or filtering is implemented at this time.

    """
    alignment_requests = AlignmentRequest.objects.all().order_by("-id")
    ar = json.dumps(alignment_requests_dto(alignment_requests))
    return HttpResponse(ar, content_type="application/json")


@csrf_exempt
def new_request(request: HttpRequest)->HttpResponse:
    """
    This will take a new DNA query and persist it with the status being set to "NEW".
    It also kicks off a request to celery to attempt to search for a match to the
    DNA query in the exiting genome files on hand.

    See app/genome_genbank_files for the genomes that will be searched.

    This is expecting a POST request.

    """
    if request.method == "POST":
        data = json.loads(request.body)
        dna_string = data["dna_string"]
        alignment_request = AlignmentRequest(dna_string=dna_string, status="new")
        alignment_request.save()
        alignment_match.delay(
            dna_string, "/user/src/app/genome_genbank_files", alignment_request.id
        )
        ar = json.dumps(alignment_request_dto(alignment_request))
        return HttpResponse(ar, content_type="application/json")
    return HttpResponse("nothing to do")


@csrf_exempt
def alignment_detail(request: HttpRequest, alignment_request_id: int = -1) -> HttpResponse:
    """
    This endpoint will accept a GET or a POST.

    This GET will return the AlignmentRequest information as well as
    the AlignmentResult information if it exists.

    The POST will take a dict containing the info to find an existing
    AlignmentRequest object and create a new AlignmentResult object. It
    will set the AlignmentRequest's status to 'COMPLETED' and saves the
    associated AlignmentResult.

    The data dict passed in for the AlignmentResult should have these
    required props:

    alignment_request_id, protein_id,

    The following are optional:
    alignment_detail, protein_dna_seq, organism, filename
    """
    if request.method == "GET":
        alignment_request = AlignmentRequest.objects.get(pk=alignment_request_id)
        alignment_results = AlignmentResult.objects.filter(
            alignment_request=alignment_request
        )
        data = {
            "alignment_request": alignment_request_dto(alignment_request),
        }
        if len(alignment_results) > 0:
            data["alignment_results"] = alignment_result_dto(alignment_results[0])
        return HttpResponse(json.dumps(data), content_type="application/json")

    if request.method == "POST":
        data = json.loads(request.body)
        alignment_request = AlignmentRequest.objects.get(
            pk=data["alignment_request_id"]
        )
        if data:
            alignment_result = AlignmentResult(
                protein_id=data["protein_id"],
                alignment_detail=data.get("alignment_detail", "No alignment details"),
                alignment_request=alignment_request,
                protein_dna_seq=data.get("protein_dna_seq"),
                organism=data.get("organism"),
                filename=data.get("filename"),
            )
            alignment_result.save()
        alignment_request.status = "COMPLETE"
        alignment_request.save()
        ar = json.dumps(alignment_request_dto(alignment_request))
        return HttpResponse(ar, content_type="application/json")
