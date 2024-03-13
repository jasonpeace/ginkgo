from django.shortcuts import render, HttpResponse
from protein_app.models import AlignmentRequest, AlignmentResult
from django.core import serializers
from django.views.decorators.csrf import csrf_exempt
import json
from .tasks import find_first_match


def home(request, id=0):
    alignment_requests = AlignmentRequest.objects.all()
    return render(request, "home.html", {"alignment_requests": alignment_requests})


# todo work on these csrf stuff
@csrf_exempt
def get_requests(request):
    alignment_requests = AlignmentRequest.objects.all()
    ar = json.dumps(alignment_requests_dto(alignment_requests))
    return HttpResponse(ar, content_type="application/json")


# todo move to a utility module
def alignment_request_dto(alignment_request):
    return {
        "id": alignment_request.id,
        "dna_string": alignment_request.dna_string,
        "status": alignment_request.status,
        "date_submitted": alignment_request.date_submitted.strftime(
            "%m/%d/%Y, %H:%M:%S"
        ),
        "date_updated": alignment_request.date_updated.strftime("%m/%d/%Y, %H:%M:%S"),
    }


# todo move to a utility module
def alignment_result_dto(alignment_result):
    return {
        "id": alignment_result.id,
        "protein_id": alignment_result.protein_id,
        "alignment_detail": alignment_result.alignment_detail,
        "protein_dna_seq": alignment_result.protein_dna_seq,
        "organism": alignment_result.organism,
        "filename": alignment_result.filename,
    }


# todo move to a utility module
def alignment_requests_dto(alignment_requests):
    r = []
    for ar in alignment_requests:
        r.append(alignment_request_dto(ar))
    return r


@csrf_exempt
def new_request(request):
    if request.method == "POST":
        data = json.loads(request.body)
        dna_string = data["dna_string"]
        alignment_request = AlignmentRequest(dna_string=dna_string, status="new")
        alignment_request.save()
        find_first_match.delay(
            dna_string, "/user/src/app/genome_genbank_files", alignment_request.id
        )
        ar = json.dumps(alignment_request_dto(alignment_request))
        return HttpResponse(ar, content_type="application/json")
    return HttpResponse("nothing to do")


@csrf_exempt
def status(request, data):
    alignment_request = AlignmentRequest.objects.get(pk=alignment_request_id)
    if request.method == "POST":
        data = json.loads(request.body)
        alignment_request.status = data["status"]
        alignment_request.save()
    status = serializers.serialize("json", alignment_request.status)
    return HttpResponse(status, content_type="application/json")


@csrf_exempt
def alignment_detail(request, alignment_request_id=-1):
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
                protein_dna_seq = data.get("protein_dna_seq"),
                organism = data.get("organism"),
                filename = data.get("filename"),
            )
            alignment_result.save()
        alignment_request.status = "COMPLETE"
        alignment_request.save()
        ar = json.dumps(alignment_request_dto(alignment_request))
        return HttpResponse(ar, content_type="application/json")
