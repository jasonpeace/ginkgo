from django.db import models

class AlignmentRequest(models.Model):
    dna_string = models.TextField(null=False)
    status = models.CharField()
    date_submitted = models.DateTimeField(auto_now_add=True)
    date_updated = models.DateTimeField(auto_now=True)


class AlignmentResult(models.Model):
    alignment_request = models.ForeignKey(AlignmentRequest, on_delete=models.CASCADE)
    protein_id = models.CharField()
    alignment_detail = models.TextField()
    protein_dna_seq = models.TextField()
    alignment_index = models.IntegerField()
    genome_name = models.CharField()
    filename = models.CharField()
