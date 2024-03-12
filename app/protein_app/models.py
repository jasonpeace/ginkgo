from django.db import models

class AlignmentRequest(models.Model):
    dna_string = models.TextField(null=False)
    status = models.CharField(max_length=100)
    date_submitted = models.DateTimeField(auto_now_add=True)
    date_updated = models.DateTimeField(auto_now=True)


class AlignmentResult(models.Model):
    alignment_request = models.ForeignKey(AlignmentRequest, on_delete=models.CASCADE)
    protein_id = models.CharField(max_length=100)
    alignment_detail = models.TextField()
    protein_dna_seq = models.TextField()
    alignment_index = models.IntegerField()
    genome_name = models.CharField(max_length=100)
    filename = models.CharField(max_length=200)