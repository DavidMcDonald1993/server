from django.db import models

# Create your models here.
class Compound(models.Model):

    name = models.CharField(max_length=200, default="")
    canonical_smiles = models.CharField(max_length=200, default="")
    added_date = models.DateTimeField("date_added")

    def __str__(self):
        return "{}-{}".format(
            self.name, 
            self.canonical_smiles)

# class PASSRequest(models.Model):
    # sdf_file = models.FileField

# class Question(models.Model):
#     question_text = models.CharField(max_length=200)
#     pub_date = models.DateTimeField('date published')

# class Choice(models.Model):
#     question = models.ForeignKey(Question, on_delete=models.CASCADE)
#     choice_text = models.CharField(max_length=200)
#     votes = models.IntegerField(default=0)