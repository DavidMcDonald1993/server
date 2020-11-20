from django.db import models

# Create your models here.
class NaturalProduct(models.Model):
    coconut_id = models.CharField(max_length=200)