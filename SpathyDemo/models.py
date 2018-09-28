from django.db import models

# Create your models here.
class SpathyFiguresModel(models.Model):
    figure_1 = models.CharField(default='name',max_length=255)
    figure_2 = models.CharField(default='name',max_length=255)
    
