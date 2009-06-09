from django.db import models

class Disease(models.Model):
    """Model for a particular disease"""
    name = models.CharField(max_length=200)
    params_json = models.TextField(blank=True)
    def __unicode__(self):
        return self.name

    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3' 
