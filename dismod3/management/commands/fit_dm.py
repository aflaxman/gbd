"""
script to generate the various fits for a disease model
"""

from django.core.management.base import BaseCommand

from dismod3.models import DiseaseModel
from dismod3.bayesian_models import fit_disease_model

class Command(BaseCommand):
    help = 'generate the map fit, and mcmc fit for a disease model.'
    args = '[disease_model_id_1 [dm_id_2 ...]]'

    def handle(self, *id_list, **options):
        dmods = DiseaseModel.objects.filter(id__in=id_list)

        for dm in dmods:
            print "\n\nFitting %s" % dm.id
            fit_disease_model.mcmc_fit(dm)
            
