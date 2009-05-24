"""
script to generate the various fits for a disease model
"""

from django.core.management.base import BaseCommand

from dismod3.models import AgeSpecificRateFunction, DiseaseModel
from dismod3.bayesian_models import fit_disease_model

class Command(BaseCommand):
    help = 'generate the map fit, and mcmc fit for a disease model.'
    args = '[incidence_rate_id] [prevalence_rate_id] [remission_rate_id] [case_fatality_rate_id]'

    def handle(self, *id_list, **options):
        in_rfs = AgeSpecificRateFunction.objects.filter(id__in=id_list)

        rfs = []
        for rf in in_rfs:
            nrf = rf.clone(priors=rf.fit.get('priors', 'smooth 1.0'))
            rfs.append(nrf)
        
        rf = rfs[0]

        
        
        dm = DiseaseModel(disease=rf.disease, region=rf.region, sex=rf.sex)
        dm.save()
        dm.rates = rfs
        dm.save()

        print "Fitting %s\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())
        fit_disease_model.mcmc_fit(dm)

        print "\nModel %d is fit\n  http://winthrop.gs.washington.edu:5432%s" % (dm.id, dm.get_absolute_url())

