"""
script to (re)generate the various fits for an age-specific rate function
"""

from django.core.management.base import BaseCommand

from dismod3.models import AgeSpecificRateFunction
from dismod3.bayesian_models import fit_rate_function

class Command(BaseCommand):
    help = '(re)generate the normal approximation,map fit, and mcmc fit for an age-specific rate function.'
    args = '[age_specific_rate_function_id_1 [asrf_id_2 ...]]'

    def handle(self, *id_list, **options):
        asrfs = AgeSpecificRateFunction.objects.filter(id__in=id_list)

        for asrf in asrfs:
            print "\n\nFast Fitting %s" % asrf.pk
            fit_rate_function.mcmc_fit(asrf, speed='fast')

        for asrf in asrfs:
            print "\n\nFitting %s" % asrf.id
            fit_rate_function.mcmc_fit(asrf)
            
