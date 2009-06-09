"""
script to (re)generate the various fits for an age-specific rate function
"""

from django.core.management.base import BaseCommand, CommandError

from dismod3.models import AgeSpecificRateFunction
from dismod3.bayesian_models import fit_multiple_rates

class Command(BaseCommand):
    help = 'fit asrfs together, assuming a common underlying rate function'
    args = 'asrf_id_1 asrf_id_2 ...'

    def handle(self, *id_list, **options):
        
        if len(id_list) > 1:
            asrfs = AgeSpecificRateFunction.objects.filter(id__in=id_list)
        else:
            raise CommandError('At least 2 asrf ids are necessary.')

        #fit_multiple_rates.map_fit(asrfs)
        fit_multiple_rates.mcmc_fit(asrfs)
