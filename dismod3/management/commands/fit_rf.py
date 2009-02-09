"""
script to (re)generate the various fits for an age-specific rate function
"""

from django.core.management.base import BaseCommand

class Command(BaseCommand):
    help = '(re)generate the normal approximation,map fit, and mcmc fit for an age-specific rate function.'
    args = '[age_specific_rrate_ffunction_id_1 [asrf_id_2 ...]]'

    def handle(self, *id_list, **options):
        from dismod3.models import AgeSpecificRateFunction
        from dismod3.models.probabilistic_utils import map_fit, mcmc_fit
        
        if id_list:
            asrfs = AgeSpecificRateFunction.objects.filter(id__in=id_list)
        else:
            asrfs = AgeSpecificRateFunction.objects.all()

        if asrfs.count() > 5:
            for asrf in asrfs:
                print "\n\nFast Fitting %s" % asrf.pk
                try:
                    map_fit(asrf)
                except ValueError:
                    print "failed"
        else:
            for asrf in asrfs:
                print "\n\nFast Fitting %s" % asrf.pk
                mcmc_fit(asrf, speed='fast')

        for asrf in asrfs:
            print "\n\nFitting %s" % asrf.id
            try:
                mcmc_fit(asrf)
            except KeyError:
                print "failed"
            
