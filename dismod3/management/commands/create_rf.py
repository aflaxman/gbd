from django.core.management.base import BaseCommand, CommandError

from dismod3.tests.bayesian_probability_test import create_test_asrf

RATE_TYPE = { 'i': 'incidence data',
              'p': 'prevalence data',
              'r': 'remission data',
              'cf': 'case fatality data',
              }

class Command(BaseCommand):
    help = 'Generate rate function data for testing purposes.\n  rate_function example: (age/100.)**2/3.'
    args = 'rate_type=(i|p|r|cf) rate_function'
                  

    def handle(self, *args, **options):
        rate_type_str, rate_function_str = args
        rf = create_test_asrf(rate_function_str,
                              rate_type=RATE_TYPE[rate_type_str])
        print "ASRF %d created:\n  http://winthrop.gs.washington.edu:5432%s" % (rf.id, rf.get_absolute_url())
        try:
            pass
        except:
            raise CommandError('invalid input arguments')


