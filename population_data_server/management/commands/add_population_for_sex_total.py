"""
Command to make a sex=total population entry out of each sex=male/female

Use from the project directory as follows::

    $ python2.5 manage.py add_population_for_sex_total
"""

from django.core.management.base import BaseCommand
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import numpy as np
import simplejson as json

from gbd.population_data_server.models import Population

class Command(BaseCommand):
    help = 'make a sex=total population entry out of each sex=male/female'

    def handle(self, *args, **options):
        for male_pop in Population.objects.filter(sex='male'):
            female_pop = Population.objects.filter(region=male_pop.region, year=male_pop.year, sex='female')[0]

            total_params = male_pop.params
            total_params['vals'] = list(np.array(male_pop.params['vals']) + np.array(female_pop.params['vals']))

            total_pop, is_new = Population.objects.get_or_create(region=male_pop.region, year=male_pop.year, sex='total')
            total_pop.params_json = json.dumps(total_params)
            total_pop.save()
            print total_pop

