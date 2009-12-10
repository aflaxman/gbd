"""
Command to create a csv of population data

Use from the project directory as follows::

    $ python2.5 manage.py create_population_csv
"""

from django.core.management.base import BaseCommand, CommandError
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import csv

from gbd.population_data_server.models import Population
from gbd.dismod3.utils import MAX_AGE

class Command(BaseCommand):
    help = 'Create .csv file from population data.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        if len(fnames) == 0:
            fname = 'population.csv'
        elif len(fnames) == 1:
            fname = fnames[0]
        else:
            raise CommandError('a single .csv file is required as input.')

        f = open(fname, 'w')
        fields = ['Country Code', 'Year', 'Sex'] + ['Age %d Population' % a for a in range(MAX_AGE)]
        csv_f = csv.DictWriter(f, fields)
        csv_f.writer.writerow(fields)

        for j, pop in enumerate(Population.objects.filter(year__in=[1990, 2005], sex__in=['male', 'female'])):
            M,C = pop.gaussian_process()
            d = {'Country Code': pop.region,
                 'Year': pop.year,
                 'Sex': pop.sex}
            for i,p in enumerate(M(range(MAX_AGE))):
                d['Age %d Population'%i] = p
                
            csv_f.writerow(d)

            print '%d: %s' % (j, pop)
        f.close()
