"""
Command to load covariate data from a csv

Use from the project directory as follows::

    $ python2.5 manage.py load_covariate_csv [options] filename.csv

"""

from django.core.management.base import BaseCommand, CommandError

import csv
import re

from gbd.covariate_data_server.models import *
from gbd.dismod3.settings import gbd_regions
from gbd.dismod3.utils import clean

class Command(BaseCommand):
    help = 'Import covariate data from a .csv file.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        if len(fnames) != 1:
            raise CommandError('a single .csv file is required as input.')
        fname = fnames[0]

        print "adding population data from %s" % fname

        csv_f = csv.DictReader(open(fname))
        data = [d for d in csv_f]
        headings = csv_f.fieldnames

        type_slug = 'GDPpc'
        type_desc = 'GDP per capita'
        type, is_new = CovariateType.objects.get_or_create(slug=type_slug, defaults={'description': type_desc})

        added = 0
        modified = 0

        for d in data:
            iso3 = d['iso3']
            sex = 'total'
            for key in d.keys():
                for year in re.findall('\d+', key):
                    year = int(year)
                    try:
                        value = float(d[key])
                    except ValueError:
                        continue
                    if year > 1900 and year < 2050:
                        cov, is_new = Covariate.objects.get_or_create(type=type, iso3=iso3, year=year, sex=sex, defaults={'value': value})
                        cov.value = value
                        cov.save()
                        added += is_new
                        modified += not is_new

        print 'added %d country-years of covariate data' % added
        print 'modified %d country-years of covariate data' % modified
        
