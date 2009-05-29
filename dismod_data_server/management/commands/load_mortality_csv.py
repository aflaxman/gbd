"""
load mortality csv
==================

Django command that loads the contents of command-line specified
mortality csv as rates for dismod 3.
"""

from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode
from django.core.management.base import BaseCommand, CommandError

import optparse
import simplejson as json

from gbd.dismod_data_server.models import Data, DiseaseModel

year_col = 1
age_col = 2
sex_col = 3
region_col = 4
rate_col = 5


def try_int(str):
    try:
        ret_val = int(str)
    except ValueError:
        ret_val = -99
    return ret_val

def try_float(str):
    try:
        ret_val = float(str)
    except ValueError:
        ret_val = -99
    return ret_val


class Command(BaseCommand):
    help = "Import all-cause mortality data from a .csv file."
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        import csv

        if len(fnames) != 1:
            raise CommandError('a single .csv file is required as input.')
        fname = fnames[0]

        print "adding rates from %s" % fname
        rate_counter = 0

        csv_file = csv.reader(open(fname))
        headings = csv_file.next()
        # check if the headings are as expected, complain if not
        if headings != ['', 'year', 'agegroup', 'sex', 'gbdregion', 'mortalityrate']:
            raise CommandError('csv file headings not as expected.')
        rates_for_region = {}

        for ii,row in enumerate(csv_file):
            params = {}
            params['data_type'] = 'all-cause mortality data'
            params['condition'] = 'all-cause mortality'

            # only load mortality data from 1990 and 2005
            if not (row[year_col] in ['1990', '2005']):
                continue

            params['year_start'] = int(row[year_col])
            params['year_end'] = int(row[year_col])

            if row[age_col] == '0_1':
                params['age_start'] = 0
                params['age_end'] = 0
            elif row[age_col] == '85above':
                params['age_start'] = 85
                params['age_end'] = -99
            else:
                ages = row[age_col].split('-')
                params['age_start'] = int(ages[0])
                params['age_end'] = int(ages[1])

            params['sex'] = row[sex_col]

            params['region'] = smart_unicode(row[region_col].strip(), errors='ignore')

            if row[rate_col] == 'NA':
                continue
            
            params['value'] = float(row[rate_col])
            params['standard_error'] = -99

            rate_data, created = Data.objects.get_or_create(**params)
            rate_counter += created

            key = str(params['region']) + params['sex'] + str(params['year_start'])
            if not rates_for_region.has_key(key):
                rates_for_region[key] = []
            rates_for_region[key] += [rate_data]

        print "added %d rates" % rate_counter

        for rates in rates_for_region.values():
            if len(rates) > 0:
                r = rates[0]
                dm, created = DiseaseModel.objects.get_or_create(condition=r.condition,
                                                                 region=r.region,
                                                                 sex=r.sex,
                                                                 year=r.year_start)
                dm.data = rates
                dm.save()
                print 'created: %s' % dm
