"""
load mortality csv
==================

Django command that loads the contents of command-line specified
mortality csv as rates for dismod 3.
"""

from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode
from django.core.management.base import BaseCommand, CommandError

import sys
import optparse
import csv
import simplejson as json

from gbd.dismod_data_server.models import Data, DiseaseModel, User
from gbd.dismod3.utils import MISSING, MAX_AGE

# collect all the column number information here,
# in case we decide to change it in the future.
"""
col_headings = ['', 'year', 'agegroup', 'sex', 'gbdregion', 'mortalityrate']
year_col = 1
age_col = 2
sex_col = 3
region_col = 4
rate_col = 5
"""
col_headings = ['gbd_region', 'sex', 'year', 'age', 'mx', 'agegroup']
year_col = 2
age_col = 5
sex_col = 1
region_col = 0
rate_col = 4

class Command(BaseCommand):
    help = 'Import all-cause mortality data from a .csv file.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        if len(fnames) != 1:
            raise CommandError('a single .csv file is required as input.')
        fname = fnames[0]

        data = [d for d in csv.DictReader(open(fname))]

        # check if the headings are as expected, complain if not
        if set(data[0].keys()) != set(col_headings):
            raise CommandError('csv file headings not as expected.')

        print 'adding rates from %s' % fname
        rate_counter = 0
        rates_for_region = {}

        for ii, d in enumerate(data):
            params = {}
            params['data_type'] = 'all-cause mortality data'
            params['condition'] = 'all-cause_mortality'

            # only load mortality data from 1990 and 2005
            if not (d['year'] in ['1990', '2005']):
                continue

            # all mortality data is specific to a single year
            params['year_start'] = int(d['year'])
            params['year_end'] = int(d['year'])

            # age ranges have several cases
            if d['agegroup'] == '0_1':
                params['age_start'] = 0
                params['age_end'] = 0
            #elif row[age_col] == '85above':
            elif row[age_col] == '80 plus':
                #params['age_start'] = 85
                params['age_start'] = 80
                params['age_end'] = MAX_AGE
            else:
                #ages = row[age_col].split('-')
                ages = row[age_col].split('_')
            elif d['agegroup'] == '80 plus':
                params['age_start'] = 80
                params['age_end'] = MAX_AGE
            else:
                ages = d['agegroup'].split('_')
                params['age_start'] = int(ages[0])
                params['age_end'] = int(ages[1])

            # sex field may have capitalization
            params['sex'] = d['sex'].lower()

            # gbd regions don't have non-ascii chars, but maybe some
            # future region will
            params['region'] = smart_unicode(d['gbd_region'].strip(), errors='ignore')
            params['gbd_region'] = params['region']

            # skip rows without rates
            if d['mx'] == 'NA':
                continue
            
            params['value'] = float(d['mx'])
            params['standard_error'] = MISSING  # no standard error is given for mortality data

            rate_data, created = Data.objects.get_or_create(**params)

            # cache data values in params_json
            rate_data.cache_params()
            rate_data.save()
            
            rate_counter += created

            key = str(params['region']) + params['sex'] + str(params['year_start'])
            if not rates_for_region.has_key(key):
                rates_for_region[key] = []
            rates_for_region[key] += [rate_data]

            if ii % 100 == 0:
                print '.',
                sys.stdout.flush()

        print 'added %d rates' % rate_counter

        all_rates = []
        for rates in rates_for_region.values():
            all_rates += rates
            
        dm = DiseaseModel(condition='all-cause_mortality',
                          region='World',
                          sex='all',
                          year='1990-2005',
                          creator=User.objects.get(id=1))
        #dm.cache_params()
        dm.save()
        for r in all_rates:
            dm.data.add(r)
        dm.save()
        print 'created: %s' % dm
        print 'url: %s' % dm.get_absolute_url()

