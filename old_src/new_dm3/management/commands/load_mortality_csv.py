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

from new_dm3.models import Data, DiseaseModel
from dismod3.models.fields import MISSING

def try_int(str):
    try:
        ret_val = int(str)
    except ValueError:
        ret_val = MISSING
    return ret_val

def try_float(str):
    try:
        ret_val = float(str)
    except ValueError:
        ret_val = MISSING
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
            # only load mortality data from specific years
            if not (row[1] in ['2000']):
                continue
            
            opts = {}

            opts['params_json'] = json.dumps({'units': 'per one', 'source': fname})

            opts['condition'] = '(all)'
            opts['data_type'] = 'all-cause mortality data'
                
            region_name = smart_unicode(row[4].strip(), errors='ignore')
            opts['region'] = region_name
            opts['gbd_region'] = region_name

            opts['year_start'] = int(row[1])
            opts['year_end'] = int(row[1])

            opts['sex'] = row[3]
            
            if row[2] == '0_1':
                opts['age_start'] = 0
                opts['age_end'] = 0
            elif row[2] == '85above':
                opts['age_start'] = 85
                opts['age_end'] = MISSING
            else:
                ages = row[2].split('-')
                opts['age_start'] = int(ages[0])
                opts['age_end'] = int(ages[1])


            if row[5] == 'NA':
                continue

            opts['value'] = float(row[5])
            opts['standard_error'] = 0.

            d, is_new = Data.objects.get_or_create(**opts)
            d.cache_params()
            d.save()

            if is_new:
                rate_counter += 1

            key = str(opts['region']) + opts['sex'] + str(opts['year_start'])
            if not rates_for_region.has_key(key):
                rates_for_region[key] = []
            rates_for_region[key] += [d]

        print "added %d rates" % rate_counter


        # collect data together into models
        for rates in rates_for_region.values():
            if len(rates) > 0:

                d = rates[0]
                
                args = {}
                args['condition'] = d.condition
                args['sex'] = d.sex
                args['region'] = d.region
                args['year'] = str(d.year_start)

                dm = DiseaseModel.objects.create(**args)

                dm.data = rates
                dm.cache_params()
                dm.save()

