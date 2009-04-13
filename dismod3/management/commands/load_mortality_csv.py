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

from dismod3.models import Rate, Region, Disease, Population

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
        disease, created = Disease.objects.get_or_create(name='all causes')
        
        rate_counter = 0

        csv_file = csv.reader(open(fname))
        headings = csv_file.next()
        # check if the headings are as expected, complain if not
        if headings != ['', 'year', 'agegroup', 'sex', 'gbdregion', 'mortalityrate']:
            raise CommandError('csv file headings not as expected.')

        for ii,row in enumerate(csv_file):
            params = {}
            params['rate_type'] = 'all-cause mortality data'
            params['disease'], created = Disease.objects.get_or_create(name='all-cause mortality')

            # only load mortality data from 1990 and 2005
            if not (row[1] in ['1990', '2005']):
                continue
            
            params['epoch_start'] = int(row[1])
            params['epoch_end'] = int(row[1])

            if row[2] == '0_1':
                params['age_start'] = 0
                params['age_end'] = 0
            elif row[2] == '85above':
                params['age_start'] = 85
                params['age_end'] = -99
            else:
                ages = row[2].split('-')
                params['age_start'] = int(ages[0])
                params['age_end'] = int(ages[1])

            params['sex'] = row[3]

            region_name = smart_unicode(row[4].strip(), errors='ignore')
            params['country'] = 'all in %s' % region_name
            params['region'], created = Region.objects.get_or_create(name=region_name)

            if row[5] == 'NA':
                continue
            
            params['numerator'] = float(row[5])*10.e7
            params['denominator'] = 10.e7

            rate_data, created = Rate.objects.get_or_create(**params)
            rate_counter += 1

        print "added %d rates" % rate_counter
