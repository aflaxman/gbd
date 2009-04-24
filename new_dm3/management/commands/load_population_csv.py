"""
Command to load the GBD population data from a csv

Use from the project directory as follows
$ python2.5 manage.py load_population_csv USABLE_IHME_GBD_POPULATION_1950-2050.csv

The file USABLE_IHME_GBD_POPULATION_1950-2050.csv is big,
so it takes 2 minutes to load.
"""

from django.core.management.base import BaseCommand
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import simplejson as json
import re

from new_dm3.models import Data
from dismod3.models.fields import MISSING

class Command(BaseCommand):
    """
    load the global population tables from a csv file formatted to
    match the gbd format, which has header:
    gbd_country, iso3, gbd_region, year, sex, \
    pop_0to1, pop_1to4, pop_5to9, pop_10to14, ...
    """
    help = "Import population data from a .csv file."
    args = '[ filename.csv [ additional_filenames.csv ] ]'

    def handle(self, *fnames, **options):
        import csv
        import numpy as np
        
        if not fnames:
            fnames = ['../USABLE_IHME_GBD_POPULATION_1950-2050.csv']

        for fname in fnames:
            print "adding population data from %s" % fname
            pop_counter = 0
            
            csv_file = csv.reader(open(fname))
            headings = csv_file.next()
            heading_nums = [[int(x) for x in re.findall('\d+', col)] for col in headings]
            
            # TODO:  check if the headings are as expected, complain if not
            for x in csv_file:
                opts = {}

                opts['params_json'] = json.dumps({'units': 'Thousands', 'source': fname})

                opts['condition'] = '(all)'
                opts['data_type'] = 'population count'
                
                opts['region'] = smart_unicode(x[0].strip(), errors='ignore')
                #opts['iso_code'] = x[1].strip()
                opts['gbd_region'] = x[2].strip()
                opts['year_start'] = int(x[3])
                opts['year_end'] = int(x[3])
                opts['sex'] = x[4].strip().lower()

                for age_list, val_str in zip(heading_nums,x)[5:]:
                    if len(age_list) == 2:
                        a0, a1 = age_list
                        if a1 == 1:  # handle difference in 0-1 and 10-14, etc
                            a1 = 0
                    else:
                        a0 = age_list[0]
                        a1 = MISSING
                    try:
                        opts['age_start'] = a0
                        opts['age_end'] = a1
                        opts['value'] = float(val_str)
                        opts['standard_error'] = 0.0
                        d, is_new = Data.objects.get_or_create(**opts)

                        if is_new:
                            pop_counter += 1

                    except ValueError:
                        pass

            print "added %d rates" % pop_counter

            
