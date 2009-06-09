"""
Command to load the GBD population data from a csv

Use from the project directory as follows
$ python2.5 manage.py load_population

The file USABLE_IHME_GBD_POPULATION_1950-2050.csv is big,
so it takes 2 minutes to load.
"""

from django.core.management.base import BaseCommand
from dismod3.models import Region, Population
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import simplejson as json
import re

MAX_AGE=100.

class Command(BaseCommand):
    """
    load the global population tables from a csv file formatted to
    match the gbd format, which has header:
    gbd_country,iso3,gbd_region,year,sex,pop_0to1,pop_1to4,pop_5to9,pop_10to14,pop_15to19,pop_20to24,pop_25to29,pop_30to34,pop_35to39,pop_40to44,pop_45to49,pop_50to54,pop_55to59,pop_60to64,pop_65to69,pop_70to74,pop_75to79,pop_80to84,pop_85to89,pop_80plus,pop_90plus
    """
    help = "Import population data from a .csv file."
    args = 'filename.csv'

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
                opts['country'] = smart_unicode(x[0].strip(), errors='ignore')
                opts['iso_code'] = x[1].strip()
                opts['region'], created = Region.objects.get_or_create(name=x[2].strip())
                opts['year'] = int(x[3])
                opts['sex'] = x[4].strip().lower()

                mesh = []
                vals = []
                for age_list, val_str in zip(heading_nums,x)[5:]:
                    if len(age_list) == 2:
                        a0, a1 = age_list
                        if a1 > 1:  # handle difference in 0-1 and 10-14, etc
                            a1 += 1
                    else:
                        a0 = age_list[0]
                        a1 = MAX_AGE
                    try:
                        vals.append(float(val_str) / float(a1 - a0))
                        mesh.append(.5 * float(a0 + a1))
                    except ValueError:
                        pass

                opts['data_json'] = json.dumps({'mesh': list(mesh),
                                                'vals': list(vals)})

                pop_by_age = Population(**opts)
                pop_by_age.save()
                pop_counter += 1

            print "added %d rates" % pop_counter

            
