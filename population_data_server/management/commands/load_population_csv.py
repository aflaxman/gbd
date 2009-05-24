"""
Command to load the GBD population data from a csv

Use from the project directory as follows::

    $ python2.5 manage.py load_population_csv USABLE_IHME_GBD_POPULATION_1950-2050.csv

The file USABLE_IHME_GBD_POPULATION_1950-2050.csv is big,
so it takes 2 minutes to load.

the global population table csv is a file formatted to
match the gbd format, which has first row::

    gbd_country, iso3, gbd_region, year, sex, pop_0to1, pop_1to4, pop_5to9, pop_10to14, \
pop_15to19, pop_20to24, pop_25to29, pop_30to34, pop_35to39, pop_40to44, pop_45to49, \
pop_50to54, pop_55to59, pop_60to64, pop_65to69, pop_70to74, pop_75to79, pop_80to84, \
pop_85to89, pop_80plus, pop_90plus
"""

from django.core.management.base import BaseCommand
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import simplejson as json
import re

from gbd.population_data_server.models import Population

MAX_AGE=100

class Command(BaseCommand):
    help = 'Import population data from a .csv file.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        import csv
        import numpy as np
        
        if not fnames:
            fnames = ['../USABLE_IHME_GBD_POPULATION_1950-2050.csv']

        gbd_region_pop = {}
        
        for fname in fnames:
            print "adding population data from %s" % fname
            pop_counter = 0
            
            csv_file = csv.reader(open(fname))
            headings = csv_file.next()
            heading_nums = [[int(x) for x in re.findall('\d+', col)] for col in headings]
            
            # TODO:  check if the headings are as expected, complain if not
            for x in csv_file:
                opts = {}
                opts['region'] = smart_unicode(x[0].strip(), errors='ignore')
                opts['year'] = int(x[3])
                opts['sex'] = x[4].strip().lower()

                mesh = []
                vals = []
                interval_start = []
                interval_length = []
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
                        interval_start.append(a0)
                        interval_length.append(a1-a0)
                    except ValueError:
                        pass

                opts['params_json'] = json.dumps({'mesh': list(mesh),
                                                  'vals': list(vals),
                                                  'interval_start': list(interval_start),
                                                  'interval_length': list(interval_length)})
                pop, is_new = Population.objects.get_or_create(**opts)
                if not is_new:
                    pop.save()
                    pop_counter += 1

                region = x[2]
                key = opts_to_key(opts)
                
                M,C = pop.gaussian_process()

                if not gbd_region_pop.has_key(region):
                    gbd_region_pop[region] = {}

                if not gbd_region_pop[region].has_key(key):
                    gbd_region_pop[region][key] = np.zeros(MAX_AGE)

                gbd_region_pop[region][key] += M(range(MAX_AGE))

            print "added %d rates" % pop_counter

            for region in gbd_region_pop.keys():
                for key in gbd_region_pop[region].keys():
                    opts = keys_to_opts(region, key)
                    
                    opts['params_json'] = json.dumps({'mesh': range(MAX_AGE),
                                                      'vals': list(gbd_region_pop[region][key])})

                    pop, is_new = Population.objects.get_or_create(**opts)

def opts_to_key(opts):
    return (opts['year'], opts['sex'])

def keys_to_opts(region, key):
    return {'region': region, 'year': key[0], 'sex': key[1]}
