"""
Command to load the GBD population data from a csv

Use from the project directory as follows::

    $ python2.5 manage.py load_population_csv USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv

The file USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv is big,
so it takes a few minutes to load.

the global population table csv has first row::
    variant, notes, countrycode, gbd_region, gbd_country, iso3, year,
    sex, poptotal, pop_0to1, pop_1to4, pop_5to9, pop_10to14, pop_15to19,
    pop_20to24, pop_25to29, pop_30to34, pop_35to39, pop_40to44,
    pop_45to49, pop_50to54, pop_55to59, pop_60to64, pop_65to69,
    pop_70to74, pop_75to79, pop_80plus, pop_80to84, pop_85to89,
    pop_90to94, pop_95to99, pop_100plus
"""

from django.core.management.base import BaseCommand
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import numpy as np
import re
import csv
import simplejson as json

from gbd.population_data_server.models import Population
from gbd.dismod3.utils import MAX_AGE

gbd_region_col = 3
gbd_country_col = 5
year_col = 6
sex_col = 7
pop_val_start_col = 9

sex_str = ['total', 'male', 'female']

class Command(BaseCommand):
    help = 'Import population data from a .csv file.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        if not fnames:
            fnames = ['../USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv']

        gbd_region_pop = {}
        
        for fname in fnames:
            print "adding population data from %s" % fname
            pop_counter = 0

            csv_file = csv.reader(open(fname))
            headings = csv_file.next()

            assert headings[gbd_country_col] == 'iso3'
            assert headings[year_col] == 'year'
            assert headings[sex_col] == 'sex'
            
            heading_nums = [[int(x) for x in re.findall('\d+', col)] for col in headings]
            
            for x in csv_file:
                opts = {}
                opts['region'] = x[gbd_country_col]
                opts['year'] = int(x[year_col])
                opts['sex'] = x[sex_col]
                assert opts['sex'] in sex_str
                Population.objects.filter(**opts).delete()
                pop = Population.objects.create(**opts)

                mesh = []
                vals = []
                interval_start = []
                interval_length = []
                for age_list, val_str in zip(heading_nums,x)[pop_val_start_col:]:
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

                pop.params = {'mesh': list(mesh),
                              'vals': list(vals),
                              'interval_start': list(interval_start),
                              'interval_length': list(interval_length)}
                pop.cache_params()
                pop.save()
                pop_counter += 1

                region = x[gbd_region_col]
                
                if not gbd_region_pop.has_key(region):
                    gbd_region_pop[region] = {}

                key = opts_to_key(opts)
                if not gbd_region_pop[region].has_key(key):
                    gbd_region_pop[region][key] = np.zeros(MAX_AGE)

                gbd_region_pop[region][key] += pop.interpolate(range(MAX_AGE))
                print pop
                
            print "added %d rates" % pop_counter

            for region in gbd_region_pop.keys():
                for key in gbd_region_pop[region].keys():
                    opts = keys_to_opts(region, key)
                    
                    Population.objects.filter(**opts).delete()
                    pop = Population.objects.create(**opts)
                    pop.params_json = json.dumps({'mesh': range(MAX_AGE),
                                                  'vals': list(gbd_region_pop[region][key])})
                    pop.save()


def opts_to_key(opts):
    return (opts['year'], opts['sex'])

def keys_to_opts(region, key):
    return {'region': region, 'year': key[0], 'sex': key[1]}
