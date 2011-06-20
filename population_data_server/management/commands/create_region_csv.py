"""
Command to extract the countries in each GBD region from a population csv

Use from the project directory as follows::

    $ python2.5 manage.py create_region_csv USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv

The file USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv is big,
so it takes a few minutes to load.

the global population table csv has first row::
    variant, notes, countrycode, gbd_region, gbd_country, iso3, year,
    sex, pop_0to1, pop_1to4, pop_5to9, pop_10to14, pop_15to19,
    pop_20to24, pop_25to29, pop_30to34, pop_35to39, pop_40to44,
    pop_45to49, pop_50to54, pop_55to59, pop_60to64, pop_65to69,
    pop_70to74, pop_75to79, pop_80plus, pop_80to84, pop_85to89,
    pop_90to94, pop_95to99, pop_100plus
"""

from django.core.management.base import BaseCommand, CommandError
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode

import csv

from gbd.population_data_server.models import Population

gbd_region_col = 3
gbd_country_col = 5

class Command(BaseCommand):
    help = 'Create country-region .csv from population file.'
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        fnames = list(fnames) + ['', '']
        fin = fnames[0] or '../USABLE_IHME_GBD_POPULATION_1950-2010_vWPP2008.csv'
        fout = fnames[1] or 'country_region.csv'

        countries = {}
        all_countries = set()

        f = open(fout, 'w')
        csv_f = csv.writer(f)
        
        data = [d for d in csv.DictReader(open(fin))]
        for d in data:
            if not d['year'] in ['1990', '2005']:
                continue
            if len(d['iso3']) != 3:
                continue
            
            if not countries.has_key(d['gbd_region']):
                countries[d['gbd_region']] = set()

            countries[d['gbd_region']].add(d['iso3'])
            all_countries.add(d['iso3'])

        countries['World'] = all_countries - set(['WLD'])
            
        print 'finished building dict'

        for r in countries.keys():
            csv_f.writerow([r] + sorted(countries[r]))
        f.close()
        
