"""
load data csv
=============

Django command that loads the contents of command-line specified csv
files as rates for dismod 3.

command options are:
  -r RATE_TYPE, --rate_type=RATE_TYPE
                        The rate_type for the data in the csv file,
                        where RATE_TYPE = incidence | prevalence |
                                          remission | case_fatality |
                                          all_cause
  -d DISEASE, --disease=DISEASE
                        The disease for the data in the csv file

"""
from django.utils.encoding import DjangoUnicodeDecodeError, smart_unicode
from django.core.management.base import BaseCommand, CommandError

import optparse
import simplejson as json

from dismod3.models import Rate, Region, Disease, Population

rate_dict = {
    'incidence': 'incidence data',
    'prevalence': 'prevalence data',
    'remission': 'remission data',
    'case_fatality': 'case fatality data',
    'all_cause': 'all-cause mortality data',
}

def try_int(str):
    try:
        ret_val = int(str)
    except ValueError:
        ret_val = -1
    return ret_val

def try_float(str):
    try:
        ret_val = float(str)
    except ValueError:
        ret_val = -1.
    return ret_val


class Command(BaseCommand):
    option_list = BaseCommand.option_list + (
        optparse.make_option('-r', '--rate_type',
                             default='prevalence',
                             help="The rate_type for the data in the csv file, where RATE_TYPE = incidence | prevalence | remission | case_fatality",
                             ),
        optparse.make_option('-d', '--disease',
                             help='The disease for the data in the csv file',
                             ),
        )
                             
    help = "Import rate data from a .csv file."
    args = 'filename.csv'

    def handle(self, *fnames, **options):
        import csv

        if len(fnames) != 1:
            raise CommandError('a single .csv file is required as input.')
        fname = fnames[0]

        try:
            rate_type = rate_dict[options.get('rate_type')]
        except KeyError:
            raise CommandError('unrecognized rate type.')
        
        
        print "adding rates from %s" % fname
        disease, created = Disease.objects.get_or_create(name=options.get('disease') or fname)
        
        rate_counter = 0

        csv_file = csv.reader(open(fname))
        headings = csv_file.next()
        # TODO:  check if the headings are as expected, complain if not

        for ii,row in enumerate(csv_file):
            country = smart_unicode(row[2].strip(), errors='ignore')
            pop = Population.objects.filter(country=country)
            if pop.count() > 0:
                region = pop[0].region
            else:
                region, created = Region.objects.get_or_create(name='World')

            params = {}
            for key, val in zip(headings, row):
                params[key] = val

            offset=0
            for sex in ['male', 'female', 'total']:
                try:
                    numerator, denominator, offset = get_rates_from_row(row, sex)
                    rate_data, created = Rate.objects.get_or_create(
                        disease=disease,
                        region=region,
                        rate_type=rate_type,
                        sex=sex,
                        country=country,
                        epoch_start=int(row[6]),
                        epoch_end=int(row[7]), 
                        age_start=int(row[19]),
                        age_end=int(row[20]),
                        numerator=numerator,
                        denominator=denominator,
                        defaults={'params_json': json.dumps(params)}
                        )
                    rate_counter += 1
                        
                except ValueError:
                    pass

        print "added %d rates" % rate_counter
                    
def get_rates_from_row(row, sex):
    from numpy import round
    
    offset_dict = { 'female': 21, 'male': 31, 'total': 41 }
    offset = offset_dict[sex]

    numerator   = try_int(row[offset+5])
    denominator = try_int(row[offset+6])

    if numerator == -1 or denominator == -1:
        rate_estimate = try_float(row[offset+0])
        if rate_estimate == -1.:
            raise ValueError, "couldn't find rate data."
        
        rate_radix = try_float(row[offset+7])
        if rate_radix == -1.:
            rate_radix = 100.

        if numerator == -1 and denominator != -1:
            numerator = int(round(denominator * (rate_estimate / rate_radix)))

        elif numerator != -1 and denominator == -1:
            denominator = int(round(numerator / (rate_estimate / rate_radix)))

        #TODO:  make more reasonable estimate in this case
        else:  # numerator == -1 and denominator == -1:
            denominator      = 50
            numerator        = int(round(denominator * rate_estimate/rate_radix))


    return numerator, denominator, offset

            
