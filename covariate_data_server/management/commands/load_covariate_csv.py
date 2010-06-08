"""
Command to load covariate data from a csv

Use from the project directory as follows::

    $ python2.5 manage.py load_covariate_csv [options] filename.csv

"""

from django.core.management.base import BaseCommand, CommandError

import csv
import re
import pylab as pl

from gbd.covariate_data_server.models import *
from gbd.dismod3.settings import gbd_regions
from gbd.dismod3.utils import clean

class Command(BaseCommand):
    help = 'Import covariate data from a .csv file.'
    args = 'filename.csv'

    def handle(self, *params, **options):
        if len(params) != 2:
            raise CommandError('a single .csv file amd type_slug are required as input.')
        fname = params[0]
        type_slug = params[1]

        print "adding population data from %s" % fname

        csv_f = csv.DictReader(open(fname))
        data = [d for d in csv_f]
        
        if data and type_slug:
            for i, d in enumerate(data):
                try:
                    d['value'] = float(d[type_slug])
                except KeyError:
                    print 'Could not find column %s (is it spelled correctly?)' % type_slug
                    return
                except ValueError:
                    print 'Could not interpret value for %s in line %d' % (type_slug, i+2)
                    return

                if d.has_key('year'):
                    try:
                        d['year'] = int(d['year'])
                    except ValueError:
                        print 'Could not interpret year in line %d' % (ii+2)
                        return
                else:
                    d['year'] = gbd.fields.ALL_YEARS
                        

                d['sex'] = d.get('sex', '')
                if not d['sex'] in ['male', 'female', 'total', '']:
                    print 'Could not interpret sex in line %d' % (ii+2)
                    return

            cov_type, is_new = CovariateType.objects.get_or_create(slug=type_slug)

            # make rates from rate_list
            vals = [d['value'] for d in data]
            mu = pl.mean(vals)
            std = pl.std(vals)

            added = 0
            modified = 0
            for d in data:
                # if sex == '' add a covariate for male, female, and total
                if d['sex'] == '':
                    sex_list = ['male', 'female', 'total']
                else:
                    sex_list = [d['sex']]
                for sex in sex_list:
                    # add a data point, save it on the data list
                    cov, is_new = Covariate.objects.get_or_create(type=cov_type,
                                                                  iso3=d['iso3'],
                                                                  year=d['year'],
                                                                  sex=sex,
                                                                  defaults={'value': 0.})
                    cov.value = (d['value'] - mu) / std
                    cov.save()
                    added += is_new
                    modified += not is_new

        """
        type_slug = 'GDPpc'
        type_desc = 'log(GDP per capita) - mu_log(GDPpc)'
        type, is_new = CovariateType.objects.get_or_create(slug=type_slug, defaults={'description': type_desc})

        vals = []
        for d in data:
            for key in d.keys():
                for year in re.findall('\d+', key):
                    year = int(year)
                    try:
                        value = float(d[key])
                    except ValueError:
                        continue
                    if year > 1900 and year < 2050:
                        vals += [pl.log(value)]
        mu = pl.mean(vals)
        std = pl.std(vals)
        print '%d data points, mean=%.2f, std=%.2f' % (len(vals), mu, std)

        added = 0
        modified = 0
        for d in data:
            iso3 = d['iso3']
            sex = 'total'
            for key in d.keys():
                for year in re.findall('\d+', key):
                    year = int(year)
                    try:
                        value = float(d[key])
                    except ValueError:
                        continue
                    if year > 1900 and year < 2050:
                        cov, is_new = Covariate.objects.get_or_create(type=type, iso3=iso3, year=year, sex=sex, defaults={'value': value})
                        import pdb;pdb.set_trace()
                        cov.value = (pl.log(value) - mu) / std
                        cov.save()
                        added += is_new
                        modified += not is_new
            try:
                print str(cov), cov.value
            except:
                pass
            """
        print 'added %d country-years of covariate data' % added
        print 'modified %d country-years of covariate data' % modified
        
