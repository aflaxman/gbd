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
        if len(params) == 1:
            fname= params[0]
            type_slug = params[0].split('/')[-1].replace('.csv','')
        elif len(params) == 2:
            fname = params[0]
            type_slug = params[1]
        else:
            raise CommandError('a single .csv file amd type_slug are required as input.')

        print "adding %s covariate data from %s" % (type_slug, fname)

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
            for i, d in enumerate(data):
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
                if i % (len(data) / 100) == 0:
                    print cov_type, (i * 100) / len(data), '% done'
        print 'added %d country-years of covariate data' % added
        print 'modified %d country-years of covariate data' % modified
        
