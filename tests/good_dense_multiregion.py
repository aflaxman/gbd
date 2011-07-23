""" DisMod III Simulation Study - Good, Dense Data
"""

GBD_PATH = '/home/abie/gbd/'
OUTPUT_PATH = '/home/j/temp/'
OUTPUT_FILE = 'good_dense_data.csv'

import sys
import os
import optparse
import csv

import pylab as pl
import numpy as np
import pymc as mc
import simplejson as json

sys.path.append(GBD_PATH)
import dismod3
import dismod3.utils
from dismod3.disease_json import DiseaseJson
import dismod3.gbd_disease_model as model

from covariate_data_server.models import Covariate

age_len = dismod3.MAX_AGE
ages = np.arange(age_len, dtype='float')

def generate_and_append_data(data, data_type, truth, age_intervals,
                             gbd_region, country, year, sex):
    """ create simulated data"""
    for a0, a1 in age_intervals:
        d = { 'condition': 'type_2_diabetes',
              'data_type': data_type,
              'gbd_region': gbd_region,
              'region': country,
              'year_start': year,
              'year_end': year,
              'sex': sex,
              'age_start': a0,
              'age_end': a1,
              'age_weights': list(np.ones(a1 + 1 - a0)),
              'id': len(data)}

        p0 = dismod3.utils.rate_for_range(truth, range(a0, a1 + 1), np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))
    
        d['value'] = p0
        if p0 == 0.:
            d['standard_error'] = .000001
        elif p0 < 1.:
            d['standard_error'] = p0 * (1-p0) / np.sqrt(1000)
        else:
            d['standard_error'] = p0 * .05

        data.append(d)

data = []
mort = dismod3.get_disease_model('all-cause_mortality')

age_intervals = [[a, a+4] for a in range(0, dismod3.MAX_AGE-4, 5)]
for sex in ['male', 'female']:
    for year in [1990, 2005]:
        for region, country in \
                [['Asia, Southeast', 'THA'],
                 ['Asia, East', 'CHN'],
                 ['Asia, South', 'IND'],
                 ['Europe, Central', 'ALB']]:

            offset = 0.
            if sex == 'male':
                offset += .5
            if year == 2005:
                offset += .5
            if region == 'Asia, South':
                offset -= .1
            if region == 'Asia, East':
                offset -= .2
            if region == 'Europe, Central':
                offset += .3

            gdp = Covariate.objects.get(iso3=country, year=year).value
            offset += .3*gdp
            
            # incidence rate
            i = .012 * mc.invlogit((ages - 44) / 3) * (1 + offset)

            # remission rate
            r = 0. * ages

            # excess-mortality rate
            f = .085 * (ages / 100) ** 2.5

            # all-cause mortality-rate
            mort_data = [d for d in mort.data if d['data_type'] == 'all-cause mortality data'
                    and d['region'] == region and d['sex'] == sex and d['year_start'] == year]
            m_all_cause = mort.mortality('all_cause', mort_data)

            # tweak excess-mortality rate to make rr start at 3.5
            f += m_all_cause * 2.5 * np.maximum((40-ages)/40, 0)


            ## compartmental model (bins S, C, D, M)
            import scipy.linalg
            from dismod3 import NEARLY_ZERO
            from dismod3.utils import trim

            SCDM = np.zeros([4, age_len])
            p = np.zeros(age_len)
            m = np.zeros(age_len)

            SCDM[0,0] = 1.
            SCDM[1,0] = 0.
            SCDM[2,0] = NEARLY_ZERO
            SCDM[3,0] = NEARLY_ZERO

            p[0] = SCDM[1,0] / (SCDM[0,0] + SCDM[1,0] + NEARLY_ZERO)
            m[0] = trim(m_all_cause[0] - f[0] * p[0], NEARLY_ZERO, 1-NEARLY_ZERO)

            for a in range(age_len - 1):
                A = [[-i[a]-m[a],  r[a]          , 0., 0.],
                     [ i[a]     , -r[a]-m[a]-f[a], 0., 0.],
                     [      m[a],       m[a]     , 0., 0.],
                     [        0.,            f[a], 0., 0.]]

                SCDM[:,a+1] = np.dot(scipy.linalg.expm(A), SCDM[:,a])

                p[a+1] = SCDM[1,a+1] / (SCDM[0,a+1] + SCDM[1,a+1] + NEARLY_ZERO)
                m[a+1] = trim(m_all_cause[a+1] - f[a+1] * p[a+1], .1*m_all_cause[a+1], 1-NEARLY_ZERO)


            # duration = E[time in bin C]
            pr_exit = 1 - r - m - f
            X = np.empty(len(pr_exit))
            t = 1.
            for a in xrange(len(X) - 1, -1, -1):
                X[a] = t * pr_exit[a]
                t = 1 + X[a]

            generate_and_append_data(data, 'prevalence data', p, age_intervals, year=year, sex=sex, gbd_region=region, country=country)
            generate_and_append_data(data, 'incidence data', i, age_intervals, year=year, sex=sex, gbd_region=region, country=country)
            generate_and_append_data(data, 'relative-risk data', (m+f)/m, age_intervals, year=year, sex=sex, gbd_region=region, country=country)
            generate_and_append_data(data, 'remission data', r, age_intervals, year=year, sex=sex, gbd_region=region, country=country)
            generate_and_append_data(data, 'duration data', X, age_intervals, year=year, sex=sex, gbd_region=region, country=country)

def data_dict_for_csv(d):
    c = {
        'GBD Cause': d['condition'],
        'Parameter': d['data_type'].replace('-', ' '),
        'Country ISO3 Code': d['region'],
        'Region': d['gbd_region'],
        'Parameter Value': d['value'],
        'Standard Error': d['standard_error'],
        'Units': 1.0,
        'Sex': d['sex'],
        'Age Start': d['age_start'],
        'Age End': d['age_end'],
        'Year Start': d['year_start'],
        'Year End': d['year_end'],
        }
    return c
        
f_file = open(OUTPUT_PATH + 'simulated_data.tsv', 'w')
csv_f = csv.writer(f_file, dialect=csv.excel_tab)

col_names = sorted(data_dict_for_csv(data[0]).keys())
        
csv_f.writerow(col_names)
for d in data:
    dd = data_dict_for_csv(d)
    csv_f.writerow([dd[c] for c in col_names])
f_file.close()
