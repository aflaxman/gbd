""" DisMod III Simulation Study - Good, Dense Data
"""

GBD_PATH = '/home/abie/gbd/'
OUTPUT_PATH = '/home/abie/'

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

from dismod3.neg_binom_model import countries_for
import random

def generate_and_append_data(data, data_type, truth, age_intervals, condition,
                             gbd_region, country, year, sex, effective_sample_size):
    """ create simulated data"""
    for a0, a1 in age_intervals:
        d = { 'condition': condition,
              'data_type': data_type,
              'gbd_region': gbd_region,
              'region': country,
              'year_start': year,
              'year_end': year,
              'sex': sex,
              'effective_sample_size': effective_sample_size,
              'age_start': a0,
              'age_end': a1,
              'age_weights': list(np.ones(a1 + 1 - a0)),
              'id': len(data)}

        p0 = dismod3.utils.rate_for_range(truth, range(a0, a1 + 1), np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))

        if p0 < 1:
            d['value'] = mc.rbinomial(d['effective_sample_size'], p0) / float(d['effective_sample_size'])
        else:
            d['value'] = p0

        
        data.append(d)

def data_dict_for_csv(d):
    c = {
        'GBD Cause': d['condition'],
        'Parameter': d['data_type'].replace('-', ' '),
        'Country ISO3 Code': d['region'],
        'Region': d['gbd_region'],
        'Parameter Value': d['value'],
        'Standard Error': '',
        'Effective Sample Size': d['effective_sample_size'],
        'Units': 1.0,
        'Sex': d['sex'],
        'Age Start': d['age_start'],
        'Age End': d['age_end'],
        'Year Start': d['year_start'],
        'Year End': d['year_end'],
        }
    return c


def generate_disease_data(condition='test_disease_1'):
    """ Generate csv files with gold-standard disease data,
    and somewhat good, somewhat dense disease data, as might be expected from a
    condition that is carefully studied in the literature
    """
    
    age_len = dismod3.MAX_AGE
    ages = np.arange(age_len, dtype='float')

    # incidence rate
    i = .012 * mc.invlogit((ages - 44) / 3)

    # remission rate
    r = 0. * ages

    # excess-mortality rate
    f_init = .085 * (ages / 100) ** 2.5

    # all-cause mortality-rate
    mort = dismod3.get_disease_model('all-cause_mortality')

    age_intervals = [[a, a+4] for a in range(0, dismod3.MAX_AGE-4, 5)]
    sparse_intervals = dict([[region, random.sample(age_intervals, (ii**2 * len(age_intervals)) / len(countries_for)**2)] for ii, region in enumerate(countries_for)])

    gold_data = []
    noisy_data = []

    for region in countries_for:
        if region == 'world':
            continue
        
        print region
        sys.stdout.flush()
        for year in [1990, 2005]:
            for sex in ['male', 'female']:
                data = [d for d in mort.data if d['data_type'] == 'all-cause mortality data'
                        and d['region'] == region and d['sex'] == sex and d['year_start'] == year]
                m_all_cause = mort.mortality('all_cause', data)

                # tweak excess-mortality rate to make rr start at 3.5
                f = f_init + m_all_cause * 2.5 * np.maximum((40-ages)/40, 0)


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
                pr_exit = np.exp(- r - m - f)
                X = np.empty(len(pr_exit))
                t = 1.
                for a in xrange(len(X) - 1, -1, -1):
                    X[a] = t * pr_exit[a]
                    t = 1 + X[a]

                params = dict(age_intervals=age_intervals, condition=condition, gbd_region=region,
                              country=countries_for[region][0], year=year, sex=sex, effective_sample_size=1.e9)

                generate_and_append_data(gold_data, 'prevalence data', p, **params)
                generate_and_append_data(gold_data, 'incidence data', i, **params)
                generate_and_append_data(gold_data, 'excess-mortality data', f, **params)
                generate_and_append_data(gold_data, 'remission data', r, **params)
                generate_and_append_data(gold_data, 'duration data', X, **params)

                params['effective_sample_size'] = 10000.0
                params['age_intervals'] = sparse_intervals[region]
                generate_and_append_data(noisy_data, 'prevalence data', p, **params)
                generate_and_append_data(noisy_data, 'incidence data', i, **params)



    col_names = sorted(data_dict_for_csv(gold_data[0]).keys())

    f_file = open(OUTPUT_PATH + '%s_gold.csv' % condition, 'w')
    csv_f = csv.writer(f_file)
    csv_f.writerow(col_names)
    for d in gold_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()

    f_file = open(OUTPUT_PATH + '%s_data.csv' % condition, 'w')
    csv_f = csv.writer(f_file)
    csv_f.writerow(col_names)

    for d in noisy_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()

def measure_fit(id, condition='test_disease_1'):
    """
    Determine the RMSE of the fit stored in model specified by id
    """

    print 'downloading fitted model...'
    sys.stdout.flush()
    dm = dismod3.get_disease_model(id)

    print 'loading gold-standard data'
    gold_data = [d for d in csv.DictReader(open(OUTPUT_PATH + '%s_gold.csv' % condition))]


    print 'comparing values'
    abs_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    rel_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    for d in gold_data:
        t = d['Parameter'].replace(' data', '')
        r = d['Region']
        y = int(d['Year Start'])
        s = d['Sex']
        key = dismod3.gbd_key_for(t, r, y, s)

        est_by_age = dm.get_mcmc('mean', key)
        a0 = int(d['Age Start'])
        a1 = int(d['Age End'])
        est_by_age = dm.get_mcmc('mean', key)

        if len(est_by_age) == 0:
            continue
        
        est = dismod3.utils.rate_for_range(est_by_age,
                                     range(a0, a1 + 1),
                                     np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))
        d['Estimate Value'] = est
        
        val = float(d['Parameter Value'])
        err = val - est
        abs_err[t].append(err)
        rel_err[t].append(100 * err / val)
        #print key, a0, a1, err

    print
    
    for k in abs_err:
        print '%s abs RMSE = %f' % (k, np.sqrt(np.mean(np.array(abs_err[k])**2)))

    for k in abs_err:
        print '%s rel pct MAE = %f' % (k, np.median(np.abs(rel_err[k])))

    col_names = sorted(set(gold_data[0].keys()) | set(['Estimate Value']))
    f_file = open(OUTPUT_PATH + '%s_gold.csv' % condition, 'w')
    csv_f = csv.writer(f_file)
    csv_f.writerow(col_names)
    csv_f = csv.DictWriter(f_file, col_names)
    for d in gold_data:
        csv_f.writerow(d)
    f_file.close()
                                        
if __name__ == '__main__':
    generate_disease_data()
