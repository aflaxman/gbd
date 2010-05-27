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
        holdout = 0
        if a0 > 50 and random.random() < .2:
            holdout = 1
            
        p0 = dismod3.utils.rate_for_range(truth, range(a0, a1 + 1), np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))
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
              'id': len(data),
              'truth': p0,
              'ignore': holdout,
              'test_set': holdout}

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
        'Ignore': d['ignore'],
        'Test Set': d['test_set'],
        'Truth': d['truth'],
        }
    return c



def predict(type, dm, d):
    for k in d.keys():
        d[dismod3.utils.clean(k)] = d[k]
        
    t = d['parameter'].replace(' data', '').replace(' ', '-')
    r = d['region']
    y = int(d['year_start'])
    s = d['sex']
    key = dismod3.gbd_key_for(t, r, y, s)

    a0 = int(d['age_start'])
    a1 = int(d['age_end'])
    est_by_age = dm.get_mcmc(type, key)

    if len(est_by_age) == 0:
        return -99

    est = dismod3.utils.rate_for_range(est_by_age,
                                 range(a0, a1 + 1),
                                 np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))
    d['estimate %s' % type] = est

    return est

def measure_fit_against_gold(id, condition='test_disease_1'):
    """
    Determine the RMSE of the fit stored in model specified by id
    """

    print 'downloading fitted model...'
    sys.stdout.flush()
    dm = dismod3.get_disease_model(id)

    print 'loading gold-standard data'
    gold_data = [d for d in csv.DictReader(open(OUTPUT_PATH + '%s_gold.tsv' % condition), dialect='excel-tab')]


    print 'comparing values'
    abs_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    rel_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    for metric in [abs_err, rel_err, ]:
        metric['excess mortality'] = []

    for d in gold_data:
        est = predict('mean', dm, d)
        if est < 0:
            continue
        val = float(d['Parameter Value'])
        err = val - est


        if d['Age Start'] <= 50:
            continue

        t = d['Parameter'].replace(' data', '')
        abs_err[t].append(err)
        rel_err[t].append(100 * err / val)

    print
    
    for k in abs_err:
        print '%s abs RMSE = %f' % (k, np.sqrt(np.mean(np.array(abs_err[k])**2)))
        print '%s abs  MAE = %f' % (k, np.median(np.abs(abs_err[k])))
    print
    
    for k in rel_err:
        print '%s rel pct RMSE = %f' % (k, np.sqrt(np.mean(np.array(rel_err[k])**2)))
        print '%s rel pct  MAE = %f' % (k, np.median(np.abs(rel_err[k])))
    print

    
    col_names = sorted(set(gold_data[0].keys()) | set(['Estimate Value']))
    f_file = open(OUTPUT_PATH + '%s_gold.tsv' % condition, 'w')
    csv_f = csv.writer(f_file, dialect='excel-tab')
    csv_f.writerow(col_names)
    csv_f = csv.DictWriter(f_file, col_names, dialect='excel-tab')
    for d in gold_data:
        csv_f.writerow(d)
    f_file.close()

    

def measure_fit_against_test_set(id):
    """
    Determine the predictive validity of the fit against data with 1 in 'test_set' column
    """

    print 'downloading fitted model...'
    sys.stdout.flush()
    dm = dismod3.get_disease_model(id)


    print 'comparing values'
    abs_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    rel_err = dict(incidence=[], prevalence=[], remission=[], duration=[])
    coverage = dict(incidence=[], prevalence=[], remission=[], duration=[])

    for metric in [abs_err, rel_err, coverage]:
        metric['excess-mortality'] = []

    for d in dm.data:
        try:
            is_test = int(d['test_set'])
        except (ValueError, KeyError):
            is_test = 0

        if is_test:
            d['region'] = d['gbd_region']
            est = predict('mean', dm, d)
            lb = predict('lower_ui', dm, d)
            ub = predict('upper_ui', dm, d)

            if est < 0 or lb < 0 or ub < 0:
                continue
            
            val = float(d['parameter_value'])
            err = val - est

            t = d['parameter'].replace(' data', '')
            abs_err[t].append(err)
            rel_err[t].append(100 * err / val)
            coverage[t].append(val >= lb and val <= ub)
            #print key, a0, a1, err

    for k in coverage:
        print '%s coverage = %f' % (k, np.sum(coverage[k]) * 100. / len(coverage[k]))
    print
    
    for k in abs_err:
        print '%s abs RMSE = %f' % (k, np.sqrt(np.mean(np.array(abs_err[k])**2)))
        print '%s abs  MAE = %f' % (k, np.median(np.abs(abs_err[k])))
    print
    
    for k in rel_err:
        print '%s rel pct RMSE = %f' % (k, np.sqrt(np.mean(np.array(rel_err[k])**2)))
        print '%s rel pct  MAE = %f' % (k, np.median(np.abs(rel_err[k])))
    print

    



def generate_disease_data(condition='test_disease_3'):
    """ Generate csv files with gold-standard disease data,
    and somewhat good, somewhat dense disease data, as might be expected from a
    condition that is carefully studied in the literature
    """
    
    age_len = dismod3.MAX_AGE
    ages = np.arange(age_len, dtype='float')

    # incidence rate
    #i = .012 * mc.invlogit((ages - 44) / 3)
    i = .01 * (np.ones_like(ages) + ages / age_len)

    # remission rate
    #r = 0. * ages
    r = .7 * np.ones_like(ages) 

    # excess-mortality rate
    #f_init = .085 * (ages / 100) ** 2.5
    SMR = 2. * np.ones_like(ages) - ages / age_len

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

                param_type = 'all-cause_mortality'
                key = dismod3.gbd_key_for(param_type, region, year, sex)
                m_all_cause = mort.mortality(key, mort.data)

                # tweak excess-mortality rate to make rr start at 3.5
                f = (SMR - 1.) * m_all_cause


                ## compartmental model (bins S, C, D, M)
                import scipy.linalg
                from dismod3 import NEARLY_ZERO
                from dismod3.utils import trim

                SCDM = np.zeros([4, age_len])
                p = np.zeros(age_len)
                m = np.zeros(age_len)

                SCDM[0,0] = 1.
                SCDM[1,0] = 0.
                SCDM[2,0] = 0.
                SCDM[3,0] = 0.

                p[0] = SCDM[1,0] / (SCDM[0,0] + SCDM[1,0] + NEARLY_ZERO)
                m[0] = trim(m_all_cause[0] - f[0] * p[0], NEARLY_ZERO, 1-NEARLY_ZERO)

                for a in range(age_len - 1):
                    A = [[-i[a]-m[a],  r[a]          , 0., 0.],
                         [ i[a]     , -r[a]-m[a]-f[a], 0., 0.],
                         [      m[a],       m[a]     , 0., 0.],
                         [        0.,            f[a], 0., 0.]]

                    SCDM[:,a+1] = np.dot(scipy.linalg.expm(A), SCDM[:,a])

                    p[a+1] = SCDM[1,a+1] / (SCDM[0,a+1] + SCDM[1,a+1] + NEARLY_ZERO)
                    m[a+1] = m_all_cause[a+1] - f[a+1] * p[a+1]


                # duration = E[time in bin C]
                hazard = r + m + f
                pr_not_exit = np.exp(-hazard)
                X = np.empty(len(hazard))
                X[-1] = 1 / hazard[-1]
                for ii in reversed(range(len(X)-1)):
                    X[ii] = (pr_not_exit[ii] * (X[ii+1] + 1)) + (1 / hazard[ii] * (1 - pr_not_exit[ii]) - pr_not_exit[ii])

                # shift prevalence to get inconsistent data
                p *= 10

                params = dict(age_intervals=age_intervals, condition=condition, gbd_region=region,
                              country=countries_for[region][0], year=year, sex=sex, effective_sample_size=1.e9)

                generate_and_append_data(gold_data, 'prevalence data', p, **params)
                generate_and_append_data(gold_data, 'incidence data', i, **params)
                generate_and_append_data(gold_data, 'excess-mortality data', f, **params)
                generate_and_append_data(gold_data, 'remission data', r, **params)
                generate_and_append_data(gold_data, 'duration data', X, **params)

                params['effective_sample_size'] = 1000.0
                params['age_intervals'] = sparse_intervals[region]
                generate_and_append_data(noisy_data, 'prevalence data', p, **params)
                generate_and_append_data(noisy_data, 'incidence data', i, **params)
                generate_and_append_data(noisy_data, 'excess-mortality data', f, **params)
                generate_and_append_data(noisy_data, 'remission data', r, **params)



    col_names = sorted(data_dict_for_csv(gold_data[0]).keys())

    f_file = open(OUTPUT_PATH + '%s_gold.tsv' % condition, 'w')
    csv_f = csv.writer(f_file, dialect='excel-tab')
    csv_f.writerow(col_names)
    for d in gold_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()

    f_file = open(OUTPUT_PATH + '%s_data.tsv' % condition, 'w')
    csv_f = csv.writer(f_file, dialect='excel-tab')
    csv_f.writerow(col_names)

    for d in noisy_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()



if __name__ == '__main__':
    generate_disease_data()
