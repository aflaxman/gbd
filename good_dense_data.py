""" DisMod III Simulation Study - Good, Dense Data

Example Usage
-------------

# To generate and run data
In [1]: import good_dense_data
In [2]: import fit_all
In [3]: import dismod3
In [4]: for ii in range(4290, 4300):
...:     good_dense_data.generate_disease_data()
...:     dismod3.add_covariates_to_disease_model(ii)
...:     fit_all.fit_all(ii)

# To see how it worked:
In [12]: import good_dense_data
In [13]: y=[good_dense_data.measure_fit_against_gold(i, 'test_disease_07_22_2010b') for i in range(4290, 4300)]
"""

import sys
import os
import optparse
import csv

import pylab as pl
import numpy as np
import pymc as mc
import simplejson as json

import dismod3
import dismod3.utils
from dismod3.disease_json import DiseaseJson
import dismod3.gbd_disease_model as model

from dismod3 import NEARLY_ZERO
from dismod3.neg_binom_model import countries_for, population_by_age, regional_population
import random

GBD_PATH = os.getcwd() + '/../'
sys.path.append(GBD_PATH)

OUTPUT_PATH = GBD_PATH


def generate_disease_data(condition='test_disease_08_30_2010'):
    """ Generate csv files with gold-standard disease data,
    and somewhat good, somewhat dense disease data, as might be expected from a
    condition that is carefully studied in the literature
    """
    
    age_len = dismod3.MAX_AGE
    ages = np.arange(age_len, dtype='float')

    # incidence rate
    #i0 = .0012 * mc.invlogit((ages - 44) / 3)
    i0 = np.maximum(0., .001 * (-.125 + np.ones_like(ages) + (ages / age_len)**2.))

    # remission rate
    #r = 0. * ages
    r = .07 * np.ones_like(ages)

    # excess-mortality rate
    #f_init = .085 * (ages / 100) ** 2.5
    SMR = 3. * np.ones_like(ages) - ages / age_len

    # all-cause mortality-rate
    mort = dismod3.get_disease_model('all-cause_mortality')

    age_intervals = [[a, a+9] for a in range(0, dismod3.MAX_AGE-4, 10)] + [[0, 100] for ii in range(10)]
    #age_intervals = [[a, a+2] for a in range(5, dismod3.MAX_AGE-4, 10)]
    
    # TODO:  take age structure from real data
    sparse_intervals = dict([[region, random.sample(age_intervals, (ii**2 * len(age_intervals)) / len(countries_for)**2 / 1)] for ii, region in enumerate(countries_for)])
    #dense_intervals = dict([[region, random.sample(age_intervals, .5)] for ii, region in enumerate(countries_for)])

    gold_data = []
    noisy_data = []
            
    for ii, region in enumerate(sorted(countries_for)):
        if region == 'world':
            continue
        
        print region
        sys.stdout.flush()

        i = i0 * (1 + float(ii) / 21)
        
        for year in [1990, 2005]:
            for sex in ['male', 'female']:

                param_type = 'all-cause_mortality'
                key = dismod3.gbd_key_for(param_type, region, year, sex)
                m_all_cause = mort.mortality(key, mort.data)

                # calculate excess-mortality rate from smr
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

                country = countries_for[region][0]
                params = dict(age_intervals=age_intervals, condition=condition, gbd_region=region,
                              country=country, year=year, sex=sex, effective_sample_size=1.e9, snr=1.e9)

                params['age_intervals'] = [[0,99]]
                generate_and_append_data(gold_data, 'prevalence data', p, **params)
                generate_and_append_data(gold_data, 'incidence data', i, **params)
                generate_and_append_data(gold_data, 'excess-mortality data', f, **params)
                generate_and_append_data(gold_data, 'remission data', r, **params)
                generate_and_append_data(gold_data, 'duration data', X, **params)

                # TODO: use this approach to age standardize all gold data, and then change it to get iX as a direct sum
                params['age_intervals'] = [[0,99]]
                iX = i * X * regional_population(key)
                generate_and_append_data(gold_data, 'incidence_x_duration', iX, **params)
                

                params['effective_sample_size'] = 1000.
                params['snr'] = 100.
                params['age_intervals'] = sparse_intervals[region]
                generate_and_append_data(noisy_data, 'prevalence data', p, **params)
                generate_and_append_data(noisy_data, 'excess-mortality data', f, **params)
                generate_and_append_data(noisy_data, 'remission data', r, **params)

                params['age_intervals'] = age_intervals
                generate_and_append_data(noisy_data, 'incidence data', i, **params)



    col_names = sorted(data_dict_for_csv(gold_data[0]).keys())

    f_file = open(OUTPUT_PATH + '%s_gold.tsv' % condition, 'w')
    csv_f = csv.writer(f_file, dialect='excel-tab')
    csv_f.writerow(col_names)
    for d in gold_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()

    f_name = OUTPUT_PATH + '%s_data.tsv' % condition
    f_file = open(f_name, 'w')
    csv_f = csv.writer(f_file, dialect='excel-tab')
    csv_f.writerow(col_names)

    for d in noisy_data:
        dd = data_dict_for_csv(d)
        csv_f.writerow([dd[c] for c in col_names])
    f_file.close()

    # upload data file
    from dismod3.disease_json import dismod_server_login, twc, DISMOD_BASE_URL
    dismod_server_login()
    twc.go(DISMOD_BASE_URL + 'dismod/data/upload/')
    twc.formvalue(1, 'tab_separated_values', open(f_name).read())

    try:
        url = twc.submit()
    except Exception, e:
        print e


def generate_and_append_data(data, data_type, truth, age_intervals, condition,
                             gbd_region, country, year, sex, effective_sample_size, snr):
    """ create simulated data"""
    for a0, a1 in age_intervals:
        d = { 'condition': condition,
              'data_type': data_type,
              'gbd_region': gbd_region,
              'region': country,
              'year_start': year,
              'year_end': year,
              'sex': sex,
              'age_start': a0,
              'age_end': a1,
              'id': len(data),}

        holdout = 0
        d['ignore'] = holdout
        d['test_set'] = holdout

        ages = range(a0, a1 + 1)
        pop = np.array([population_by_age[(country, str(year), sex)][a] for a in ages])
        if np.sum(pop) > 0:
            pop /= float(np.sum(pop))  # normalize the pop weights to sum to 1
        else:
            pop = np.ones_like(ages) / float(len(ages))  # for countries where pop is zero, fill in constant structure
        d['age_weights'] = list(pop)

        p0 = dismod3.utils.rate_for_range(truth, ages, pop)
        d['truth'] = p0

        if p0 == 0:
            p1 = 0
        else:
            p1 = mc.rtruncnorm(p0, snr * p0**-2, 0, np.inf)
            assert not np.isnan(p1)
        
        if p1 == 0.:
            import pdb; pdb.set_trace() # are there zeros in the data?  if so why?
            print 'zeros found'
            d['value'] = p1
            d['effective_sample_size'] = effective_sample_size
        else:
            # add noise to the data
            d['value'] = p1
            d['standard_error'] = p0 / np.sqrt(snr)
        data.append(d)

def data_dict_for_csv(d):
    c = {
        'GBD Cause': d['condition'],
        'Parameter': d['data_type'].replace('-', ' '),
        'Country ISO3 Code': d['region'],
        'Region': d['gbd_region'],
        'Parameter Value': d['value'],
        'Standard Error': d.get('standard_error', ''),
        'Effective Sample Size': d.get('effective_sample_size', ''),
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

    ages = range(a0, a1 + 1)
    #pop = np.ones(a1 + 1 - a0) / float(a1 + 1 - a0))
    c = d['country_iso3_code']
    pop = [population_by_age[(c, str(y), s)][a] for a in ages]
    pop /= np.sum(pop)  # normalize the pop weights to sum to 1

    est = dismod3.utils.rate_for_range(est_by_age, ages, pop)
    d['estimate %s' % type] = est

    return est

def measure_fit_against_gold(id, condition):
    """
    Determine the RMSE of the fit stored in model specified by id
    """

    print 'downloading model %d' % id
    sys.stdout.flush()
    dm = dismod3.load_disease_model(id)

    #print 'loading gold-standard data'
    gold_data = [d for d in csv.DictReader(open(OUTPUT_PATH + '%s_gold.tsv' % condition), dialect='excel-tab')]


    #print 'comparing values'
    abs_err = dict(incidence=[], prevalence=[], remission=[], duration=[], incidence_x_duration=[])
    rel_err = dict(incidence=[], prevalence=[], remission=[], duration=[], incidence_x_duration=[])
    for metric in [abs_err, rel_err, ]:
        metric['excess mortality'] = []

    for d in gold_data:
        est = predict('mean', dm, d)
        if est < 0:
            continue
        val = float(d['Parameter Value'])
        err = val - est


        #if d['Age Start'] <= 50:
        #    continue

        t = d['Parameter'].replace(' data', '')
        abs_err[t].append(err)
        rel_err[t].append(100 * err / val)

    for k in abs_err:
        print '%s abs RMSE = %f' % (k, np.sqrt(np.mean(np.array(abs_err[k])**2)))
        print '%s abs  MAE = %f' % (k, np.median(np.abs(abs_err[k])))
    print
    
    for k in rel_err:
        print '%s rel pct RMSE = %f' % (k, np.sqrt(np.mean(np.array(rel_err[k])**2)))
        print '%s rel pct  MAE = %f' % (k, np.median(np.abs(rel_err[k])))
    print


    k = 'incidence_x_duration'
    print '%s rel pct MAE =\t%f' % (k, np.median(np.abs(rel_err[k])))
    return np.median(np.abs(rel_err[k]))


# results of simulation for n=300, cv=20
# y=[good_dense_data.measure_fit_against_gold(i, 'test_disease_7') for i in range(4022, 4036) + range(4037,4040) + [4065]]
# In [82]: sort(y)[[5,10,15]]
# Out[82]: array([  7.84227934,  10.7896475 ,  15.96722595])

# results of simulation for n=2048, cv=2
# y=[good_dense_data.measure_fit_against_gold(i, 'test_disease_8') for i in [4088, 4089, 4090, 4091, 4092, 4093, 4094, 4095, 4096, 4097, 4099, 4100, 4104, 4105, 4106, 4107, 4112, 4113, 4114, 4115]]
# In [82]: sort(y)[[5,10,15]]
# Out[82]: array([  7.84227934,  10.7896475 ,  15.96722595])


    # add estimate value as a column in the gold data tsv, for looking
    # in more detail with a spreadsheet or different code
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

if __name__ == '__main__':
    generate_disease_data()
