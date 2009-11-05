""" DisMod III Covariate Simulation Study

For more realistic simulated data, use the age, year, and country data
from a specified disease model.  For example:

    $ python2.5 generate_covariate_data.py 175

"""

import sys
import os
import optparse
import csv
import random

import pylab as pl
import numpy as np
import scipy.linalg
import pymc as mc
import simplejson as json

import dismod3
import dismod3.utils
from dismod3.utils import clean, trim, NEARLY_ZERO
from dismod3.disease_json import DiseaseJson
import dismod3.gbd_disease_model as model


# simulation parameters
usage = 'usage: %prog [options] disease_model_id'
parser = optparse.OptionParser(usage)

parser.add_option('-n', '--studysize', dest='study_size', default='10000',
                  help='number of subjects in each age-range')
parser.add_option('-d', '--dispersion', dest='dispersion', default='10000',
                  help='dispersion of study-level beta-binomial in data generation process')

(options, args) = parser.parse_args()


# check that args are correct
if len(args) == 1:
    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')
        exit()
else:
    parser.error('incorrect number of arguments')
    exit()


# fetch requested model
dm = dismod3.get_disease_model(id)


# define ground truth
age_len = dismod3.MAX_AGE
ages = np.arange(age_len, dtype='float')

print 'defining model transition parameters'

truth = {}

# all-cause mortality-rate
m = np.array(
    [ 0.03266595,  0.01114646,  0.00450302,  0.00226896,  0.00143311,
      0.00109108,  0.00094584,  0.00087981,  0.00083913,  0.0008073 ,
      0.00078515,  0.00077967,  0.00079993,  0.00085375,  0.00094349,
      0.00106717,  0.00121825,  0.00138438,  0.00154968,  0.00170171,
      0.0018332 ,  0.00194182,  0.00202949,  0.00210058,  0.00215954,
      0.00221083,  0.00225905,  0.00230878,  0.00236425,  0.00242902,
      0.00250614,  0.00259834,  0.00270792,  0.00283638,  0.00298377,
      0.00314906,  0.00333064,  0.00352692,  0.00373758,  0.00396539,
      0.00421567,  0.0044955 ,  0.00481308,  0.00517634,  0.00559085,
      0.00606009,  0.00658595,  0.00716878,  0.00780775,  0.00850146,
      0.00924804,  0.01004529,  0.01089158,  0.01178793,  0.01274115,
      0.0137633 ,  0.01487031,  0.01608018,  0.01740874,  0.01886325,
      0.02044349,  0.02214463,  0.02396039,  0.02589065,  0.0279525 ,
      0.03017836,  0.03261135,  0.03530052,  0.03828981,  0.04160153,
      0.04523777,  0.04918468,  0.05341633,  0.05790466,  0.06263516,
      0.06760523,  0.07281963,  0.07828758,  0.08401736,  0.09000903,
      0.09625542,  0.10274424,  0.10945923,  0.11638187,  0.1234935 ,
      0.13077522,  0.13820759,  0.14577067,  0.15344416,  0.16120755,
      0.16904026,  0.17692176,  0.18483165,  0.19274966,  0.20065553,
      0.20852876,  0.2163489 ,  0.22409584,  0.23174999,  0.23929245,
      0.2467051 ])

for region in dismod3.gbd_regions:
    for year in dismod3.gbd_years:
        for sex in dismod3.gbd_sexes:
            key = dismod3.gbd_key_for('%s', region, year, sex)

            if clean(region) == 'north_america_high_income':
                regional_offset = 0.
            else:
                regional_offset = -.5

            time_offset = (int(year)-1997)/10.

            if clean(sex) == 'male':
                sex_offset = .25
            else:
                sex_offset = 0.
            
            # incidence rate
            i = mc.invlogit(mc.logit(.012 * mc.invlogit((ages - 44) / 3)) + regional_offset + time_offset + sex_offset)
            truth[key % 'incidence'] = i

            # remission rate
            r = 0. * ages
            truth[key % 'remission'] = r

            # excess-mortality rate
            f = .085 * (ages / 100) ** 2.5
            truth[key % 'excess-mortality'] = f

            ## compartmental model (bins S, C, D, M)
            SCDM = np.zeros([4, age_len])
            SCDM[0,0] = 1.

            for a in range(age_len - 1):
                A = [[-i[a]-m[a],  r[a]          , 0., 0.],
                     [ i[a]     , -r[a]-m[a]-f[a], 0., 0.],
                     [      m[a],       m[a]     , 0., 0.],
                     [        0.,            f[a], 0., 0.]]

                SCDM[:,a+1] = trim(np.dot(scipy.linalg.expm2(A), SCDM[:,a]), 0, 1)

            S = SCDM[0,:]
            C = SCDM[1,:]

            # prevalence = # with condition / (# with condition + # without)
            p = C / (S + C + NEARLY_ZERO)
            truth[key % 'prevalence'] = p
            truth[key % 'relative-risk'] = (m + f) / m

            # duration = E[time in bin C]
            pr_exit = 1 - r - m - f
            X = np.empty(len(pr_exit))
            t = 1.
            for a in xrange(len(X) - 1, -1, -1):
                X[a] = t * pr_exit[a]
                t = 1 + X[a]
            truth[key % 'duration'] = X


# generate synthetic data from fictitious ground truth, using data age
# structure from selected disease model
print '\nsimulating noisy realizations'

dispersion = float(options.dispersion)
n = float(options.study_size)

def generate_synthetic_data(truth, key, d):
    """ create simulated data"""
    a0 = d['age_start']
    a1 = d['age_end']
    age_weights = d['age_weights']
        
    d.update(condition='type_2_diabetes',
             year_start=y,
             year_end=y)

    p0 = dismod3.utils.rate_for_range(truth[key], range(a0, a1 + 1), np.ones(a1 + 1 - a0)/(a1+1-a0))
    p0 = dismod3.utils.trim(p0, 1.e-6, 1. - 1.e-6)

    # TODO: make beta dispersion study level (instead of datum level)
    # p1 = mc.rbeta(p0 * dispersion, (1 - p0) * dispersion)
    p1 = p0

    # TODO: add additional covariates
    if key.find('prevalence') != -1:
        if random.random() < .5:
            d['self-reported'] = True
            p1 = mc.invlogit(mc.logit(p1) - .5)
        else:
            d['self-reported'] = False
    
    p2 = mc.rbinomial(n, p1) / n
    
    d['value'] = p2
    if p2 > 0:
        d['standard_error'] = np.sqrt(p2 * (1 - p2) / n)

    return d

data = []

for d in dm.data:
    t = d['data_type']
    r = d['gbd_region']
    y = d['year_start']
    s = d['sex']
    key = dismod3.gbd_key_for(t.replace(' data',''), r, y, s)
    if key not in truth.keys():
        continue

    data += [generate_synthetic_data(truth, key, d)]

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
        'Self Reported': d.get('self-reported', '')
        }
    return c

f_file = open('simulated_data.tsv', 'w')
csv_f = csv.writer(f_file, dialect=csv.excel_tab)

col_names = sorted(data_dict_for_csv(data[0]).keys())
        
csv_f.writerow(col_names)
for d in data:
    dd = data_dict_for_csv(d)
    csv_f.writerow([dd[c] for c in col_names])

f_file.close()

# upload a new disease model which knows ground truth (but needs to
# have the data from the csv loaded separately)

dm.data = []
dm.params.pop('id')
dm.id = -1
for key in truth:
    dm.set_truth(key, truth[key])

url = dismod3.post_disease_model(dm)
print url

