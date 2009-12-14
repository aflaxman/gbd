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

age_len = dismod3.MAX_AGE
ages = np.arange(age_len, dtype='float')

# incidence rate
i = .012 * mc.invlogit((ages - 44) / 3)

# remission rate
r = 0. * ages

# excess-mortality rate
f = .085 * (ages / 100) ** 2.5

# all-cause mortality-rate
m_all_cause = np.array(
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


def generate_and_append_data(data, data_type, truth, age_intervals,
                             gbd_region='North America, High Income', country='USA', year=2005, sex='male'):
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

        p0 = dismod3.utils.rate_for_range(truth, range(a0, a1 + 1), np.ones(a1 + 1 - a0))
    
        d['value'] = p0
        if p0 < 1.:
            d['standard_error'] = p0 * (1-p0) / np.sqrt(1000)
        else:
            d['standard_error'] = p0 * .05

        data.append(d)
    
data = []

age_intervals = [[a, a] for a in range(dismod3.MAX_AGE)]
year = 2005
sex = 'male'

generate_and_append_data(data, 'prevalence data', p, age_intervals, year=year, sex=sex)
generate_and_append_data(data, 'incidence data', i, age_intervals, year=year, sex=sex)
generate_and_append_data(data, 'relative-risk data', (m+f)/m, age_intervals, year=2005, sex=sex)

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
