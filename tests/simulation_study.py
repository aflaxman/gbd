""" DisMod III Simulation Study

This file contains scripts for validating the Generic Disease Modeling
System DisMod III on simulated data loosely based on Type II Diabetes
in Southeast Asia.
"""

from pylab import *
import numpy as np
import pymc as mc
import simplejson as json

import dismod3
import dismod3.utils
from dismod3.disease_json import DiseaseJson

import dismod3.gbd_disease_model as model

# flags for which figures to generate
PLOT_TRUTH = False
PLOT_FIT = False

age_len = 101
ages = np.arange(age_len, dtype='float')

# define model transition parameters inspired by available data for
# Type II Diabetes in Southeast Asia

# incidence rate
i = .012 * mc.invlogit((ages - 44) / 3)

# remission rate
r = 0. * ages

# case-fatality rate
f = .085 * (ages / 100) ** 2.5

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


## compartmental model (bins S, C, D, M)
S = np.zeros(age_len); C = np.zeros(age_len); D = np.zeros(age_len); M = np.zeros(age_len)
S[0] = 1.; C[0] = 0.; D[0] = 0.; M[0] = 0.

for a in range(age_len - 1):
    S[a+1] = S[a]*(1-i[a]-m[a]) + C[a]*r[a]
    C[a+1] = S[a]*i[a]          + C[a]*(1-r[a]-m[a]-f[a])
    D[a+1] =                      C[a]*f[a]               + D[a]
    M[a+1] = S[a]*m[a]          + C[a]*m[a]                      + M[a]


# prevalence = # with condition / (# with condition + # without)
p = C / (S + C)

# duration = E[time in bin C]
pr_exit = 1 - r - m - f
X = np.empty(len(pr_exit))
t = 1.
for a in xrange(len(X) - 1, -1, -1):
    X[a] = t * pr_exit[a]
    t = 1 + X[a]


# make a plot of the ground truth for the simulation
if PLOT_TRUTH:
    clf()

    subplot(2,2,1)
    plot(ages, 100 * i, alpha=.75, linewidth=2, linestyle='-', label='Incidence')
    plot(ages, 100 * r, alpha=.75, linewidth=2, linestyle='--', label='Remission')
    plot(ages, 100 * f, alpha=.75, linewidth=2, linestyle=':', label='Case Fatality')
    plot(ages, 100 * m, '.-', alpha=.75, linewidth=2, label='All-cause Mortality')
    legend(loc='upper left')
    xlabel('Age (years)')
    ylabel('Rate (per 100)')
    title('Rates')
    axis([20, 80, -.1, 5])
    
    subplot(2,2,2)
    plot(ages, (m + f) / m, color='k', alpha=.75, linewidth=2, label='Relative Risk')
    xlabel('Age (years)')
    ylabel('Relative Risk')
    title('Relative Risk')
    axis([20, 80, 1, 3.5])
    
    subplot(2,2,3)
    plot(ages, 100*p, color='m', alpha=.75, linewidth=2, label='Prevalence')
    xlabel('Age (years)')
    ylabel('Prevalence (per 100)')
    title('Prevalence')
    axis([20, 80, -.5, 25])
    
    subplot(2,2,4)
    plot(X, color='r', alpha=.75, linewidth=2)
    xlabel('Age (years)')
    ylabel('Duration (years)')
    title('Duration')
    axis([20, 80, -1, 60])
    

# simulate noisy realizations of prevalence and case-fatality from ground truth

#dispersion = 500.
#n = 10000.

dispersion = float(sys.argv[1])
n = float(sys.argv[2])
    
data = []

prev_age_intervals = [[25, 34], [35, 44], [45, 54], [55, 64],
                      [65, 74], [75, 84], [85, 100]]

for a0, a1 in prev_age_intervals:
    d = { 'condition': 'type_2_diabetes',
          'data_type': 'prevalence data',
          'gbd_region': 'Asia, Southeast',
          'region': 'Thailand',
          'year_start': 2005,
          'year_end': 2005,
          'sex': 'male',
          'age_start': a0,
          'age_end': a1,
          'age_weights': list(np.ones(a1 + 1 - a0)),
          'id': len(data)}

    p0 = dismod3.utils.rate_for_range(p, range(a0, a1 + 1), np.ones(a1 + 1 - a0))
    p1 = mc.rbeta(p0 * dispersion, (1 - p0) * dispersion)
    p2 = mc.rbinomial(n, p1) / n
    
    d['value'] = p2
    d['standard_error'] = np.sqrt(p2 * (1 - p2) / n)

    data.append(d)

for a0, a1 in [[25,59], [60,74], [75,100]]:
    d = { 'condition': 'type_2_diabetes',
          'data_type': 'case-fatality data',
          'gbd_region': 'Asia, Southeast',
          'region': 'Thailand',
          'year_start': 2005,
          'year_end': 2005,
          'sex': 'male',
          'age_start': a0,
          'age_end': a1,
          'age_weights': list(np.ones(a1 + 1 - a0)),
          'id': len(data)}

    r0 = dismod3.utils.rate_for_range(f, range(a0, a1 + 1), np.ones(a1 + 1 - a0))
    r1 = mc.rbeta(r0 * dispersion, (1 - r0) * dispersion)
    r2 = mc.rbinomial(n, r1) / n
    
    d['value'] = r2
    d['standard_error'] = np.sqrt(r2 * (1 - r2) / n)

    data.append(d)

dm = DiseaseJson(json.dumps({'params':
                                 {"id": 1,
                                  "sex": "male",
                                  "region": "asia_southeast",
                                  "year": "2005",
                                  "condition": "type_ii_diabetes"},
                             "data": data}))

key = dismod3.utils.gbd_key_for('%s', dm.get_region(), 2005, 'male')
dm.set_initial_value(key % 'all-cause_mortality', m)

dm.set_priors(key % 'remission', ' zero 0 100, ')
dm.set_priors(key % 'case-fatality', ' zero 0 10, smooth 10, ')
dm.set_priors(key % 'incidence', ' smooth 10, ')

dm.set_truth(key % 'remission', r)
dm.set_truth(key % 'incidence', i)
dm.set_truth(key % 'prevalence', p)
dm.set_truth(key % 'case-fatality', f)
dm.set_truth(key % 'relative-risk', (m + f) / m)
dm.set_truth(key % 'duration', X)
dm.set_truth(key % 'yld', X * i)

keys = model.gbd_keys(region_list=['asia_southeast'], year_list=[2005], sex_list=['male'])
model.fit(dm, method='map', keys=keys)
model.fit(dm, method='mcmc', keys=keys)
model.fit(dm, method='map', keys=keys)

total_yld = sum(i * X)
est_yld = sum(dm.get_mcmc('median', key % 'yld'))
est_yld_upper_ui = sum(dm.get_mcmc('lower_ui', key % 'yld'))
est_yld_lower_ui = sum(dm.get_mcmc('upper_ui', key % 'yld'))

if PLOT_FIT:
    dismod3.plotting.tile_plot_disease_model(dm.to_json(), keys)
    figtext(.5, .15, total_yld_str = """
Dashed = Truth
Dotted = MLE
Solid = Median

True YLD = %.2f
Est YLD = %.2f (%.2f, %.2f)
""" % (total_yld, est_yld, est_yld_upper_ui, est_yld_lower_ui))


print n, dispersion, total_yld, est_yld, est_yld_lower_ui, est_yld_upper_ui
