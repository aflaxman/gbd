""" Functions for generating synthetic data

Data will have the form::

    region,country,year,age,y,se,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9
    Caribbean,ABW,1980,0,40.65,0.0,1,0,0.6,0.6,1.5,-0.5,0.5,0.8,1.0,0.5
    ...
"""

import pylab as pl
import pymc as mc
import pymc.gp as gp
import csv

# full data goal
age_range = pl.arange(0,81,5)
time_range = pl.arange(1980, 2005)
regions = 21

def generate_fe(out_fname='data.csv'):
    """ Generate random data based on a fixed effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = beta * X_r,c,t

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    a = 20.
    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    data = col_names()
    for t in time_range:
        for r in c4:
            for c in c4[r]:
                x = [1] + [t-1990.] + list(pl.randn(8))
                y = float(pl.dot(beta, x))
                se = 0.
                data.append([r, c, t, a, y, se] + list(x))

    write(data, out_fname)

def generate_smooth_gp_re_a(out_fname='data.csv', country_variation=True):
    """ Generate random data based on a nested gaussian process random
    effects model with age, with covariates that vary smoothly over
    time (where unexplained variation in time does not interact with
    unexplained variation in age)

    This function generates data for all countries in all regions, and
    all age groups based on the model::

        Y_r,c,t = beta * X_r,c,t + f_r(t) + g_r(a) + f_c(t)

        beta = [30., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        f_r ~ GP(0, C(3.))
        g_r ~ GP(0, C(2.))
        f_c ~ GP(0, C(1.)) or 0 depending on country_variation flag
        C(amp) = Matern(amp, scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ GP(t; 0, C(1)) for k >= 2
    """
    c4 = countries_by_region()

    data = col_names()

    beta = [30., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    C0 = gp.matern.euclidean(time_range, time_range, amp=3., scale=10., diff_degree=2)
    C1 = gp.matern.euclidean(age_range, age_range, amp=3., scale=10., diff_degree=2)
    C2 = gp.matern.euclidean(time_range, time_range, amp=.1, scale=5., diff_degree=2)
    C3 = gp.matern.euclidean(time_range, time_range, amp=1., scale=10., diff_degree=2)

    g = mc.rmv_normal_cov(pl.zeros_like(age_range), C1)
    for r in c4:
        f_r = mc.rmv_normal_cov(pl.zeros_like(time_range), C0)
        g_r = mc.rmv_normal_cov(g, C1)
        for c in c4[r]:
            f_c = mc.rmv_normal_cov(pl.zeros_like(time_range), C2)

            x_gp = {}
            for k in range(2,10):
                x_gp[k] = mc.rmv_normal_cov(pl.zeros_like(time_range), C3)

            for j, t in enumerate(time_range):
                for i, a in enumerate(age_range):
                    x = [1] + [j] + [x_gp[k][j] for k in range(2,10)]
                    y = float(pl.dot(beta, x)) + f_r[j] + g_r[i]
                    if country_variation:
                        y += f_c[j]
                    se = 0.
                    data.append([r, c, t, a, y, se] + list(x))
    write(data, out_fname)

def add_sampling_error(in_fname='data.csv', out_fname='noisy_data.csv', std=1.):
    """ add normally distributed noise to data.csv y column

    Parameters
    ----------
    std : float, or array of floats
      standard deviation of noise
    """
    data = pl.csv2rec(in_fname)
    if type(std) == float:
        std = std * pl.ones(len(data))
    for i, row in enumerate(data):
        data[i].y += std[i] * pl.randn(1)
        data[i].se += std[i]
    pl.rec2csv(data, out_fname)

def knockout_uniformly_at_random(in_fname='noisy_data.csv', out_fname='missing_noisy_data.csv', pct=20.):
    """ replace data.csv y column with uniformly random missing entries

    Parameters
    ----------
    pct : float, percent to knockout
    """
    data = pl.csv2rec(in_fname)
    for i, row in enumerate(data):
        if pl.rand() < pct/100.:
            data[i].y = pl.nan
    pl.rec2csv(data, out_fname)

# helper functions
def write(data, out_fname):
    """ write data to file"""
    fout = open(out_fname, 'w')
    csv.writer(fout).writerows(data)
    fout.close()

def countries_by_region():
    """ form dictionary of countries, keyed by gbd region"""
    c4 = dict([[d[0], d[1:]] for d in csv.reader(open('../country_region.csv'))])
    c4.pop('World')

    [c4.pop(k) for k in sorted(c4.keys())[regions:]]  # keep only specified number of regions
    return c4

def col_names():
    """ generate column names for csv file"""
    return [['region', 'country', 'year', 'age', 'y', 'se'] + ['x%d'%i for i in range(10)]]
