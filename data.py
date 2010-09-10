""" Functions for generating synthetic data

Data will have the form::

    region,country,year,age,y,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9
    Oceania,FJI,1990,35,11.4,1,0.0,-0.08,-0.4,0.1,0.2,-0.8,-1.1,-1.8,1.0
    ...
"""

from pylab import randn, dot, arange, zeros, zeros_like, array
from pymc import rmv_normal_cov, gp
import csv

# full data goal
age_range = arange(0,81,5)
time_range = arange(1980, 2005)
regions=21

def write(data, out_fname):
    """ write data to file"""
    fout = open(out_fname, 'w')
    csv.writer(fout).writerows(data)
    fout.close()

def countries_by_region():
    """ form dictionary of countries, keyed by gbd region"""
    c4 = dict([[d[0], d[1:]] for d in csv.reader(open('../country_region.csv'))])
    c4.pop('World')

    [c4.pop(k) for k in sorted(c4.keys())[regions:]]  # remove more keys for faster testing
    return c4

def col_names(comment=''):
    """ generate column names for csv file"""
    if comment:
        header = [ ['# %s' % comment] ]
    else:
        header = []
    return header + [['region', 'country', 'year', 'age', 'y'] + ['x%d'%i for i in range(10)]]

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
    data = col_names(comment='fe model, beta=%s' % beta)
    for t in time_range:
        for r in c4:
            for c in c4[r]:
                x = [1] + [t-1990.] + list(randn(8))
                y = float(dot(beta, x))
                data.append([r, c, t, a, y] + list(x))

    write(data, out_fname)


def generate_re(out_fname='data.csv'):
    """ Generate random data based on a random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = (beta + u_r,c,t) * X_r,c,t

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        u_r,c,t[k] ~ N(0,.1^2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    a = 20.
    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    data = col_names(comment='re model, beta=%s'%beta)
    for t in time_range:
        for r in c4:
            for c in c4[r]:
                x = [1] + [t-1990.] + list(randn(8))
                y = float(dot(beta+.1*randn(10), x))
                data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)


def generate_gp_re(out_fname='data.csv'):
    """ Generate random data based on a gaussian process random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = beta * X_r,c,t + f_r(t)

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        f_r ~ GP(0,C)
        C = Matern(amp=3., scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    a = 20.
    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    data = col_names('gp re model')
    C = gp.matern.euclidean(time_range, time_range, amp=3., scale=20., diff_degree=2)
    for r in c4:
        f = rmv_normal_cov(zeros(15), C)
        for c in c4[r]:
            for t in time_range:
                x = [1] + [t-1990.] + list(randn(8))
                y = float(dot(beta, x)) + f[t-1990]
                data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)


def generate_nre(out_fname='data.csv'):
    """ Generate random data based on a nested random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        u_r[k] ~ N(0,.1^2)
        u_r,c,t[k] ~ N(0,.1^2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    a = 20.
    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    data = col_names('nested re model')
    for t in time_range:
        for r in c4:
            u_r = .1*randn(10)
            for c in c4[r]:
                x = [1] + [t-1990.] + list(randn(8))
                y = float(dot(beta + u_r + .1*randn(10), x))
                data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)

def generate_ngp_re(out_fname='data.csv'):
    """ Generate random data based on a nested gaussian process random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = (beta + u_r) * X_r,c,t + f_r(t) + f_c(t)

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        u_r ~ N(0, 1.^2)
        f_r ~ GP(0, C1)
        f_c ~ GP(0,C2)
        C_r = Matern(amp=3., scale=20., diff_degree=2)
        C_c = Matern(amp=1., scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    a = 20.
    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    data = col_names(comment='nested gaussian process re model')
    
    C_r = gp.matern.euclidean(time_range, time_range, amp=3., scale=20., diff_degree=2)
    C_c = gp.matern.euclidean(time_range, time_range, amp=1., scale=20., diff_degree=2)
    for r in c4:
        f_r = rmv_normal_cov(zeros_like(time_range), C_r)
        for c in c4[r]:
            f_c = rmv_normal_cov(zeros_like(time_range), C_c)
            
            for t in range(1990, 2005):
                x = [1] + [t-1990.] + list(randn(8))
                y = float(dot(beta, x)) + f_r[t-1990] + f_c[t-1990]
                data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)

def generate_ngp_re_a(out_fname='data.csv'):
    """ Generate random data based on a nested gaussian process random
    effects model with age

    This function generates data for all countries in all regions, and
    all age groups based on the model::

        Y_r,c,t = (beta + u_r + u_a) * X_r,c,t + f_r(t, a) + f_c(t, a)

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        u_r ~ N(0, 1.^2)
        u_a ~ N(0, 1.^2)
        f_r ~ GP(0, C1)
        f_c ~ GP(0,C2)
        C_r = Matern(amp=3., scale=20., diff_degree=2)
        C_c = Matern(amp=1., scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ N(0, 1) for k >= 2
    """
    c4 = countries_by_region()

    data = col_names('nested gp re model with age')

    age_time = array([[a/5., t] for a in age_range for t in time_range])

    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    for i, a in enumerate(age_range):
        for r in c4:
            C_r = gp.matern.euclidean(age_time, age_time, amp=3., scale=20., diff_degree=2)
            f_r = rmv_normal_cov(zeros(len(age_time)), C_r)
            for c in c4[r]:
                C_c = gp.matern.euclidean(age_time, age_time, amp=1., scale=20., diff_degree=2)
                f_c = rmv_normal_cov(zeros(len(age_time)), C_c)

                for j, t in enumerate(time_range):
                    x = [1] + [t-1990.] + list(randn(8))
                    y = float(dot(beta, x)) + f_r[i*len(time_range) + j] + f_c[i*len(time_range) + j]
                    data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)


# current amount for development and testing
age_range = arange(0,81,40)
time_range = arange(1980, 2005, 5)
regions = 8
def generate_smooth_gp_re_a(out_fname='data.csv'):
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
        f_c ~ GP(0, C(1.))
        C(amp) = Matern(amp, scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ GP(t; 0, C(1)) for k >= 2
    """
    c4 = countries_by_region()

    data = col_names(comment='smooth gp random effect with age')

    beta = [30., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    C0 = gp.matern.euclidean(time_range, time_range, amp=3., scale=10., diff_degree=2)
    C1 = gp.matern.euclidean(age_range, age_range, amp=3., scale=10., diff_degree=2)
    C2 = gp.matern.euclidean(time_range, time_range, amp=.1, scale=5., diff_degree=2)
    C3 = gp.matern.euclidean(time_range, time_range, amp=1., scale=10., diff_degree=2)

    g = rmv_normal_cov(zeros_like(age_range), C1)
    for r in c4:
        f_r = rmv_normal_cov(zeros_like(time_range), C0)
        g_r = rmv_normal_cov(g, C1)
        for c in c4[r]:
            f_c = rmv_normal_cov(zeros_like(time_range), C2)

            x_gp = {}
            for k in range(2,10):
                x_gp[k] = rmv_normal_cov(zeros_like(time_range), C3)

            for j, t in enumerate(time_range):
                for i, a in enumerate(age_range):
                    x = [1] + [j] + [x_gp[k][j] for k in range(2,10)]
                    y = float(dot(beta, x)) + f_r[j] + g_r[i] #+ f_c[j]
                    data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)

def generate_smooth_ngp_re_a_full(out_fname='data.csv'):
    """ Generate random data based on a nested gaussian process random
    effects model with age, with covariates that vary smoothly over
    time

    This function generates data for all countries in all regions, and
    all age groups based on the model::

        Y_r,c,t = (beta + u_r + u_a) * X_r,c,t + f_r(t, a) + f_c(t, a)

        beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
        u_r ~ N(0, 1.^2)
        u_a ~ N(0, 1.^2)
        f_r ~ GP(0, C1)
        f_c ~ GP(0,C2)
        C_r = Matern(amp=3., scale=20., diff_degree=2)
        C_c = Matern(amp=1., scale=20., diff_degree=2)

        X_r,c,t[0] = 1
        X_r,c,t[1] = t - 1990.
        X_r,c,t[k] ~ GP(t; 0, C3) for k >= 2
        C3 = Matern(amp=1., scale=20., diff_degree=2)
    """
    c4 = countries_by_region()

    data = col_names()

    age_range = arange(0,81,10)
    time_range = arange(1990, 2005)
    age_time = array([[a/5., t] for a in age_range for t in time_range])

    beta = [10., -.5, .1, .1, -.1, 0., 0., 0., 0., 0.]
    C3 = gp.matern.euclidean(time_range, time_range, amp=1., scale=20., diff_degree=2)
    for i, a in enumerate(age_range):
        for r in c4:
            C_r = gp.matern.euclidean(age_time, age_time, amp=3., scale=20., diff_degree=2)
            f_r = rmv_normal_cov(zeros(len(age_time)), C_r)
            for c in c4[r]:
                C_c = gp.matern.euclidean(age_time, age_time, amp=1., scale=20., diff_degree=2)
                f_c = rmv_normal_cov(zeros(len(age_time)), C_c)

                x_gp = {}
                for k in range(2,10):
                    x_gp[k] = rmv_normal_cov(zeros(len(time_range)), C3)

                for j, t in enumerate(time_range):
                    x = [1] + [t-1990.] + [x_gp[k][t-1990.] for k in range(2,10)]
                    y = float(dot(beta, x)) + f_r[i*len(time_range) + j] + f_c[i*len(time_range) + j]
                    data.append([r, c, t, a, y] + list(x))
    write(data, out_fname)

def add_sampling_error(in_fname='data.csv', out_fname='noisy_data.csv', std=1.):
    """ add normally distributed noise to data.csv y column

    Parameters
    ----------
    std : float, standard deviation of noise
    """
    from pylab import csv2rec, rec2csv, randn

    data = csv2rec(in_fname, skiprows=1)  # skiprows hack, for this old version of csv2rec
    for i, row in enumerate(data):
        data[i].y += std * randn(1)
    rec2csv(data, out_fname)

def knockout_uniformly_at_random(in_fname='noisy_data.csv', out_fname='missing_noisy_data.csv', pct=20.):
    """ replace data.csv y column with uniformly random missing entries

    Parameters
    ----------
    pct : float, percent to knockout
    """
    from pylab import csv2rec, rec2csv, nan, rand

    data = csv2rec(in_fname)
    for i, row in enumerate(data):
        if rand() < pct/100.:
            data[i].y = nan
    rec2csv(data, out_fname)

