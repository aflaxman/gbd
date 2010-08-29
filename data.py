""" Functions for generating synthetic data"""

from pylab import randn, dot
import csv

def write(data):
    """ write data.csv file"""
    fout = open('data.csv', 'w')
    csv.writer(fout).writerows(data)
    fout.close()

def countries_by_region():
    c4 = dict([[d[0], d[1:]] for d in csv.reader(open('../country_region.csv'))])
    c4.pop('World')
    return c4

def col_names():
    return [['region', 'country', 'year', 'y'] + ['x%d'%i for i in range(10)]]

def generate_fe():
    """ replace data.csv with random data based on a fixed effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = beta . X_r,c,t + e_r,c,t
        e_r,c,t ~ N(0,1)

    """
    c4 = countries_by_region()

    data = col_names()
    beta = randn(10)
    for t in range(1990, 2005):
        for r in c4:
            for c in c4[r]:
                x = [1] + list(randn(9))
                y = float(dot(beta, x) + randn(1))
                data.append([r, c, t, y] + list(x))

    write(data)


def generate_re():
    """ replace data.csv with random data based on a random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = (beta + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r,c,t[k] ~ N(0,1)
        e_r,c,t ~ N(0,1)

        beta ~ N(0, 10^2)
        X_r,c,t[k] ~ N(0, 1) for k >= 1

    """
    c4 = countries_by_region()

    data = col_names()
    beta = 10.*randn(10)
    for t in range(1990, 2005):
        for r in c4:
            for c in c4[r]:
                x = [1] + list(randn(9))
                y = float(dot(beta+randn(10), x) + randn(1))
                data.append([r, c, t, y] + list(x))
    write(data)


def generate_nre():
    """ replace data.csv with random data based on a nested random effects model

    This function generates data for all countries in all regions, based on the model::

        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r[k] ~ N(0,2^2)
        u_r,c,t[k] ~ N(0,1)
        e_r,c,t ~ N(0,1)

        beta ~ N(0, 10^2)
        X_r,c,t[k] ~ N(0, 1) for k >= 1

    """
    c4 = countries_by_region()

    data = col_names()
    beta = 10.*randn(10)
    for t in range(1990, 2005):
        for r in c4:
            u_r = .2*randn(10)
            for c in c4[r]:
                x = [1] + list(randn(9))
                y = float(dot(beta + u_r + randn(10), x) + randn(1))
                data.append([r, c, t, y] + list(x))
    write(data)
