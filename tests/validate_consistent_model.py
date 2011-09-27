""" Test data model

These tests are use randomized computation, so they might fail
occasionally due to stochastic variation
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import consistent_model
reload(consistent_model)
reload(data_model)


def validate_consistent_model_sim():
    # generate simulated data
    n = 50
    sigma_true = .025

    # start with truth
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    # choose age intervals to measure
    age_start = pl.array(mc.runiform(0, 100, n), dtype=int)
    age_start.sort()  # sort to make it easy to discard the edges when testing
    age_end = pl.array(mc.runiform(age_start+1, pl.minimum(age_start+10,100)), dtype=int)

    # find truth for the integral across the age intervals
    import scipy.integrate
    pi_interval_true = [scipy.integrate.trapz(pi_age_true[a_0i:(a_1i+1)]) / (a_1i - a_0i) 
                        for a_0i, a_1i in zip(age_start, age_end)]

    # simulate the noisy measurement of the rate in each interval
    p = mc.rnormal(pi_interval_true, 1./sigma_true**2.)

    # store the simulated data in a pandas DataFrame
    data = pandas.DataFrame(dict(value=p, age_start=age_start, age_end=age_end))
    data['effective_sample_size'] = pl.maximum(p*(1-p)/sigma_true**2, 1.)
    data['standard_error'] = pl.nan
    data['year_start'] = 2005.  # TODO: make these vary
    data['year_end'] = 2005.
    data['sex'] = 'total'
    data['area'] = 'all'
    data['data_type'] = 'p'

    # generate a simple hierarchy graph for the model
    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')
    

    # create model and priors
    vars = consistent_model.consistent_model(data, {}, hierarchy, 'all')

    # fit model
    mc.MAP(vars).fit(method='fmin_powell', verbose=1, iterlim=30)
    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [vars[k]['gamma_bar'] for k in 'irf'] + [vars[k]['gamma'] for k in 'irf'])
    m.sample(30000, 15000, 15)

    return m

if __name__ == '__main__':
    m = validate_consistent_model_sim()
