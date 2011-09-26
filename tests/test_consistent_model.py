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

def test_consistent_model_forward():
    d = pandas.DataFrame(dict(data_type=[''], value=[0.], age_start=[0], age_end=[0], effective_sample_size=[1.],
                              standard_error=[pl.nan], upper_ci=[pl.nan], lower_ci=[pl.nan]))

    # generate a simple hierarchy graph for the model
    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')

    vars = consistent_model.consistent_model(d, hierarchy, 'all')

    vars['i']['gamma_bar'].value = pl.log(.01)
    vars['r']['gamma_bar'].value = pl.log(.0001)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)

    vars['i']['gamma_bar'].value = pl.log(.02)
    vars['r']['gamma_bar'].value = pl.log(.0001)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)

    vars['i']['gamma_bar'].value = pl.log(2.)
    vars['r']['gamma_bar'].value = pl.log(30.)
    vars['f']['gamma_bar'].value = pl.log(.0001)
    print vars['p']['mu_age'].value[::10].round(3)


def test_consistent_model_sim():
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
    data['data_type'] = 'p'
    data['effective_sample_size'] = pl.maximum(p*(1-p)/sigma_true**2, 1.)

    data = data.append(pandas.DataFrame(dict(value=[0.], age_start=[0], age_end=[100],
                                             data_type=['r'], effective_sample_size=[1000])), ignore_index=True)

    data['standard_error'] = pl.nan
    data['upper_ci'] = pl.nan
    data['lower_ci'] = pl.nan

    data['year_start'] = 2005.
    data['year_end'] = 2005.
    data['sex'] = 'total'
    data['area'] = 'all'


    # generate a simple hierarchy graph for the model
    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')
    

    # create model and priors
    vars = consistent_model.consistent_model(data, hierarchy, 'all')

    # fit model
    m = mc.MCMC(vars)
    m.sample(1)

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
