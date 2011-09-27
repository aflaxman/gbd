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
import data_simulation

def test_consistent_model_forward():
    d = pandas.DataFrame(dict(data_type=[''], value=[0.], age_start=[0], age_end=[0], effective_sample_size=[1.],
                              standard_error=[pl.nan], upper_ci=[pl.nan], lower_ci=[pl.nan]))

    # generate a simple hierarchy graph for the model
    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')

    vars = consistent_model.consistent_model(d, {}, hierarchy, 'all')

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
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    data = data_simulation.simulated_age_intervals(n, a, pi_age_true, sigma_true)
    data['data_type'][-1] = 'r'  # make sure that there are multiple data types in the data set

    # generate simple hierarchy and priord for the model
    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')
    priors = {}

    # create model and priors
    vars = consistent_model.consistent_model(data, priors, hierarchy, 'all')

    # fit model
    m = mc.MCMC(vars)
    m.sample(1)

    return vars

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
