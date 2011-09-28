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
import covariate_model
import data_simulation
reload(data_model)

def test_data_model_sim():
    # generate simulated data
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    data = data_simulation.simulated_age_intervals(n, a, pi_age_true, sigma_true)
    hierarchy, output_template = data_simulation.small_output()
    

    # create model and priors
    vars = data_model.data_model('test', data, {}, hierarchy, 'all')


    # fit model
    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [m.gamma_bar, m.gamma, m.beta])
    m.sample(3)

    # check estimates
    pi_usa = covariate_model.predict_for(output_template, hierarchy, 'all', 'USA', 'male', 1990, vars)

    

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
