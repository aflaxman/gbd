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
    data_type = 'p'
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    d = data.ModelData()
    d.input_data = data_simulation.simulated_age_intervals(data_type, n, a, pi_age_true, sigma_true)
    d.hierarchy, d.output_template = data_simulation.small_output()
    
    # create model and priors
    vars = data_model.data_model('test', d, data_type,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None)


    # fit model
    m = mc.MCMC(vars)
    m.sample(3)

    # check estimates
    pi_usa = covariate_model.predict_for(d.output_template, d.hierarchy, 'all', 'total', 'all', 'USA', 'male', 1990, vars)


    # create model w/ emp prior
    # create model and priors
    vars = data_model.data_model('test', d, data_type,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=pi_usa.mean(0))

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
