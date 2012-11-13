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
import ism
import covariate_model
import data_simulation
reload(ism)

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
    vars = ism.age_specific_rate(d, data_type,
                                 reference_area='all', reference_sex='total', reference_year='all',
                                 mu_age=None, mu_age_parent=None, sigma_age_parent=None)


    # fit model
    m = mc.MCMC(vars)
    m.sample(3)

    # check estimates
    pi_usa = covariate_model.predict_for(d, d.parameters, 'all', 'total', 'all', 'USA', 'male', 1990, 0., vars[data_type], -pl.inf, pl.inf)


    # create model w/ emp prior
    # create model and priors
    vars = ism.age_specific_rate(d, data_type,
                                 reference_area='all', reference_sex='total', reference_year='all',
                                 mu_age=None, mu_age_parent=pi_usa.mean(0), sigma_age_parent=pi_usa.std(0))


def test_data_model_lower_bound():
    # generate simulated data
    data_type = 'csmr'
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)

    d = data.ModelData()
    d.input_data = data_simulation.simulated_age_intervals(data_type, n, a, pi_age_true, sigma_true)
    d.input_data = d.input_data.append(data_simulation.simulated_age_intervals('pf', n, a, pi_age_true*2., sigma_true),
                                       ignore_index=True)
    d.hierarchy, d.output_template = data_simulation.small_output()
    
    # create model and priors
    vars = ism.age_specific_rate(d, 'pf',
                                 reference_area='all', reference_sex='total', reference_year='all',
                                 mu_age=None, mu_age_parent=None, sigma_age_parent=None, lower_bound='csmr')


    # fit model
    m = mc.MCMC(vars)
    m.sample(3)

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
