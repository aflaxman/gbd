""" Test age integrating model

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

import data
import rate_model
import age_pattern
import age_integrating_model
reload(age_integrating_model)
reload(age_pattern)
import data_simulation

def test_age_standardizing_approx():
    # simulate data
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)
    ages=pl.arange(101)
    d = data_simulation.simulated_age_intervals(n, a, pi_age_true, sigma_true)

    # create model and priors
    vars = {}
    vars.update(age_pattern.pcgp('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(age_integrating_model.age_standardize_approx('test', pl.ones_like(vars['mu_age'].value), vars['mu_age'], d['age_start'], d['age_end'], ages))
    vars['pi'] = vars['mu_interval']
    vars.update(rate_model.normal_model('test', pi=vars['pi'], sigma=0, p=d['value'], s=sigma_true))

    # fit model
    m = mc.MCMC(vars)
    m.sample(3)

def test_age_integrating_midpoint_approx():
    # simulate data
    n = 50
    sigma_true = .025
    a = pl.arange(0, 100, 1)
    pi_age_true = .0001 * (a * (100. - a) + 100.)
    ages = pl.arange(101)
    d = data_simulation.simulated_age_intervals(n, a, pi_age_true, sigma_true)

    # create model and priors
    vars = {}
    vars.update(age_pattern.pcgp('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(age_integrating_model.midpoint_approx('test', vars['mu_age'], d['age_start'], d['age_end'], ages))
    vars['pi'] = vars['mu_interval']
    vars.update(rate_model.normal_model('test', pi=vars['pi'], sigma=0, p=d['value'], s=sigma_true))

    # fit model
    m = mc.MCMC(vars)
    m.sample(3)

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
