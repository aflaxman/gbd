""" Test Age Pattern Model

These tests are use randomized computation, so they might fail
occasionally due to stochastic variation
"""

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc

import data
import rate_model
import age_pattern
reload(age_pattern)

def test_age_pattern_model_sim():
    # simulate normal data
    a = pl.arange(0, 100, 5)
    pi_true = .0001 * (a * (100. - a) + 100.)
    sigma_true = .025

    p = pl.maximum(0., mc.rnormal(pi_true, 1./sigma_true**2.))

    # create model and priors
    vars = {}

    vars.update(age_pattern.age_pattern('test', ages=pl.arange(101), knots=pl.arange(0,101,5), sigma=.1))

    vars['pi'] = mc.Lambda('pi', lambda mu=vars['mu_age'], a=a: mu[a])
    vars.update(rate_model.normal_model('test', vars['pi'], 0., p, sigma_true))

    # fit model
    m = mc.MCMC(vars)
    m.sample(2)


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
