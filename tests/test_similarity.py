""" Test similarity prior

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
import similarity_prior_model
reload(similarity_prior_model)
reload(age_pattern)

def test_similarity_prior():
    # generate prior age pattern
    a = pl.arange(101)
    pi_parent = .0001 * (a * (100. - a) + 100.)

    w = .1

    # create model and priors
    vars = {}

    vars.update(age_pattern.age_pattern('test', ages=pl.arange(101), knots=pl.arange(0,101,5), sigma=.1))

    vars['pi'] = mc.Lambda('pi', lambda mu=vars['mu_age'], a=a: mu[a])
    vars.update(similarity_prior_model.similar('test', vars['pi'], pi_parent, 0., w))

    # fit model
    m = mc.MCMC(vars)
    m.sample(2)


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
