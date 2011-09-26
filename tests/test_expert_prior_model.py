""" Test expert prior model

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
import age_integrating_model
import expert_prior_model
reload(age_integrating_model)
reload(age_pattern)

def test_expert_model_level_value():
    d = data.ModelData()

    # create model with no priors
    vars = {}
    vars.update(age_pattern.pcgp('test', ages=pl.arange(101), knots=pl.arange(0,101,5), rho=40.))
    vars.update(expert_prior_model.level_constraints('test', d.parameters, vars['mu_age']))

    # create model with level value priors
    d.parameters['p']['level_value'] = dict(value=.1, age_below=15, age_above=95)
    vars = {}
    vars.update(age_pattern.pcgp('test', ages=pl.arange(101), knots=pl.arange(0,101,5), rho=40.))
    vars.update(expert_prior_model.level_constraints('test', d.parameters, vars['mu_age']))


    # fit model
    m = mc.MCMC(vars)
    m.sample(300)


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
