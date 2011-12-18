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
    ages=pl.arange(101)

    # create model with no priors
    vars = {}
    vars.update(age_pattern.age_pattern('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(expert_prior_model.level_constraints('test', {}, vars['mu_age'], ages))

    # fit model
    m = mc.MCMC(vars)
    m.sample(3)


    # create model with expert priors
    parameters = {}
    parameters['level_value'] = dict(value=.1, age_below=15, age_above=95)
    parameters['level_bound'] = dict(upper=.01, lower=.001)
    vars = {}
    vars.update(age_pattern.age_pattern('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(expert_prior_model.level_constraints('test', parameters, vars['mu_age'], ages))

    # fit model
    m = mc.MCMC(vars)
    m.sample(3)


def test_expert_model_derivative_sign():
    d = data.ModelData()
    ages=pl.arange(101)

    # create model with no priors
    vars = {}
    vars.update(age_pattern.age_pattern('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(expert_prior_model.derivative_constraints('test', {}, vars['mu_age'], ages))

    # create model with expert priors
    parameters = {}
    parameters['increasing'] = dict(age_start=15, age_end=95)
    parameters['decreasing'] = dict(age_start=0, age_end=0)
    vars = {}
    vars.update(age_pattern.age_pattern('test', ages, knots=pl.arange(0,101,5), sigma=.01))
    vars.update(expert_prior_model.derivative_constraints('test', parameters, vars['mu_age'], vars['knots']))


    # fit model
    m = mc.MCMC(vars)
    m.sample(3)


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
