""" Test covariate model

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
import covariate_model
reload(covariate_model)

def test_covariate_model_sim():
    # simulate normal data
    X = mc.rnormal(0., 1.**2, size=(128,3))
    beta_true = [-.1, .1, .2]
    Y_true = pl.dot(X, beta_true)

    pi_true = pl.exp(Y_true)
    sigma_true = .1

    p = mc.rnormal(pi_true, 1./sigma_true**2.)

    # create model and priors
    vars = dict(beta=mc.Uninformative('beta', value=pl.zeros_like(beta_true)))
    vars.update(covariate_model.covariate_model('test', 1, X, vars['beta']))
    vars.update(rate_model.normal_model('test', vars['pi'], 0., p, sigma_true))

    # fit model
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [m.beta])
    m.sample(20000, 10000, 10)

    # compare estimate to ground truth (skip endpoints, because they are extra hard to get right)
    assert pl.allclose(m.beta.stats()['mean'], beta_true, rtol=.2)
    lb, ub = m.beta.stats()['95% HPD interval'].T
    assert pl.mean((lb <= beta_true) & (beta_true <= ub)) > .5


# TODO: test similarity priors

# TODO: test predicting age pattern

# TODO: test covariates for predicting heterogeneity

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
