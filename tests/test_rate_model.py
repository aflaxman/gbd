""" Test Rate Model

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
reload(rate_model)

def test_neg_binom_model_sim():
    # simulate negative binomial data
    pi_true = .01
    delta_true = 50

    n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
    k = pl.array(mc.rnegative_binomial(n*pi_true, delta_true), dtype=float)
    p = k/n

    # create NB model and priors
    vars = dict(pi=mc.Uniform('pi', 0., 1., value=.01),
                delta=mc.Uniform('delta', 0., 10000., value=1000.))
    vars.update(rate_model.neg_binom_model('sim', vars['pi'], vars['delta'], p, n))

    # fit NB model
    m = mc.MCMC(vars)
    m.sample(1)



def test_normal_model_sim():
    # simulate data
    pi_true = 10.
    sigma_true = 1.

    n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
    p = mc.rnormal(pi_true, 1./(sigma_true**2. + 1./n))

    # create model and priors
    vars = dict(pi=mc.Uniform('pi', 0., 1000., value=.01),
                sigma=mc.Uniform('sigma', 0., 10000., value=1000.))
    vars.update(rate_model.normal_model('sim', vars['pi'], vars['sigma'], p, 1./pl.sqrt(n)))

    # fit model
    m = mc.MCMC(vars)
    m.sample(1)


def test_log_normal_model_sim():
    # simulate negative binomial data
    pi_true = 2.
    sigma_true = .1

    n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
    p = pl.exp(mc.rnormal(pl.log(pi_true), 1./(sigma_true**2 + 1./n)))

    # create model and priors
    vars = dict(pi=mc.Uniform('pi', 0., 1000., value=.01),
                sigma=mc.Uniform('sigma', 0., 10000., value=1000.))
    vars.update(rate_model.log_normal_model('sim', vars['pi'], vars['sigma'], p, 1./pl.sqrt(n)))

    # fit model
    m = mc.MCMC(vars)
    m.sample(1)

if __name__ == '__main__':
    import nose
    nose.runmodule()
    
