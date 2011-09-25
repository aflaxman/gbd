""" Test Age Pattern Model

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
reload(age_pattern)

def test_age_pattern_model_sim():
    # simulate normal data
    a = pl.arange(0, 100, 5)
    pi_true = .0001 * (a * (100. - a) + 100.)
    sigma_true = .025

    p = mc.rnormal(pi_true, 1./sigma_true**2.)

    # create model and priors
    vars = {}

    vars.update(age_pattern.pcgp('test', knots=pl.arange(0,101,5), rho=25.))

    vars['pi'] = mc.Lambda('pi', lambda mu=vars['mu_age'], a=a: mu[a])
    vars.update(rate_model.normal_model('test', vars['pi'], 0., p, sigma_true))

    # fit model
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [m.gamma_bar, m.gamma])
    m.sample(20000, 10000, 10)

    # plot results
    pl.plot(pl.arange(100), m.mu_age.stats()['95% HPD interval'], 'k', linestyle='steps-post:')
    pl.plot(pl.arange(100), m.mu_age.stats()['mean'], 'k-', drawstyle='steps-post')
    pl.plot(a, pi_true, 'g-')
    pl.plot(a, p, 'ro')

    # compare estimate to ground truth (skip endpoints, because they are extra hard to get right)
    assert pl.allclose(m.pi.stats()['mean'][2:-2], pi_true[2:-2], rtol=.2)
    lb, ub = m.pi.stats()['95% HPD interval'].T
    assert pl.mean((lb <= pi_true)[2:-2] & (pi_true <= ub)[2:-2]) > .75


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
