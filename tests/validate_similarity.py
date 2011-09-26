""" Test similarity prior

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

    vars.update(age_pattern.pcgp('test', ages=pl.arange(101), knots=pl.arange(0,101,5), rho=25.))

    vars['pi'] = mc.Lambda('pi', lambda mu=vars['mu_age'], a=a: mu[a])
    vars.update(similarity_prior_model.similar('test', vars['pi'], pi_parent, w))

    # fit model
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [m.gamma_bar, m.gamma])
    m.sample(20000, 10000, 10)

    print 'pi mc error:', m.pi.stats()['mc error'].round(2)

    # plot results
    pl.plot(pl.arange(101), m.mu_age.stats()['95% HPD interval'], 'k', linestyle='steps-post:')
    pl.plot(pl.arange(101), m.mu_age.stats()['mean'], 'k-', drawstyle='steps-post')
    pl.plot(a, pi_parent, 'g-')

    # compare estimate to ground truth (skip endpoints, because they are extra hard to get right)
    assert pl.allclose(m.pi.stats()['mean'], pi_parent, atol=.05)
    lb, ub = m.pi.stats()['95% HPD interval'].T
    assert pl.mean((lb <= pi_parent) & (pi_parent <= ub)) > .75


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
