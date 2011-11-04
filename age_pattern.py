""" Several rate models"""

import pylab as pl
import pymc as mc


def pcgp(name, ages, knots, sigma):
    """ Generate PyMC objects for a piecewise constant Gaussian process (PCGP) model

    Parameters
    ----------
    name : str
    knots : array, locations of the discontinuities in the piecewise constant function
    ages : array, points to interpolate to
    sigma : pymc.Node, smoothness parameter for Matern covariance of GP

    Results
    -------
    Returns dict of PyMC objects, including 'gamma' and 'mu_age'
    the observed stochastic likelihood and data predicted stochastic
    """
    gamma = [mc.Normal('gamma_%s_%d'%(name,k), 0., 10.**-2, value=0) for k in knots]
    #gamma = [mc.Uniform('gamma_%s_%d'%(name,k), -20., 20., value=-10.) for k in knots]

    # TODO: fix AdaptiveMetropolis so that this is not necessary
    flat_gamma = mc.Lambda('flat_gamma_%s'%name, lambda gamma=gamma: pl.array([x for x in pl.flatten(gamma)]))


    import scipy.interpolate
    @mc.deterministic(name='mu_age_%s'%name)
    def mu_age(gamma=flat_gamma, knots=knots, ages=ages):
        mu = scipy.interpolate.interp1d(knots, pl.exp(gamma), 'linear', bounds_error=False, fill_value=0.)
        return mu(ages)

    vars = dict(gamma=gamma, mu_age=mu_age, ages=ages, knots=knots)

    if (sigma > 0) and (not pl.isinf(sigma)):
        print 'adding smoothing of', sigma
        @mc.potential(name='smooth_mu_%s'%name)
        def smooth_gamma(gamma=flat_gamma, knots=knots, tau=sigma**-2):
            # the following is to include a "noise floor" so that level value
            # zero prior does not exert undue influence on age pattern
            # smoothing
            gamma = gamma.clip(pl.log(pl.exp(gamma).mean()/10.), pl.inf)  # only include smoothing on values within 10x of mean

            return mc.normal_like(pl.sum(pl.diff(gamma)**2) / (knots[-1] - knots[0]), 0, tau)
        vars['smooth_gamma'] = smooth_gamma

    return vars
