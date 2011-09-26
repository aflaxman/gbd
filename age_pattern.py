""" Several rate models"""

import pylab as pl
import pymc as mc


def pcgp(name, knots, ages, rho):
    """ Generate PyMC objects for a piecewise constant Gaussian process (PCGP) model

    Parameters
    ----------
    name : str
    knots : array, locations of the discontinuities in the piecewise constant function
    ages : array, points to interpolate to
    rho : pymc.Node, smoothness parameter for Matern covariance of GP

    Results
    -------
    Returns dict of PyMC objects, including 'gamma' and 'mu_age'
    the observed stochastic likelihood and data predicted stochastic
    """
    gamma_bar = mc.Uniform('gamma_bar_%s'%name, -20., 20., value=-5.)

    # TODO: trim the edges of the knots based on the level value
    # parameters, so we don't waste time optimizing something that
    # doesn't matter

    C_func = mc.gp.FullRankCovariance(mc.gp.matern.euclidean, amp=1., scale=rho, diff_degree=2)
    C_chol = pl.cholesky(C_func(knots, knots))
    gamma = mc.MvNormalChol('gamma_%s'%name, mu=pl.zeros_like(knots), sig=C_chol, value=pl.zeros_like(knots))

    import scipy.interpolate
    @mc.deterministic(name='mu_age_%s'%name)
    def mu_age(gamma_bar=gamma_bar, gamma=gamma, knots=knots, ages=ages):
        mu = scipy.interpolate.interp1d(knots, pl.exp(gamma_bar + gamma), 'zero', bounds_error=False, fill_value=0.)
        return mu(ages)

    return dict(gamma_bar=gamma_bar, gamma=gamma, mu_age=mu_age)
