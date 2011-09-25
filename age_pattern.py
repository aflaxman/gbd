""" Several rate models"""

import pylab as pl
import pymc as mc


def pcgp(name, gamma_bar, knots, rho):
    """ Generate PyMC objects for a piecewise constant Gaussian process (PCGP) model

    Parameters
    ----------
    name : str
    gamma_bar : pymc.Node, expected values of PCGP
    knots : array, locations of the discontinuities in the piecewise constant function
    rho : pymc.Node, smoothness parameter for Matern covariance of GP

    Results
    -------
    Returns dict of PyMC objects, including 'gamma' and 'mu'
    the observed stochastic likelihood and data predicted stochastic
    """

    C_func = mc.gp.FullRankCovariance(mc.gp.matern.euclidean, amp=1., scale=rho, diff_degree=2)
    C_chol = pl.cholesky(C_func(knots, knots))
    gamma = mc.MvNormalChol('gamma_%s'%name, mu=pl.zeros_like(knots), sig=C_chol, value=pl.zeros_like(knots))

    import scipy.interpolate
    all_ages = pl.arange(100)
    @mc.deterministic(name='mu_%s'%name)
    def mu(gamma_bar=gamma_bar, gamma=gamma, knots=knots):
        mu = scipy.interpolate.interp1d(knots, pl.exp(gamma_bar + gamma), 'zero', bounds_error=False, fill_value=0.)
        return mu(all_ages)

    return dict(gamma=gamma, mu=mu)
