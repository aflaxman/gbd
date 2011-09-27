""" Several rate models"""

import pylab as pl
import pymc as mc


def pcgp(name, ages, knots, rho):
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
    gamma = mc.Uniform('gamma_%s'%name, -12., 6., value=pl.zeros_like(knots))

    import scipy.interpolate
    @mc.deterministic(name='mu_age_%s'%name)
    def mu_age(gamma_bar=gamma_bar, gamma=gamma, knots=knots, ages=ages):
        mu = scipy.interpolate.interp1d(knots, pl.exp(gamma_bar + gamma), 'zero', bounds_error=False, fill_value=0.)
        return mu(ages)

    C_func = mc.gp.FullRankCovariance(mc.gp.matern.euclidean, amp=10., scale=rho, diff_degree=2)
    C_chol = pl.cholesky(C_func(ages, ages))
    @mc.potential(name='smooth_mu_%s'%name)
    def smooth_gamma(mu_age=mu_age, gamma_bar=gamma_bar, C_chol=C_chol):
        return mc.mv_normal_chol_like(pl.log(mu_age).clip(-12., 6.) - gamma_bar, pl.zeros_like(mu_age), C_chol)

    return dict(gamma_bar=gamma_bar, gamma=gamma, mu_age=mu_age, smooth_gamma=smooth_gamma, ages=ages)
