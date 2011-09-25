""" Age integrating models"""

import pylab as pl
import pymc as mc


def midpoint_approx(name, mu_age, age_start, age_end):
    """ Generate PyMC objects for approximating the integral of gamma from age_start[i] to age_end[i]

    Parameters
    ----------
    name : str
    mu_age : pymc.Node with values of PCGP
    age_start, age_end : array

    Results
    -------
    Returns dict of PyMC objects, including 'mu_interval'
    the approximate integral of gamma  data predicted stochastic
    """
    @mc.deterministic(name='mu_interval_%s'%name)
    def mu_interval(mu_age=mu_age, age_start=age_start, age_end=age_end):
        age_mid = (age_start + age_end) / 2
        return mu_age.take(age_mid)

    return dict(mu_interval=mu_interval)
