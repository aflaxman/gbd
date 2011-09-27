""" Age integrating models"""

import pylab as pl
import pymc as mc


def age_standardize_approx(name, age_weights, mu_age, age_start, age_end, ages):
    """ Generate PyMC objects for approximating the integral of gamma from age_start[i] to age_end[i]

    Parameters
    ----------
    name : str
    age_weights : array, len == len(ages)
    mu_age : pymc.Node with values of PCGP
    age_start, age_end : array

    Results
    -------
    Returns dict of PyMC objects, including 'mu_interval'
    the approximate integral of gamma  data predicted stochastic
    """
    cum_sum_weights = pl.cumsum(age_weights)

    @mc.deterministic(name='cum_sum_mu_%s'%name)
    def weighted_sum_mu(mu_age=mu_age, age_weights=age_weights):
        return pl.cumsum(mu_age*age_weights)

    age_start = age_start.__array__().clip(ages[0], ages[-1]) - ages[0]  # FIXME: Pandas bug, makes clip require __array__()
    age_end = age_end.__array__().clip(ages[0], ages[-1]) - ages[0]
    @mc.deterministic(name='mu_interval_%s'%name)
    def mu_interval(weighted_sum_mu=weighted_sum_mu, cum_sum_weights=cum_sum_weights, mu_age=mu_age, age_start=age_start, age_end=age_end):
        mu = (weighted_sum_mu[age_end] - weighted_sum_mu[age_start]) / (cum_sum_weights[age_end] - cum_sum_weights[age_start])
        
        # correct cases where age_start == age_end
        i = age_start == age_end
        if pl.any(i):
            mu[i] = mu_age[age_start[i]]

        return mu

    return dict(mu_interval=mu_interval)


def midpoint_approx(name, mu_age, age_start, age_end, ages):
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
        return mu_age.take(age_mid-ages[0])

    return dict(mu_interval=mu_interval)
