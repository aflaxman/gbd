""" Expert prior models"""

import pylab as pl
import pymc as mc

import similarity_prior_model

def level_constraints(name, parameters, unconstrained_mu_age, ages):
    """ Generate PyMC objects implementing priors on the value of the rate function

    Parameters
    ----------
    name : str
    parameters : dict
    unconstrained_mu_age : pymc.Node with values of PCGP

    Results
    -------
    Returns dict of PyMC objects, including 'unconstrained_mu_age' and 'mu_age'
    """
    if 'level_value' not in parameters or 'level_bounds' not in parameters:
        return {}

    @mc.deterministic(name='value_constrained_mu_age_%s'%name)
    def mu_age(unconstrained_mu_age=unconstrained_mu_age,
               value=parameters['level_value']['value'],
               age_before=pl.clip(parameters['level_value']['age_before']-ages[0], 0, len(ages)),
               age_after=pl.clip(parameters['level_value']['age_after']-ages[0], 0, len(ages)),
               lower=parameters['level_bounds']['lower'],
               upper=parameters['level_bounds']['upper']):
        mu_age = unconstrained_mu_age.copy()
        mu_age[:age_before] = value
        if age_after < len(mu_age)-1:
            mu_age[(age_after+1):] = value
        return mu_age.clip(lower, upper)

    mu_sim = similarity_prior_model.similar('value_constrained_mu_age_%s'%name, mu_age, unconstrained_mu_age, .01)

    return dict(mu_age=mu_age, unconstrained_mu_age=unconstrained_mu_age, mu_sim=mu_sim)


def derivative_constraints(name, parameters, mu_age, ages):
    """ Generate PyMC objects implementing priors on the value of the rate function

    Parameters
    ----------
    name : str
    parameters : dict
    mu_age : pymc.Node with values of PCGP
    ages : array

    Results
    -------
    Returns dict of PyMC objects, including 'mu_age_derivative_potential'
    """
    if 'increasing' not in parameters or 'decreasing' not in parameters:
        return {}

    @mc.potential(name='mu_age_derivative_potential_%s'%name)
    def mu_age_derivative_potential(mu_age=mu_age,
                                    increasing_a0=pl.clip(parameters['increasing']['age_start']-ages[0], 0, len(ages)),
                                    increasing_a1=pl.clip(parameters['increasing']['age_end']-ages[0], 0, len(ages)),
                                    decreasing_a0=pl.clip(parameters['decreasing']['age_start']-ages[0], 0, len(ages)),
                                    decreasing_a1=pl.clip(parameters['decreasing']['age_end']-ages[0], 0, len(ages))):
        mu_prime = pl.diff(mu_age)
        inc_violation = mu_prime[increasing_a0:increasing_a1].clip(-1., 0.).sum()
        dec_violation = mu_prime[decreasing_a0:decreasing_a1].clip(0., 1.).sum()
        return mc.normal_like([inc_violation, dec_violation], 0., 1.e-6**-2)

    return dict(mu_age_derivative_potential=mu_age_derivative_potential)

