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

