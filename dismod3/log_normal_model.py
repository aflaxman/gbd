import numpy as np
import pymc as mc

from dismod3.utils import trim, interpolate, rate_for_range, indices_for_range, generate_prior_potentials
from dismod3.settings import NEARLY_ZERO, MISSING

def setup(dm, key, data_list, rate_stoch):
    """ Generate the PyMC variables for a log-normal model of
    a function of age

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    key : str
      the name of the key for everything about this model (priors,
      initial values, estimations)

    data_list : list of data dicts
      the observed data to use in the beta-binomial liklihood function

    rate_stoch : pymc.Stochastic
      a PyMC stochastic (or deterministic) object, with
      len(rate_stoch.value) == len(dm.get_estimation_age_mesh()).

    Results
    -------
    vars : dict
      Return a dictionary of all the relevant PyMC objects for the
      log-normal model.  vars['rate_stoch'] is of particular
      relevance, for details see the beta_binomial_model
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    vars['rate_stoch'] = rate_stoch

    @mc.deterministic(name='%s_max' % key)
    def mu_max(mu=rate_stoch):
        return max(mu)
    @mc.deterministic(name='%s_min' % key)
    def mu_min(mu=rate_stoch):
        return min(mu)
    vars.update(mu_max=mu_max, mu_min=mu_min)
        
    # set up priors and observed data
    prior_str = dm.get_priors(key)
    vars['priors'] = generate_prior_potentials(prior_str, est_mesh, rate_stoch, mu_max, mu_min)

    vars['observed_rates'] = []
    for d in data_list:
        age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
        age_weights = d.get('age_weights', np.ones(len(age_indices)) / len(age_indices))

        lb, ub = dm.bounds_per_1(d)
        if np.isnan(lb) or lb == ub:
            se = 1.
        else:
            se = (np.log(ub) - np.log(lb)) / (2. * 1.96)
            
        @mc.observed
        @mc.stochastic(name='obs_%d' % d['id'])
        def obs(f=rate_stoch,
                age_indices=age_indices,
                age_weights=age_weights,
                value=np.log(dm.value_per_1(d)),
                tau=se**-2, data=d):
            f_i = rate_for_range(f, age_indices, age_weights)
            return mc.normal_like(value, np.log(f_i), tau)
        vars['observed_rates'].append(obs)
        
    return vars
