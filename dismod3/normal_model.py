import numpy as np
import pymc as mc
import random

from dismod3.utils import trim, interpolate, rate_for_range, indices_for_range, generate_prior_potentials
from dismod3.settings import NEARLY_ZERO, MISSING

def setup(dm, key, data_list, rate_stoch):
    """ Generate the PyMC variables for a normal model of
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
      normal model.  vars['rate_stoch'] is of particular
      relevance, for details see the beta_binomial_model
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    if np.any(np.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    vars['rate_stoch'] = rate_stoch

    # set up priors and observed data
    prior_str = dm.get_priors(key)
    vars['priors'] = generate_prior_potentials(prior_str, est_mesh, rate_stoch)

    vars['observed_rates'] = []
    for d in data_list:
        # set up observed stochs for all relevant data
        id = d['id']
        
        if d['value'] == MISSING:
            print 'WARNING: data %d missing value' % id
            continue

        # ensure all rate data is valid
        d_val = dm.value_per_1(d)
        d_se = dm.se_per_1(d)

        if d['age_start'] < est_mesh[0] or d['age_end'] > est_mesh[-1]:
            raise ValueError, 'Data %d is outside of estimation range---([%d, %d] is not inside [%d, %d])' \
                % (d['id'], d['age_start'], d['age_end'], est_mesh[0], est_mesh[-1])

        age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
        age_weights = d.get('age_weights', np.ones(len(age_indices)) / len(age_indices))

        # data must have standard error to use normal model
        if d_se == 0:
            raise ValueError, 'Data %d has invalid standard error' % d['id']

        @mc.observed
        @mc.stochastic(name='obs_%d' % id)
        def obs(f=rate_stoch,
                age_indices=age_indices,
                age_weights=age_weights,
                value=d_val,
                tau=1./(d_se)**2):
            f_i = rate_for_range(f, age_indices, age_weights)
            return mc.normal_like(value, f_i, tau)
        vars['observed_rates'].append(obs)
        
    return vars
