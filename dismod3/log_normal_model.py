import pylab as pl
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

    # set up priors and observed data
    prior_str = dm.get_priors(key)
    generate_prior_potentials(vars, prior_str, est_mesh)

    vars['observed_rates'] = []
    for d in data_list:
        age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
        age_weights = d.get('age_weights', pl.ones(len(age_indices)) / len(age_indices))

        lb, ub = dm.bounds_per_1(d)
        se = (pl.log(ub) - pl.log(lb)) / (2. * 1.96)
        if pl.isnan(se) or se <= 0.:
            se = 1.
        print 'data %d: log(value) = %f, se = %f' % (d['id'], pl.log(dm.value_per_1(d)), se)
        
        @mc.observed
        @mc.stochastic(name='obs_%d' % d['id'])
        def obs(f=vars['rate_stoch'],
                age_indices=age_indices,
                age_weights=age_weights,
                value=pl.log(dm.value_per_1(d)),
                tau=se**-2, data=d):
            f_i = rate_for_range(f, age_indices, age_weights)
            return mc.normal_like(value, pl.log(f_i), tau)
        vars['observed_rates'].append(obs)
        
    return vars
