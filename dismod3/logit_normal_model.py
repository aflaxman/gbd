import numpy as np
import pymc as mc

from dismod3.utils import debug, interpolate, rate_for_range, indices_for_range, generate_prior_potentials
from dismod3.settings import NEARLY_ZERO, MISSING

# re-use the beta_binomial_model's store_mcmc_fit function
# (might need to override this in the future)
from dismod3.beta_binomial_model import store_mcmc_fit

def setup(dm, key, data_list, rate_stoch=None):
    """ Generate the PyMC variables for a logit-normal model of
    a single rate function

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    key : str
      the name of the key for everything about this model (priors,
      initial values, estimations)

    data_list : list of data dicts
      the observed data to use in the logit-normal liklihood function

    rate_stoch : pymc.Stochastic, optional
      a PyMC stochastic (or deterministic) object, with
      len(rate_stoch.value) == len(dm.get_estimation_age_mesh()).
      This is used to link rate stochs into a larger model,
      for example.

    Results
    -------
    vars : dict
      Return a dictionary of all the relevant PyMC objects for the
      rate model.  vars['rate_stoch'] is of particular
      relevance; this is what is used to link the rate model
      into more complicated models, like the generic disease model.

    Details
    -------
    The logit-normal model parameters are the following:
      * the mean age-specific rate function
      * systematic error in this mean
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    if np.any(np.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    # create varible for interpolated logit rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        @mc.deterministic(name='interp_logit(%s)' % key)
        def interp_logit_rate(rate_stoch=rate_stoch):
            return mc.logit(rate_stoch)
    else:
        param_mesh = dm.get_param_age_mesh()
        initial_value = dm.get_initial_value(key)

        # find the logit of the initial values, which is a little bit
        # of work because initial values are sampled from the est_mesh,
        # but the logit_initial_values are needed on the param_mesh
        logit_initial_value = mc.logit(
            interpolate(est_mesh, initial_value, param_mesh))
        
        logit_rate = mc.Normal('logit(%s)' % key,
                               mu=-5.*np.ones(len(param_mesh)),
                               tau=1.e-2,
                               value=logit_initial_value)
        #logit_rate = [mc.Normal('logit(%s)_%d' % (key, a), mu=-5., tau=1.e-2) for a in param_mesh]
        vars['logit_rate'] = logit_rate

        @mc.deterministic(name='interp_logit(%s)' % key)
        def interp_logit_rate(logit_rate=logit_rate):
            return interpolate(param_mesh, logit_rate, est_mesh)

        @mc.deterministic(name=key)
        def rate_stoch(interp_logit_rate=interp_logit_rate):
            return mc.invlogit(interp_logit_rate)
    vars['interp_logit_rate'] = interp_logit_rate
    vars['rate_stoch'] = rate_stoch

    # create stochastic variable for over-dispersion "random effect"
    overdispersion = mc.Gamma('overdispersion_%s' % key, alpha=10., beta=10. / .001)
    vars['overdispersion'] = overdispersion

    # create potentials for priors
    vars['priors'] = generate_prior_potentials(dm.get_priors(key), est_mesh, rate_stoch)

    # create observed stochastics for data
    vars['data'] = data_list
    vars['observed_rates'] = []

    min_val = min([1.e-9] + [dm.value_per_1(d) for d in data_list if dm.value_per_1(d) > 0]) # TODO: assess validity of this minimum value
    max_se = max([.000001] + [dm.se_per_1(d) for d in data_list if dm.se_per_1(d) > 0])  # TODO: assess validity of this maximum std err

    for d in data_list:
        try:
            age_indices, age_weights, logit_val, logit_se = values_from(dm, d, min_val, max_se)
        except ValueError:
            continue

        @mc.observed
        @mc.stochastic(name='data_%d' % d['id'])
        def obs(value=logit_val,
                logit_rate=interp_logit_rate,
                logit_se=logit_se,
                overdispersion=overdispersion,
                age_indices=age_indices,
                age_weights=age_weights):
            mean_val = rate_for_range(logit_rate, age_indices, age_weights)
            return mc.normal_like(x=value, mu=mean_val, tau=1. / (overdispersion**2 + logit_se**2))
            
        vars['observed_rates'].append(obs)
        
    return vars

def values_from(dm, d, min_val=1.e-5, max_se=.1):
    """ Extract the normalized values from a piece of data

    Parameters
    ----------
    dm : disease model

    d : data dict

    min_val : float, optional
      the value to use instead of zero, since logit cannot model true zero

    max_se : float, optional
      the standard error to use for data with missing or zero standard error
    """
    est_mesh = dm.get_estimate_age_mesh()

    # get the index vector and weight vector for the age range
    age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
    age_weights = d['age_weights']

    # ensure all rate data is valid
    d_val = dm.value_per_1(d)
    if d_val < 0 or d_val > 1:
        debug('WARNING: data %d not in range (0,1)' % d[id])
        raise ValueError
    elif d_val == 0.:
        d_val = min_val / 10.  # TODO: determine if this is an acceptible way to deal with zero
    elif d_val == 1.:
        d_val = 1. - min_val / 10.

    logit_val = mc.logit(d_val)

    d_se = dm.se_per_1(d)
    if d_se == MISSING:
        d_se = max_se * 10. #TODO: determine if this is an acceptible way to deal with missing
    elif d_se == 0.:
        d_se = max_se

    logit_se = (1/d_val + 1/(1-d_val)) * d_se

    return age_indices, age_weights, logit_val, logit_se
