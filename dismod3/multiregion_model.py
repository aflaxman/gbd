import numpy as np
import pymc as mc

from model_utils import *
from bayesian_models import probabilistic_utils
import beta_binomial_model as rate_model

def fit(dm, method='map', data_type='prevalence data'):
    """ Generate an estimate of multiregion beta binomial model parameters
    using maximum a posteriori liklihood (MAP) or Markov-chain Monte
    Carlo (MCMC)

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    method : string, optional
      the parameter estimation method, either 'map' or 'mcmc'

    data_type : str, optional
      Only data in dm.data with d['data_type'] == data_type will be
      included in each beta-binomial liklihood function

    Example
    -------
    >>> import dismod3
    >>> import dismod3.multiregion_model as model
    >>> dm = dismod3.get_disease_model(849)
    >>> model.fit(dm, method='map')
    >>> model.fit(dm, method='mcmc')
    """
    if not hasattr(dm, 'vars'):
        initialize(dm, data_type)

    if method == 'map':
        if not hasattr(dm, 'map'):
            dm.map = mc.MAP(dm.vars)
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        for r in dm.data_by_region.keys() + ['World']:
            dm.set_map(rate_key(data_type,r),
                       dm.vars[rate_key(data_type,r)]['rate_stoch'].value)
    elif method == 'mcmc':
        if not hasattr(dm, 'mcmc'):
            dm.mcmc = mc.MCMC(dm.vars)
            for r in dm.data_by_region.keys() + ['World']:
                logit_p_stochs = dm.vars[rate_key(data_type,r)]['logit_p_stochs']
                if len(logit_p_stochs) > 0:
                    dm.mcmc.use_step_method(
                        mc.AdaptiveMetropolis, logit_p_stochs)
                    
        dm.mcmc.sample(iter=40000, burn=10000, thin=30, verbose=1)
        for r in dm.data_by_region:
            rate_model.store_mcmc_fit(dm, dm.vars[rate_key(data_type,r)]['rate_stoch'],
                                      rate_key(data_type,r))


def initialize(dm, data_type='prevalence data'):
    """ Initialize the stochastic and deterministic random variables
    for the multiregion beta binomial model of age-specific rate functions

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    data_type : str, optional
      Only data in dm.data with d['data_type'] == data_type will be
      included in the beta-binomial liklihood functions
    
    Results
    -------
    * Sets dm's param_age_mesh and estimate_age_mesh, if they are not
      already set.

    * Sets the units of dm

    * Create PyMC variables for the multiregion beta binomial model, and stores
      them in dm.vars
    """
    if dm.get_param_age_mesh() == []:
        dm.set_param_age_mesh([0.0, 10.0, 20.0, 30.0, 40.0,
                               50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
    if dm.get_estimate_age_mesh() == []:
        dm.set_estimate_age_mesh(range(MAX_AGE))

    # sort the data by GBD regions
    dm.data_by_region = {}
    for d in dm.data:
        if d['data_type'] != data_type:
            continue
        r = d['gbd_region']
        dm.data_by_region[r] = dm.data_by_region.get(r, []) + [d]

    # find initial values for the data from each region
    for r, r_data in dm.data_by_region.items():
        # use a random subset of the data if there is a lot of it,
        # to speed things up
        data = [d for d in r_data]
        if len(data) > 25:
            import random
            data = random.sample(data,25)

        dm.set_units(rate_key(data_type, r), '(per person-year)')
        dm.fit_initial_estimate(rate_key(data_type, r), data)

    # set initial world value to average of regional values
    avg_value = np.mean(
        [dm.get_initial_value(rate_key(data_type, r)) \
             for r in dm.data_by_region.keys()],
        axis=0)
    dm.set_units(rate_key(data_type, 'World'), '(per person-year)')
    dm.set_initial_value(rate_key(data_type, 'World'), avg_value)

    dm.vars = setup(dm, data_type)


def setup(dm, data_type='prevalence data'):
    """ Generate the PyMC variables for a multiregion beta binomial
    model of a rate function that varies by GBD region

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    data_type : str, optional
      Name stochs according to this string
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      beta binomial model.  dm.vars['rate_stochs'] is itself a
      dictionary of stochs, keyed by rate_key(data_type,region).

    Details
    -------
    The beta binomial model models for each regional rate are linked
    together by assuming that they are all realizations of a single
    underlying rate function
    """
    
    vars = {}
    rate_stochs = {}

    # set up individual rate stochs
    for r in dm.data_by_region.keys() + ['World']:
        stoch_key = rate_key(data_type, r)
        vars[stoch_key] = rate_model.setup(dm, dm.data_by_region.get(r, []), stoch_key)

    world_key = rate_key(data_type, 'World')
    world_rate = vars[world_key]['rate_stoch']

    # link regional estimates together through a hierarchical model,
    # where each region rate is a realization of the world
    # rate with gaussian noise
    for r in dm.data_by_region.keys():
        stoch_key = rate_key(data_type, r)
        @mc.potential(name='hierarchical_potential_%s'%stoch_key)
        def hier_potential(x=vars[stoch_key]['rate_stoch'], y=world_rate):
            return mc.normal_like(x-y, 0., 100.)
        vars[stoch_key]['h_potential'] = hier_potential

    return vars


def rate_key(data_type, region):
    """ Make a human-readable dictionary key"""
    return '%s+%s' % (data_type, region)


def plot(dm, data_type='prevalence data'):
    """ Plot the results of the multiregion beta binomial
    model fit by GBD region

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    """
    est = {}
    for r in dm.data_by_region.keys() + ['World']:
        est[r] = dm.vars[rate_key(data_type, r)]['rate_stoch'].value

    
    for r in sorted(est.keys(), key=lambda est: est[k][-1], reverse=True):
        plot(dm.get_estimate_age_mesh(),
             est[r],
             linewidth=(r == 'World') and 6 or 3,
             alpha=(r == 'World') and .75 or .5,
             label=r)

    legend(loc=(1.1,-.1), pad=0, handletextsep=0)
