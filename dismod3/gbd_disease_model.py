import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import clean

import generic_disease_model as submodel
import beta_binomial_model as rate_model

output_data_types = ['Incidence', 'Remission', 'Case-fatality', 'Prevalence', 'Duration']

def fit(dm, method='map', type='all', region='all', year='all', sex='all'):
    """ Generate an estimate of the generic disease model parameters
    using maximum a posteriori liklihood (MAP) or Markov-chain Monte
    Carlo (MCMC)

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    method : string, optional
      the parameter estimation method, either 'map' or 'mcmc'

    type : str, optional, one of 'incidence', 'remission', 'case-fatality'
    region : str, optional, one of 21 GBD regions or 'all'
    year : str, optional, one of '1995', '2005', 'all'
    sex : str, optional, one of 'male', 'female', 'all'
      speed up computation time by only fitting the submodel given by type, region, year, sex
      
    Example
    -------
    >>> import dismod3
    >>> import dismod3.gbd_disease_model as model
    >>> dm = dismod3.get_disease_model(1)
    >>> model.fit(dm, method='map')
    >>> model.fit(dm, method='mcmc')
    """
    
    if not hasattr(dm, 'vars'):
        initialize(dm)

    if type == 'all':
        types = output_data_types
    elif type == 'prevalence':
        # prevalence is controlled by incidence, remission, and case-fatality
        types = output_data_types
    else:
        types = [ type ]

    if region == 'all':
        regions = dismod3.gbd_regions
    else:
        regions = [ region ]

    if year == 'all':
        years = dismod3.gbd_years
    else:
        years = [ year ]

    if sex == 'all':
        sexes = dismod3.gbd_sexes
    else:
        sexes = [ sex ]

    sub_var_list = []
    for t in types:
        for r in regions:
            for y in years:
                for s in sexes:
                    key = dismod3.gbd_key_for(t, r, y, s)
                    sub_var_list.append(dm.vars[key])

    if method == 'map':
        dm.map = mc.MAP(sub_var_list)
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        for t in types:
            for r in regions:
                for y in years:
                    for s in sexes:
                        key = dismod3.gbd_key_for(t, r, y, s)
                        dm.set_map(key, dm.vars[key]['rate_stoch'].value)
                        # set initial values from map estimate to save
                        # time in the future
                        dm.set_initial_value(key, dm.vars[key]['rate_stoch'].value)
                        
    elif method == 'mcmc':
        # TODO: make MAP object for selected submodel
        dm.mcmc = mc.MCMC(sub_var_list)
        for v in sub_var_list:
            if len(v.get('logit_p_stochs', [])) > 0:
                dm.mcmc.use_step_method(
                    mc.AdaptiveMetropolis, v['logit_p_stochs'])
                    
        dm.mcmc.sample(iter=60*1000, burn=10*1000, thin=50, verbose=1)
        for t in types:
            for r in regions:
                for y in years:
                    for s in sexes:
                        key = dismod3.gbd_key_for(t, r, y, s)
                        rate_model.store_mcmc_fit(dm, key, dm.vars[key]['rate_stoch'])
                        # set initial values from map estimate to save
                        # time in the future
                        dm.set_initial_value(key, dm.vars[key]['rate_stoch'].stats()['mean'])


def initialize(dm):
    """ Initialize the stochastic and deterministic random variables
    for the multi-region/year/sex generic disease model

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    * Sets the units of all estimates in the dm

    * Create PyMC variables for the generic disease model, and store
      them in dm.vars
    """

    # find initial values for the rates that can be set
    data = {}
    for t in ['incidence', 'remission', 'case-fatality']:
        for r in dismod3.gbd_regions:
            for y in dismod3.gbd_years:
                for s in dismod3.gbd_sexes:
                    key = dismod3.gbd_key_for(t, r, y, s)
                    data[key] = \
                        [d for d in dm.data if
                         clean(d['data_type']).find(clean(t)) != -1 and
                         clean(d['gbd_region']) == clean(r) and
                         ((y == 2005 and d['year_end'] >= 2000) or
                          (y == 1995 and d['year_start'] < 2000)) and
                         d['sex'] == s]

                    # use a random subset of the data if there is a lot of it,
                    # to speed things up
                    if len(data[key]) > 25:
                        dm.fit_initial_estimate(key, random.sample(data[key], 25))
                    else:
                        dm.fit_initial_estimate(key, data[key])
                    dm.set_units(key, '(per person-year)')
                    
                    
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                dm.set_units(dismod3.gbd_key_for('prevalence', r, y, s), '(per person)')
                dm.set_units(dismod3.gbd_key_for('duration', r, y, s), '(years)')

    dm.vars = setup(dm)


def setup(dm):
    """ Generate the PyMC variables for a multi-region/year/sex generic
    disease model.

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      multi-region/year/sex generic disease model.
    """
    
    vars = {}
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for('%s', r, y, s)
                data = [d for d in dm.data if
                        clean(d['gbd_region']) == clean(r) and
                        ((y == 2005 and d['year_end'] >= 2000) or
                         (y == 1995 and d['year_start'] < 2000)) and
                        d['sex'] == s]
                sub_vars = submodel.setup(dm, key, data)
                vars.update(sub_vars)

    # TODO: setup potentials on region/year/sex similarities
    
    return vars
