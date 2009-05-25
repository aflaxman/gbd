import numpy as np
import pymc as mc

import dismod3
import generic_disease_model as submodel
import beta_binomial_model as rate_model

output_data_types = ['Incidence', 'Remission', 'Case fatality', 'Prevalence', 'Duration']

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

    type : str, optional, one of 'incidence', 'remission', 'case_fatality'
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
    
    if method == 'map':
        dm.map = mc.MAP([dm.vars[dismod3.gbd_key_for(t, r, y, s)]
                         for t in types for r in regions for y in years for s in sexes])
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
        for t in types:
            for r in regions:
                for y in years:
                    for s in sexes:
                        key = dismod3.gbd_key_for(t, r, y, s)
                        dm.set_map(key, dm.vars[key]['rate_stoch'])
                        # TODO: set initial values from map estimate
                        # to save time in the future
                        
    elif method == 'mcmc':
        # TODO: make MAP object for selected submodel
        sub_var_list = [dm.vars[dismod3.gbd_key_for(t, r, y, s)]
                         for t in types for r in regions for y in years for s in sexes]
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
                        rate_model.store_mcmc_fit(dm, dm.vars[key]['rate_stoch'], key)


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
    * Sets dm's param_age_mesh and estimate_age_mesh, if they are not
      already set.

    * Sets the units of all estimates in the dm

    * Create PyMC variables for the generic disease model, and store
      them in dm.vars
    """

    # TODO: move these age_mesh defaults into DiseaseModel
    if dm.get_param_age_mesh() == []:
        dm.set_param_age_mesh([0.0, 10.0, 20.0, 30.0, 40.0,
                               50.0, 60.0, 70.0, 80.0, 90.0, 100.0])
    if dm.get_estimate_age_mesh() == []:
        dm.set_estimate_age_mesh(range(MAX_AGE))

    # find initial values for the rates that can be set
    for t in ['incidence', 'remission', 'case_fatality']:
        for r in dismod3.gbd_regions:
            for y in dismod3.gbd_years:
                for s in dismod3.gbd_sexes:
                    key = dismod3.gbd_key_for(t, r, y, s)
                    rate_model.initialize(dm, key) # TODO: include only appropriate data here
                    
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                t = 'prevalence'
                dm.set_units(dismod3.gbd_key_for(t, r, y, s), '(per person)')
                t = 'duration'
                dm.set_units(dismod3.gbd_key_for(t, r, y, s), '(years)')

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
                vars.update(submodel.setup(dm, r, y, s))

    # TODO: setup potentials on region/year/sex similarities
    
    return vars
