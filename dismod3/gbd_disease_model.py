import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import clean, gbd_keys
from dismod3.logit_gp_step import *

import generic_disease_model as submodel
#import beta_binomial_model as rate_model
import logit_normal_model as rate_model

def fit(dm, method='map', keys=gbd_keys(), iter=50000, burn=25000, thin=1, verbose=1):
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

    keys : list, optional
      a list of gbd keys for the parameters to fit; it can speed up
      computation to holding some parameters constant while allowing
      others to vary

    iter : int, optional
    burn : int, optional
    thin : int, optional
      parameters for the MCMC, which control how long it takes, and
      how accurate it is
      
    Example
    -------
    >>> import dismod3
    >>> import dismod3.gbd_disease_model as model
    >>> dm = dismod3.get_disease_model(1)
    >>> dm.params['estimate_type'] = 'fit region-year-sex individually'
    >>> keys = model.gbd_keys(region_list=['australasia'], year_list=[1990], sex_list=['male'])
    >>> keys += model.gbd_keys(region_list=['north_america_high_income'], year_list=[1990], sex_list=['male'])
    >>> keys += model.gbd_keys(region_list=['world'], year_list=['total'], sex_list=['total'])
    >>> model.fit(dm, method='map', keys=keys)
    >>> model.fit(dm, method='norm_approx', keys=keys)
    >>> model.fit(dm, method='mcmc', keys=keys)
    """
    if not hasattr(dm, 'vars'):
        print 'initializing model vars... ',
        initialize(dm, keys)
        print 'finished'

    sub_var_list = [dm.vars[k] for k in keys]

    # remove similarity potentials, if fitting submodels individually
    if dm.params.get('estimate_type').find('individually') != -1:
        for vl in sub_var_list:
            for k in vl.keys():
                if k.find('similarity') != -1:
                    vl.pop(k)

    if method == 'map':
        print 'making MAP object... ',
        dm.map = mc.MAP(sub_var_list)
        print 'finished'
        try:
            dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=verbose)
            #dm.map.fit(method='fmin_l_bfgs_b', iterlim=500, tol=.00001, verbose=verbose)  # ~ twice as fast as powell's method
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        
        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                val = dm.vars[k]['rate_stoch'].value
                dm.set_map(k, val)
                dm.set_initial_value(k, val)  # better initial value may save time in the future

    if method == 'norm_approx':
        dm.na = mc.NormApprox(sub_var_list, eps=.0001)

        try:
            dm.na.fit(method='fmin_powell', iterlim=500, tol=.00001, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                dm.set_map(k, dm.vars[k]['rate_stoch'].value)

        try:
            dm.na.sample(1000, verbose=verbose)
            for k in keys:
                # TODO: rename 'rate_stoch' to something more appropriate
                if dm.vars[k].has_key('rate_stoch'):
                    rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

                        
    elif method == 'mcmc':
        dm.mcmc = mc.MCMC(sub_var_list)
        try:
            dm.mcmc.sample(iter=100, burn=0, thin=1, verbose=1)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])


def initialize(dm, keys):
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

                    if not key in keys:
                        continue

                    dm.set_units(key, '(per person-year)')

                    if dm.has_initial_value(key):
                        continue

                    data[key] = [d for d in dm.data if relevant_to(d, t, r, y, s)]

                    # use a subset of potentially relevant data if there is a lot of it,
                    # to speed things up
                    initialization_data = random_shuffle(data[key]) \
                                          + random_shuffle([d for d in dm.data if relevant_to(d, t, 'all', 'all', 'all') and not d in data[key]])
        
                    if len(initialization_data) > 25:
                        dm.fit_initial_estimate(key, initialization_data[:25])
                    else:
                        dm.fit_initial_estimate(key, initialization_data)
                    
                    
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                dm.set_units(dismod3.gbd_key_for('prevalence', r, y, s), '(per person)')
                dm.set_units(dismod3.gbd_key_for('duration', r, y, s), '(years)')

    dm.vars = setup(dm, keys)

def random_shuffle(x):
    import copy, random
    y = copy.copy(x)
    random.shuffle(y)
    return y

def similarity_prior(name, v1, v2):
    """ Generate a PyMC potential for the similarity of to age-specific rate functions
    """
    @mc.potential(name=name)
    def similarity(r1=v1['rate_stoch'], r2=v2['rate_stoch'],
                        d1=v1['dispersion'], d2=v2['dispersion']):
        return mc.normal_like(np.diff(mc.logit(r1)) - np.diff(mc.logit(r2)), 0., 100. / (d1**2 + d2**2))
    return similarity

def setup(dm, keys):
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

    # for each region-year-sex triple, create stochastic vars for a
    # generic disease submodel
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for('%s', r, y, s)
                if not key%'prevalence' in keys:
                    continue
                data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]
                sub_vars = submodel.setup(dm, key, data)
                vars.update(sub_vars)

    # link regional estimates together through a hierarchical model,
    # which models the difference between delta(region rate) and delta(world rate)
    # as a mean-zero gaussian, with precision = conf(region rate) + conf(world rate)
    world_key = dismod3.gbd_key_for('%s', 'world', '1997', 'total')
    sub_vars = submodel.setup(dm, world_key, [])
    vars.update(sub_vars)
#      for t in ['incidence', 'remission', 'case-fatality']:
#         vars[world_key % t]['h_potentials'] = []

#     for t in ['incidence', 'remission', 'case-fatality']:
#         for r in dismod3.gbd_regions:
#             for s in dismod3.gbd_sexes:
#                     k1 = dismod3.gbd_key_for(t, r, '1990', s)
#                     k2 = dismod3.gbd_key_for(t, r, '2005', s)
#                     vars[k1]['time_similarity'] = similarity_prior('time_similarity_%s_%s' % (k1, k2), vars[k1], vars[k2])

#             for y in dismod3.gbd_years:
#                     k1 = dismod3.gbd_key_for(t, r, y, 'male')
#                     k2 = dismod3.gbd_key_for(t, r, y, 'female')
#                     vars[k1]['sex_similarity'] = similarity_prior('sex_similarity_%s_%s' % (k1, k2), vars[k1], vars[k2])

    
    return vars

def relevant_to(d, t, r, y, s):
    """ Determine if data is relevant to specified type, region, year, and sex
    
    Parameters
    ----------
    d : data hash
    t : str, one of 'incidence data', 'prevalence data', etc... or 'all'
    r : str, one of 21 GBD regions or 'all'
    y : int, one of 1990, 2005 or 'all'
    s : sex, one of 'male', 'female' or 'all'
    """
    # check if data is of the correct type
    if t != 'all':
        if clean(d['data_type']).find(clean(t)) == -1:
            return False

    # check if data is from correct region
    if r != 'all':
        if clean(d['gbd_region']) != clean(r) and clean(d['gbd_region']) != 'all':
            return False

    # check if data is from relevant year
    if y != 'all':
        y = int(y)
        if not y in [1990, 2005]:
            raise KeyError, 'GBD Year must be 1990 or 2005'
        if y == 2005 and d['year_end'] < 1997:
            return False
        if y == 1990 and d['year_start'] > 1997:
            return False

    # check if data is for relevant sex
    if y != 'all':
        if clean(d['sex']) != clean('total') and clean(d['sex']) != clean(s):
            return False

    # if code makes it this far, the data is relevent
    return True
