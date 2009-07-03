import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import clean, gbd_keys
from dismod3.logit_gp_step import *

import generic_disease_model as submodel
import beta_binomial_model as rate_model
#import logit_normal_model as rate_model

def fit(dm, method='map', keys=gbd_keys(), iter=1000, burn=10*1000, thin=50, verbose=0):
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
    >>> keys = model.gbd_keys(region_list=['australasia'], year_list=[1990], sex_list=['male'])
    >>> keys += model.gbd_keys(region_list=['north_america_high_income'], year_list=[1990], sex_list=['male'])
    >>> keys += model.gbd_keys(region_list=['world'], year_list=['total'], sex_list=['total'])
    >>> model.fit(dm, method='map', keys=keys)
    >>> model.fit(dm, method='mcmc', keys=keys)
    """
    if not hasattr(dm, 'vars'):
        initialize(dm)

    sub_var_list = [dm.vars[k] for k in keys]

    # remove similarity potentials, if fitting submodels individually
    if dm.params.get('estimate_type').find('individually') != -1:
        for vl in sub_var_list:
            for k in vl.keys():
                if k.find('similarity') != -1:
                    vl.pop(k)

    if method == 'map':
        dm.map = mc.MAP(sub_var_list)
        try:
            dm.map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
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

        dm.na.fit(method='fmin_powell', verbose=verbose)
        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                dm.set_map(k, dm.vars[k]['rate_stoch'].value)

        dm.na.sample(1000, verbose=verbose)
        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])

                        
    elif method == 'mcmc':
        dm.mcmc = mc.MCMC(sub_var_list)
        for v in sub_var_list:
            if len(v.get('latent_p', [])) > 0:
                dm.mcmc.use_step_method(mc.AdaptiveMetropolis, v['latent_p'], verbose=verbose)
            if v.get('logit_rate'):
                lr = v['logit_rate']

                ## logit_gp_step is a variant of hit-and-run that uses
                ## metropolis rejection instead of a Gibbs step
                dm.mcmc.use_step_method(LogitGPStep, lr, dm=dm, key=v['rate_stoch'].__name__, data_list=v['data'], verbose=verbose)

                ## adaptive metropolis also works fine, possibly it is slower to mix (possibly faster...)
                #dm.mcmc.use_step_method(mc.AdaptiveMetropolis, lr, verbose=verbose)

                # pick a smooth initial value
                sm = LogitGPStep(lr, dm=dm, key=v['rate_stoch'].__name__, data_list=v['data'])
                lr.value = sm.random()  # FIXME:  is this doing anything?

        try:
            dm.mcmc.sample(iter=thin*iter+burn, burn=burn, thin=thin, verbose=1)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        try:
            for k in keys:
                if dm.vars[k].has_key('rate_stoch'):
                    rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])
                    # better initial value may save time in the future
                    dm.set_initial_value(k, dm.vars[k]['rate_stoch'].stats()['mean'])
        except IndexError:
            # if user cancels with cntl-c before burn-in is completed,
            # save attempt will raise an IndexError, because trace is
            # empty
            pass


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

                    dm.set_units(key, '(per person-year)')

                    if dm.has_initial_value(key):
                        continue
                    
                    data[key] = [d for d in dm.data if relevant_to(d, t, r, y, s)]

                    # use a subset of potentially relevant data if there is a lot of it,
                    # to speed things up
                    initialization_data = random_shuffle(data[key]) \
                                          + random_shuffle([d for d in dm.data if relevant_to(d, t, 'all', 'all', 'all')])
                    if len(initialization_data) > 25:
                        dm.fit_initial_estimate(key, initialization_data[:25])
                    else:
                        dm.fit_initial_estimate(key, initialization_data)
                    
                    
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                dm.set_units(dismod3.gbd_key_for('prevalence', r, y, s), '(per person)')
                dm.set_units(dismod3.gbd_key_for('duration', r, y, s), '(years)')

    dm.vars = setup(dm)

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
                        d1=v1['overdispersion'], d2=v2['overdispersion']):
        return mc.normal_like(r1 - r2, 0., 1. / .01**2)
        return mc.normal_like(np.diff(np.log(r1)) - np.diff(np.log(r2)), 0., 1. / .1**2)
    return similarity

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

    # for each region-year-sex triple, create stochastic vars for a
    # generic disease submodel
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for('%s', r, y, s)
                data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]
                sub_vars = submodel.setup(dm, key, data)
                vars.update(sub_vars)

    # link regional estimates together through a hierarchical model,
    # which models the difference between delta(region rate) and delta(world rate)
    # as a mean-zero gaussian, with precision = conf(region rate) + conf(world rate)
    world_key = dismod3.gbd_key_for('%s', 'world', 'total', 'total')
    sub_vars = submodel.setup(dm, world_key, [])
    vars.update(sub_vars)
    for t in ['incidence', 'remission', 'case-fatality']:
        vars[world_key % t]['h_potentials'] = []

    for t in ['incidence', 'remission', 'case-fatality']:
        for r in dismod3.gbd_regions:
            for s in dismod3.gbd_sexes:
                    k1 = dismod3.gbd_key_for(t, r, '1990', s)
                    k2 = dismod3.gbd_key_for(t, r, '2005', s)
                    vars[k1]['time_similarity'] = similarity_prior('time_similarity_%s_%s' % (k1, k2), vars[k1], vars[k2])

            for y in dismod3.gbd_years:
                    k1 = dismod3.gbd_key_for(t, r, y, 'male')
                    k2 = dismod3.gbd_key_for(t, r, y, 'female')
                    vars[k1]['sex_similarity'] = similarity_prior('sex_similarity_%s_%s' % (k1, k2), vars[k1], vars[k2])

    
    return vars

def relevant_to(d, t, r, y, s):
    """ Determine if data is relevant to specified type, region, year, and sex
    
    Parameters
    ----------
    d : data hash
    t : str, one of 'incidence data', 'prevalence data', etc... or 'all'
    r : str, one of 21 GBD regions
    y : int, one of 1990, 2005
    s : sex, one of 'male', 'female'
    """
    if t != 'all' and clean(d['data_type']).find(clean(t)) == -1:
        return False
    if clean(d['gbd_region']) != clean(r):
        return False
    if y == 2005 and d['year_end'] < 1997:
        return False
    if y == 1990 and d['year_start'] > 1997:
        return False
    if clean(d['sex']) != clean('total') and clean(d['sex']) != clean(s):
        return False
    
    return True
