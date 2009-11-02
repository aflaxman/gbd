import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import clean, gbd_keys
from dismod3.logit_gp_step import *

import generic_disease_model as submodel
import neg_binom_model as rate_model

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
        dm.vars = setup(dm, keys)
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
        map_method = 'fmin_powell'
        #map_method = 'fmin_l_bfgs_b'
        mc.MAP([dm.vars[k] for k in keys if k.find('incidence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if k.find('remission') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if k.find('case-fatality') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if
                k.find('incidence') != -1 or
                k.find('bins') != -1 or
                k.find('prevalence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if
                k.find('case-fatality') != -1 or
                k.find('bins') != -1 or
                k.find('prevalence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        
        dm.map = mc.MAP(sub_var_list)
        print 'finished'
        try:
            dm.map.fit(method=map_method, iterlim=500, tol=.001, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        
        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                val = dm.vars[k]['rate_stoch'].value
                dm.set_map(k, val)
                #dm.set_initial_value(k, val)  # better initial value may save time in the future

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
            dm.mcmc.sample(iter=5000, burn=0, thin=1, verbose=1)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])


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

                dm.set_units(key%'prevalence', '(per person)')
                dm.set_units(key%'duration', '(years)')
                for t in 'incidence', 'remission', 'case-fatality':
                    dm.set_units(key%t, '(per person-year)')
                    dm.fit_initial_estimate(key%t, [d for d in dm.data if relevant_to(d, t, r, y, s)])

                if not key%'prevalence' in keys:
                    continue
                data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]
                sub_vars = submodel.setup(dm, key, data)
                vars.update(sub_vars)

    world_key = dismod3.gbd_key_for('%s', 'world', '1997', 'total')
    sub_vars = submodel.setup(dm, world_key, [])
    vars.update(sub_vars)
    
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
