import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import clean, gbd_keys

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
    """
    if not hasattr(dm, 'vars'):
        print 'initializing model vars... ',
        dm.calc_effective_sample_size(dm.data)
        dm.vars = setup(dm, keys)
        print 'finished'

    if method == 'map':
        print 'initializing MAP object... ',
        map_method = 'fmin_powell'
        #map_method = 'fmin_l_bfgs_b'

        mc.MAP([dm.vars[k] for k in keys if k.find('incidence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if k.find('remission') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if k.find('excess-mortality') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if
                k.find('incidence') != -1 or
                k.find('bins') != -1 or
                k.find('prevalence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)
        mc.MAP([dm.vars[k] for k in keys if
                k.find('excess-mortality') != -1 or
                k.find('bins') != -1 or
                k.find('prevalence') != -1]).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)

        dm.map = mc.MAP(dm.vars)
        print 'finished'

        try:
            dm.map.fit(method=map_method, iterlim=500, tol=.001, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        
        for k in keys:
            try:
                val = dm.vars[k]['rate_stoch'].value
                dm.set_map(k, val)
            except KeyError:
                pass

    if method == 'norm_approx':
        dm.na = mc.NormApprox(dm.vars, eps=.0001)

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
        # make pymc warnings go to stdout
        import sys
        mc.warnings.warn = sys.stdout.write
        
        dm.mcmc = mc.MCMC(dm.vars)
        try:
            dm.mcmc.sample(iter=iter*thin+burn, thin=thin, burn=burn, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for k in keys:
            try:
                if dm.vars[k].has_key('rate_stoch'):
                    rate_model.store_mcmc_fit(dm, k, dm.vars[k]['rate_stoch'])
            except KeyError:
                pass


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

    # for each region-year-sex triple among the keys
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for('%s', r, y, s)
                if not key%'prevalence' in keys:
                    continue

                dm.set_units(key%'prevalence', '(per person)')
                dm.set_units(key%'duration', '(years)')
                for t in 'incidence', 'remission', 'excess-mortality':
                    dm.set_units(key%t, '(per person-year)')
                    #dm.get_initial_estimate(key%t, [d for d in dm.data if relevant_to(d, t, r, y, s)])

                data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]
                sub_vars = submodel.setup(dm, key, data)
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
    # ignore data if requested
    if d.get('ignore'):
        return False
    
    # check if data is of the correct type
    if t != 'all':
        if clean(d['data_type']).find(clean(t)) != 0:
            return False

    # check if data is from correct region
    if r != 'all' and r != 'world':
        if clean(d['gbd_region']) != clean(r) and clean(d['gbd_region']) != 'all':
            return False

    # check if data is from relevant year
    if y != 'all':
        y = int(y)
        if not y in [1990, 1997, 2005]:
            raise KeyError, 'GBD Year must be 1990 or 2005 (or 1997 for all years)'
        if y == 2005 and d['year_end'] < 1997:
            return False
        if y == 1990 and d['year_start'] > 1997:
            return False

    # check if data is for relevant sex
    if s != 'all':
        if clean(d['sex']) != clean(s) and clean(d['sex']) != 'all':
            return False

    # if code makes it this far, the data is relevent
    return True
