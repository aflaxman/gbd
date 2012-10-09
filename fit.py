""" Module for DisMod-MR model fitting methods"""

import time
import pymc as mc

import dismod3

import fit_model
reload(fit_model)

def fit_asr(model, data_type, iter=2000, burn=1000, thin=1, tune_interval=100, verbose=False):
    """ Fit data model for one epidemiologic parameter using MCMC
    
    :Parameters:
      - `model` : data.ModelData
      - `data_type` : str, one of 'i', 'r', 'f', 'p', or 'pf'
      - `iter` : int, number of posterior samples fit
      - `burn` : int, number of posterior samples to discard as burn-in
      - `thin` : int, samples thinned by this number
      - `tune_interval` : int
      - `verbose` : boolean

    :Results:
      - returns a pymc.MCMC object created from vars, that has been fit with MCMC

    .. note::
      - `burn` must be less than `iter`
      - `thin` must be less than `iter` minus `burn`

    """
    assert burn < iter, 'burn must be less than iter'
    assert thin < iter - burn, 'thin must be less than iter-burn'

    vars = model.vars[data_type]
    
    start_time = time.time()
    map = mc.MAP(vars)
    m = mc.MCMC(vars)

    ## use MAP to generate good initial conditions
    try:
        method='fmin_powell'
        tol=.001

        fit_model.logger.info('finding initial values')
        fit_model.find_asr_initial_vals(vars, method, tol, verbose)

        fit_model.logger.info('\nfinding MAP estimate')
        map.fit(method=method, tol=tol, verbose=verbose)
        
        if verbose:
            fit_model.print_mare(vars)
        fit_model.logger.info('\nfinding step covariances estimate')
        fit_model.setup_asr_step_methods(m, vars)

        fit_model.logger.info('\nresetting initial values (1)')
        fit_model.find_asr_initial_vals(vars, method, tol, verbose)
        fit_model.logger.info('\nresetting initial values (2)\n')
        map.fit(method=method, tol=tol, verbose=verbose)
    except KeyboardInterrupt:
        fit_model.logger.warning('Initial condition calculation interrupted')

    ## use MCMC to fit the model
    fit_model.print_mare(vars)

    fit_model.logger.info('sampling from posterior\n')
    m.iter=iter
    m.burn=burn
    m.thin=thin
    if verbose:
        try:
            m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
        except TypeError:
            m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=False, verbose=verbose)
    else:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=False)

    m.wall_time = time.time() - start_time
    
    model.map = map
    model.mcmc = m
    
    return model.map, model.mcmc

# TODO: move fit_model.fit_consistent_model to fit.fit_consistent
import fit_model
def fit_consistent(model, iter=2000, burn=1000, thin=1, tune_interval=100, verbose=False):
    """Fit data model for all epidemiologic parameters using MCMC
    
    :Parameters:
      - `model` : data.ModelData
      - `iter` : int, number of posterior samples fit
      - `burn` : int, number of posterior samples to discard as burn-in
      - `thin` : int, samples thinned by this number
      - `tune_interval` : int
      - `verbose` : boolean

    :Results:
      - returns a pymc.MCMC object created from vars, that has been fit with MCMC

    .. note::
      - `burn` must be less than `iter`
      - `thin` must be less than `iter` minus `burn`

    """
    model.map, model.mcmc = fit_model.fit_consistent_model(model.vars, iter, burn, thin, tune_interval, verbose)

