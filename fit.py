""" Module for DisMod-MR model fitting methods"""

import time
import pymc as mc
import pylab as pl

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
    assert burn < iter, 'burn must be less than iter'
    assert thin < iter - burn, 'thin must be less than iter-burn'

    param_types = 'i r f p pf rr smr m_with X'.split()

    vars = model.vars
    
    start_time = time.time()
    map = mc.MAP(vars)
    m = mc.MCMC(vars)

    ## use MAP to generate good initial conditions
    try:
        method='fmin_powell'
        tol=.001

        fit_model.logger.info('fitting submodels')
        fit_model.find_consistent_spline_initial_vals(vars, method, tol, verbose)

        for t in param_types:
            fit_model.find_re_initial_vals(vars[t], method, tol, verbose)
            fit_model.logger.info('.')

        fit_model.find_consistent_spline_initial_vals(vars, method, tol, verbose)
        fit_model.logger.info('.')

        for t in param_types:
            fit_model.find_fe_initial_vals(vars[t], method, tol, verbose)
            fit_model.logger.info('.')

        fit_model.find_consistent_spline_initial_vals(vars, method, tol, verbose)
        fit_model.logger.info('.')

        for t in param_types:
            fit_model.find_dispersion_initial_vals(vars[t], method, tol, verbose)
            fit_model.logger.info('.')

        fit_model.logger.info('\nfitting all stochs\n')
        map.fit(method=method, tol=tol, verbose=verbose)

        if verbose:
            from fit_posterior import inspect_vars
            print inspect_vars({}, vars)

    except KeyboardInterrupt:
        fit_model.logger.warning('Initial condition calculation interrupted')

    ## use MCMC to fit the model

    try:
        fit_model.logger.info('finding step covariances')
        vars_to_fit = [[vars[t].get('p_obs'), vars[t].get('pi_sim'), vars[t].get('smooth_gamma'), vars[t].get('parent_similarity'),
                        vars[t].get('mu_sim'), vars[t].get('mu_age_derivative_potential'), vars[t].get('covariate_constraint')] for t in param_types]
        max_knots = max([len(vars[t]['gamma']) for t in 'irf'])
        for i in range(max_knots):
            stoch = [vars[t]['gamma'][i] for t in 'ifr' if i < len(vars[t]['gamma'])]

            if verbose:
                print 'finding Normal Approx for', [n.__name__ for n in stoch]
            try:
                na = mc.NormApprox(vars_to_fit + stoch)
                na.fit(method='fmin_powell', verbose=verbose)
                cov = pl.array(pl.inv(-na.hess), order='F')
                if pl.all(pl.eigvals(cov) >= 0):
                    m.use_step_method(mc.AdaptiveMetropolis, stoch, cov=cov)
                else:
                    raise ValueError
            except ValueError:
                if verbose:
                    print 'cov matrix is not positive semi-definite'
                m.use_step_method(mc.AdaptiveMetropolis, stoch)

            fit_model.logger.info('.')

        for t in param_types:
            fit_model.setup_asr_step_methods(m, vars[t], vars_to_fit)

            # reset values to MAP
            fit_model.find_consistent_spline_initial_vals(vars, method, tol, verbose)
            fit_model.logger.info('.')
        map.fit(method=method, tol=tol, verbose=verbose)
        fit_model.logger.info('.')
    except KeyboardInterrupt:
        fit_model.logger.warning('Initial condition calculation interrupted')

    fit_model.logger.info('\nsampling from posterior distribution\n')
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