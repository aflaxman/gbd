""" Routines for fitting disease models"""
import sys
import time

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

## set number of threads to avoid overburdening cluster computers
try:
    import mkl
    mkl.set_num_threads(1)
except ImportError:
    pass

def fit_data_model(vars, iter, burn, thin, tune_interval):
    """ Fit data model using MCMC
    Input
    -----
    vars : dict

    Results
    -------
    returns a pymc.MCMC object created from vars, that has been fit with MCMC
    """
    start_time = time.time()
    ## use MAP to generate good initial conditions
    method='fmin_powell'
    tol=.001
    verbose=1
    try:
        map = mc.MAP(vars)

        ## generate initial value by fitting knots sequentially
        vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                       vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
        vars_to_fit += [vars.get('beta')]  # include fixed effects in sequential fit

        for i, n in enumerate(vars['gamma']):
            print 'fitting first %d knots of %d' % (i+1, len(vars['gamma']))
            vars_to_fit.append(n)
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        #vars_to_fit += [vars.get('eta'), vars.get('zeta')]
        #vars_to_fit += [vars.get('sigma_alpha'), vars.get('alpha'), vars.get('alpha_potentials')]
        #map = mc.MAP(vars_to_fit)
        map.fit(method=method, tol=tol, verbose=verbose)

    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)

    m.am_grouping = 'default'
    if m.am_grouping == 'alt1':
        m.use_step_method(mc.AdaptiveMetropolis, vars['beta'])
        #m.use_step_method(mc.AdaptiveMetropolis, [n for n in vars['alpha'] if isinstance(n, mc.Stochastic)])
        #m.use_step_method(mc.AdaptiveMetropolis, vars['sigma_alpha'])
        #m.use_step_method(mc.AdaptiveMetropolis, vars['gamma'][1:])

    elif m.am_grouping == 'alt2':
        #m.use_step_method(mc.AdaptiveMetropolis, vars['tau_alpha'])
        m.use_step_method(mc.AdaptiveMetropolis, vars['gamma'])
        m.use_step_method(mc.AdaptiveMetropolis, [vars[s] for s in 'alpha beta gamma_bar'.split() if isinstance(vars[s], mc.Stochastic)])

    elif m.am_grouping == 'alt3':
        #m.use_step_method(mc.AdaptiveMetropolis, vars['tau_alpha'])
        m.use_step_method(mc.AdaptiveMetropolis, [vars[s] for s in 'alpha beta gamma_bar gamma'.split() if isinstance(vars[s], mc.Stochastic)])

    m.iter=iter
    m.burn=burn
    m.thin=thin
    try:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
    except TypeError:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval)

    m.wall_time = time.time() - start_time
    return map, m



def fit_consistent_model(vars, iter, burn, thin, tune_interval):
    """ Fit data model using MCMC
    Input
    -----
    vars : dict

    Results
    -------
    returns a pymc.MCMC object created from vars, that has been fit with MCMC
    """
    start_time = time.time()
    param_types = 'i r f p pf rr smr m_with X'.split()

    ## use MAP to generate good initial conditions
    method='fmin_powell'
    tol=.001
    verbose=1
    try:
        map = mc.MAP(vars)

        ## generate initial value by fitting knots sequentially
        vars_to_fit = [vars['logit_C0']]
        for t in param_types:
            vars_to_fit += [vars[t].get('covariate_constraint'),
                            vars[t].get('mu_age_derivative_potential'), vars[t].get('mu_sim'),
                            vars[t].get('p_obs'), vars[t].get('parent_similarity'), vars[t].get('smooth_gamma'),]
        max_knots = max([len(vars[t]['gamma']) for t in 'irf'])
        for i in range(1, max_knots+1):
            print 'fitting first %d knots of %d' % (i, max_knots)
            vars_to_fit += [vars[t]['gamma'][:i] for t in 'irf']
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

            from fit_posterior import inspect_vars
            print inspect_vars({}, vars)[-10:]

        ## then fix effect coefficients for each rate separately
        for t in param_types:
            vars_to_fit = [vars[t].get('alpha'), vars[t].get('alpha_potentials'),
                           vars[t].get('beta'), vars[t].get('eta'), vars[t].get('zeta')]
            vars_to_fit += [vars[t].get('covariate_constraint'),
                            vars[t].get('mu_age_derivative_potential'), vars[t].get('mu_sim'),
                            vars[t].get('p_obs'), vars[t].get('parent_similarity'), vars[t].get('smooth_gamma'),]
            if pl.any([isinstance(n, mc.Stochastic) for n in vars_to_fit]):
                print 'fitting additional parameters for %s data' % t
                mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

                print inspect_vars({}, vars)[-10:]

        print 'fitting all stochs'
        map.fit(method=method, tol=tol, verbose=verbose)

        print inspect_vars({}, vars)

    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)

    m.am_grouping = 'alt2'

    if m.am_grouping == 'alt1':
        m.use_step_method(mc.AdaptiveMetropolis, [n for t in 'ifr' for n in vars[t]['gamma'][1:]])
        for t in param_types:
            #for node in 'tau_alpha':
            #    if isinstance(vars[t].get(node), mc.Stochastic):
            #        m.use_step_method(mc.AdaptiveMetropolis, var[t][node])

            # group all offset terms together in AdaptiveMetropolis
            print 'grouping stochastics'
            var_list = []
            for node in 'alpha beta gamma_bar gamma':
                if isinstance(vars[t].get(node), mc.Stochastic):
                    var_list.append(var[t][node])
            if len(var_list) > 0:
                m.use_step_method(mc.AdaptiveMetropolis, var_list)

    elif m.am_grouping == 'alt2':
        for i in range(max_knots):
            m.use_step_method(mc.AdaptiveMetropolis, [vars[t]['gamma'][i] for t in 'ifr' if i < len(vars[t]['gamma'])])


    m.iter=iter
    m.burn=burn
    m.thin=thin
    try:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
    except TypeError:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval)
    m.wall_time = time.time() - start_time

    return map, m


