""" Routines for fitting disease models"""
import sys

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

## set number of threads to avoid overburdening cluster computers
import mkl
mkl.set_num_threads(2)

def fit_data_model(vars, iter=15000, burn=5000, thin=90, tune_interval=1000):
    """ Fit data model using MCMC
    Input
    -----
    vars : dict

    Results
    -------
    returns a pymc.MCMC object created from vars, that has been fit with MCMC
    """

    ## use MAP to generate good initial conditions
    method='fmin_powell'
    tol=.001
    verbose=1
    try:
        vars_to_fit = [vars['gamma_bar'], vars['p_obs'], vars.get('pi_sim'), vars.get('smooth_gamma'),
                       vars.get('mu_sim'), vars.get('mu_age_derivative_potential')]
        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        for i, n in enumerate(vars['gamma'][1:]):  # skip first knot on list, since it is not a stoch
            print 'fitting first %d knots of %d' % (i+2, len(vars['gamma']))
            vars_to_fit.append(n)
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
        
        mc.MAP(vars).fit(method=method, tol=tol, verbose=verbose)
    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)

    m.am_grouping = 'default'
    if m.am_grouping == 'alt1':
        for s in 'alpha beta gamma'.split():
            m.use_step_method(mc.AdaptiveMetropolis, vars[s])

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
    m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)

    return m



def fit_consistent_model(vars, iter=50350, burn=15000, thin=350, tune_interval=1000):
    """ Fit data model using MCMC
    Input
    -----
    vars : dict

    Results
    -------
    returns a pymc.MCMC object created from vars, that has been fit with MCMC
    """
    param_types = 'i r f p pf rr'.split()

    ## use MAP to generate good initial conditions
    method='fmin_powell'
    tol=.001
    verbose=1
    try:
        vars_to_fit = [vars['logit_C0']]
        for t in param_types:
            vars_to_fit += [vars[t].get('gamma_bar'), vars[t].get('p_obs'), vars[t].get('pi_sim'), vars[t].get('smooth_gamma'),
                            vars[t].get('mu_sim'), vars[t].get('mu_age_derivative_potential')]

        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        max_knots = max([len(vars[t].get('gamma', [])) for t in param_types])
        for i in range(1, max_knots):  # skip first knot on list, since it is not a stoch
            print 'fitting first %d knots of %d' % (i+1, max_knots)
            vars_to_fit += [vars[t].get('gamma')[:i] for t in 'irf']
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        print 'fitting all stochs'
        mc.MAP(vars).fit(method=method, tol=tol, verbose=verbose)
    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)
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

    m.sample(iter, burn, thin, verbose=verbose-1, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)

    return m


