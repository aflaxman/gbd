""" Routines for fitting disease models"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx


def fit_data_model(vars):
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
    tol=1
    verbose=1
    try:
        vars_to_fit = [vars['gamma_bar'], vars['p_obs'], vars.get('pi_sim')]
        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        vars_to_fit.append(vars['gamma'])
        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        mc.MAP(vars).fit(method=method, tol=tol, verbose=verbose)
    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)
    for node in 'tau_alpha alpha beta gamma':
        if isinstance(vars.get(node), mc.Stochastic):
            m.use_step_method(mc.AdaptiveMetropolis, var[node])

    m.sample(15000, 10000, 5)

    return m



def fit_consistent_model(vars, iter=5000, burn=0, thin=1):
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
    tol=1
    verbose=1
    try:
        vars_to_fit = [vars['logit_C0']] \
            + [vars[t].get('gamma_bar') for t in param_types] \
            + [vars[t].get('p_obs') for t in param_types]
        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

        mc.MAP(vars).fit(method=method, tol=tol, verbose=verbose)
    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)
    for t in param_types:
        for node in 'tau_alpha alpha beta gamma':
            if isinstance(vars[t].get(node), mc.Stochastic):
                m.use_step_method(mc.AdaptiveMetropolis, var[t][node])

    m.sample(iter, burn, thin, verbose=verbose-1)

    return m


