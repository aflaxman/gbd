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
def print_mare(vars):
    if 'p_obs' in vars:
        are = pl.atleast_1d(pl.absolute((vars['p_obs'].value - vars['pi'].value)/vars['pi'].value))
        print 'mare:', pl.round_(pl.median(are), 2)

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
    map = mc.MAP(vars)
    m = mc.MCMC(vars)

    ## use MAP to generate good initial conditions
    try:
        method='fmin_powell'
        tol=.001
        verbose=1

        find_asr_initial_vals(vars, method, tol, verbose)
        map.fit(method=method, tol=tol, verbose=verbose)
        
        print_mare(vars)
        setup_asr_step_methods(m, vars)
    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    print_mare(vars)

    m.iter=iter
    m.burn=burn
    m.thin=thin
    try:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
    except TypeError:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval)

    m.wall_time = time.time() - start_time
    return map, m



param_types = 'i r f p pf rr smr m_with X'.split()
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
    map = mc.MAP(vars)
    m = mc.MCMC(vars)


    ## use MAP to generate good initial conditions
    try:
        method='fmin_powell'
        tol=.001
        verbose=1

        find_consistent_spline_initial_vals(vars, method, tol, verbose)

        for t in param_types:
            find_re_initial_vals(vars[t], method, tol, verbose)

        find_consistent_spline_initial_vals(vars, method, tol, verbose)

        for t in param_types:
            find_fe_initial_vals(vars[t], method, tol, verbose)

        find_consistent_spline_initial_vals(vars, method, tol, verbose)

        for t in param_types:
            find_dispersion_initial_vals(vars[t], method, tol, verbose)

        print 'fitting all stochs'
        map.fit(method=method, tol=tol, verbose=verbose)

        from fit_posterior import inspect_vars
        print inspect_vars({}, vars)

    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model

    max_knots = max([len(vars[t]['gamma']) for t in 'irf'])
    for i in range(max_knots):
        stoch = [vars[t]['gamma'][i] for t in 'ifr' if i < len(vars[t]['gamma'])]

        print 'finding Normal Approx for', [n.__name__ for n in stoch]
        vars_to_fit = [vars[t].get('p_obs'), vars[t].get('pi_sim'), vars[t].get('smooth_gamma'), vars[t].get('parent_similarity'),
                       vars[t].get('mu_sim'), vars[t].get('mu_age_derivative_potential'), vars[t].get('covariate_constraint')]
        na = mc.NormApprox(vars_to_fit + stoch)
        na.fit(method='fmin_powell', verbose=1)
        cov = pl.array(pl.inv(-na.hess), order='F')
        if pl.all(pl.eigvals(cov) >= 0):
            m.use_step_method(mc.AdaptiveMetropolis, stoch, cov=cov)
        else:
            print 'cov matrix is not positive semi-definite'
            m.use_step_method(mc.AdaptiveMetropolis, stoch)

    for t in param_types:
        setup_asr_step_methods(m, vars[t])

    m.iter=iter
    m.burn=burn
    m.thin=thin
    try:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
    except TypeError:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval)
    m.wall_time = time.time() - start_time

    return map, m



def find_consistent_spline_initial_vals(vars, method, tol, verbose):
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



def find_asr_initial_vals(vars, method, tol, verbose):
    for outer_reps in range(3):
        find_spline_initial_vals(vars, method, tol, verbose)
        find_re_initial_vals(vars, method, tol, verbose)
        find_spline_initial_vals(vars, method, tol, verbose)
        find_fe_initial_vals(vars, method, tol, verbose)
        find_spline_initial_vals(vars, method, tol, verbose)
        find_dispersion_initial_vals(vars, method, tol, verbose)

def find_spline_initial_vals(vars, method, tol, verbose):
    ## generate initial value by fitting knots sequentially
    vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                   vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]

    for i, n in enumerate(vars['gamma']):
        print 'fitting first %d knots of %d' % (i+1, len(vars['gamma']))
        vars_to_fit.append(n)
        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
        print_mare(vars)

def find_re_initial_vals(vars, method, tol, verbose):
    if 'hierarchy' not in vars:
        return

    for reps in range(3):
        for p in nx.traversal.bfs_tree(vars['hierarchy'], 'all'):
            successors = vars['hierarchy'].successors(p)
            if successors:
                print successors

                vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                               vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
                vars_to_fit += [vars.get('alpha_potentials')]

                re_vars = [vars['alpha'][vars['U'].columns.indexMap[n]] for n in successors + [p] if n in vars['U']]
                vars_to_fit += re_vars
                if len(re_vars) > 0:
                    mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

                print pl.round_([re.value for re in re_vars if isinstance(re, mc.Node)], 2)
                print_mare(vars)

    print 'sigma_alpha'
    vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                   vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
    vars_to_fit += [vars.get('sigma_alpha')]
    mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
    print pl.round_([s.value for s in vars['sigma_alpha']])
    print_mare(vars)


def find_fe_initial_vals(vars, method, tol, verbose):
    vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                   vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
    vars_to_fit += [vars.get('beta')]  # include fixed effects in sequential fit
    mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
    print_mare(vars)

def find_dispersion_initial_vals(vars, method, tol, verbose):
    vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                   vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
    vars_to_fit += [vars.get('eta'), vars.get('zeta')]
    mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
    print_mare(vars)


def setup_asr_step_methods(m, vars):
    # groups RE stochastics that are suspected of being dependent
    groups = []
    fe_group = [n for n in vars.get('beta', []) if isinstance(n, mc.Stochastic)]
    ap_group = [n for n in vars.get('gamma', []) if isinstance(n, mc.Stochastic)]
    groups += [[g_i, g_j] for g_i, g_j in zip(ap_group[1:], ap_group[:-1])] + [fe_group, ap_group, fe_group+ap_group]
    for a in vars.get('hierarchy', []):
        group = []
        if a in vars['U']:
            for b in nx.shortest_path(vars['hierarchy'], 'all', a):
                if b in vars['U']:
                    n = vars['alpha'][vars['U'].columns.indexMap[b]]
                    if isinstance(n, mc.Stochastic):
                        group.append(n)
        groups.append(group)
        #if len(group) > 0:
            #group += ap_group
            #groups.append(group)
            #group += fe_group
            #groups.append(group)
                    
    for stoch in groups:
        if len(stoch) > 0 and pl.all([isinstance(n, mc.Stochastic) for n in stoch]):
            # only step certain stochastics, for understanding convergence
            #if 'gamma_i' not in stoch[0].__name__:
            #    print 'no stepper for', stoch
            #    m.use_step_method(mc.NoStepper, stoch)
            #    continue

            print 'finding Normal Approx for', [n.__name__ for n in stoch]
            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
            na = mc.NormApprox(vars_to_fit + stoch)
            na.fit(method='fmin_powell', verbose=1)
            cov = pl.array(pl.inv(-na.hess), order='F')
            print 'opt:', pl.round_([n.value for n in stoch], 2)
            print 'cov:\n', cov.round(4)
            if pl.all(pl.eigvals(cov) >= 0):
                m.use_step_method(mc.AdaptiveMetropolis, stoch, cov=cov)
            else:
                print 'cov matrix is not positive semi-definite'
                m.use_step_method(mc.AdaptiveMetropolis, stoch)

