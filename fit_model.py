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
    are = pl.absolute((vars['p_obs'].value - vars['pi'].value)/vars['pi'].value)
    print 'mare:', pl.median(are).round(2)

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

        for outer_reps in range(3):
            ## generate initial value by fitting knots sequentially
            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]

            for i, n in enumerate(vars['gamma']):
                print 'fitting first %d knots of %d' % (i+1, len(vars['gamma']))
                vars_to_fit.append(n)
                mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
                print_mare(vars)

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
                        mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)

                        print pl.round_([re.value for re in re_vars], 2)
                        print_mare(vars)

            print 'sigma_alpha'
            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
            vars_to_fit += [vars.get('sigma_alpha')]
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
            print pl.round_([s.value for s in vars['sigma_alpha']])
            print_mare(vars)

            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]

            for i, n in enumerate(vars['gamma']):
                print 're-fitting first %d knots of %d' % (i+1, len(vars['gamma']))
                vars_to_fit.append(n)
                mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
                print_mare(vars)


            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
            vars_to_fit += [vars.get('beta')]  # include fixed effects in sequential fit
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
            print_mare(vars)

            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]

            for i, n in enumerate(vars['gamma']):
                print 'fitting first %d knots of %d' % (i+1, len(vars['gamma']))
                vars_to_fit.append(n)
                mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
                print_mare(vars)

            vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                           vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
            vars_to_fit += [vars.get('eta'), vars.get('zeta')]
            mc.MAP(vars_to_fit).fit(method=method, tol=tol, verbose=verbose)
            print_mare(vars)

        map.fit(method=method, tol=tol, verbose=verbose)
        print_mare(vars)

    except KeyboardInterrupt:
        print 'Initial condition calculation interrupted'

    ## use MCMC to fit the model
    m = mc.MCMC(vars)

    for stoch in [vars['alpha'], vars['gamma']]:
        print 'finding Normal Approx for', stoch
        vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                       vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
        na = mc.NormApprox(vars_to_fit + stoch)
        na.fit(method='fmin_powell', verbose=1)
        m.use_step_method(mc.AdaptiveMetropolis, stoch, cov=pl.array(.95*pl.inv(-na.hess)+.05*pl.eye(len(na.hess)), order='F'))

    m.iter=iter*50
    m.burn=burn*50
    m.thin=thin*50
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

    for i in range(1, max_knots):
        stoch = [vars[t]['gamma'][i] for t in 'ifr' if i < len(vars[t]['gamma'])] + [vars[t]['gamma'][i-1] for t in 'ifr' if i < len(vars[t]['gamma'])]
        print 'finding Normal Approx for', stoch
        vars_to_fit = [vars.get('p_obs'), vars.get('pi_sim'), vars.get('smooth_gamma'), vars.get('parent_similarity'),
                       vars.get('mu_sim'), vars.get('mu_age_derivative_potential'), vars.get('covariate_constraint')]
        na = mc.NormApprox(vars_to_fit + stoch)
        na.fit(method='fmin_powell', verbose=1)
        m.use_step_method(mc.AdaptiveMetropolis, stoch, cov=pl.array(.95*pl.inv(-na.hess)+.05*pl.eye(len(na.hess)), order='F'))

    m.iter=iter*5
    m.burn=burn*5
    m.thin=thin*5
    try:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval, progress_bar=True, progress_bar_fd=sys.stdout)
    except TypeError:
        m.sample(m.iter, m.burn, m.thin, tune_interval=tune_interval)
    m.wall_time = time.time() - start_time

    return map, m


