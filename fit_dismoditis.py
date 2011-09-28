""" Fit a simulation data set"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import consistent_model
reload(consistent_model)
reload(data_model)

def plot_model_data(vars):
    """ 2x2 tile plot of data intervals for consistent model"""
    for j, t in enumerate('i r f p pf'.split()):
        pl.subplot(2, 3, j+1)
        pl.title(t)
        age_start = vars[t]['data']['age_start']
        age_end = vars[t]['data']['age_end']
        p = vars[t]['data']['value']
        for a_0i, a_1i, p_i in zip(age_start, age_end, p):
            pl.plot([a_0i, a_1i], [p_i,p_i], 'ks-', mew=1, mec='w', ms=4)


def plot_model_params(vars, i, ages=pl.arange(101)):
    """ 2x2 tile plot of params for consistent model"""
    for j, t in enumerate('irfp'):
        pl.subplot(2, 3, j+1)
        pl.plot(ages, vars[t]['mu_age'].value, color=pl.cm.spectral((i+.01)/10))


def demo_model_fit(vars, n_maps=0, n_mcmcs=2, ages=pl.arange(101)):
    """ fit model, showing the process
    don't try to run to completion, this is just for testing
    """

    map = mc.MAP(vars)
    for i  in range(n_maps):
        map.fit(method='fmin_powell', verbose=1, iterlim=1)
        plot_model_params(vars, i, ages)

    m = mc.MCMC(vars)
    m.use_step_method(mc.AdaptiveMetropolis, [vars[k]['gamma_bar'] for k in 'irf'] + [vars[k]['gamma'] for k in 'irf'])
    for i in range(n_mcmcs):
        m.sample(1001)
        plot_model_params(vars, i, ages)
       
def fit_model(vars, ages=pl.arange(101)):
    try:
        map = mc.MAP([vars[t]['gamma_bar'] for t in 'irf'] +
                     [vars[t]['p_obs'] for t in 'i r f p pf'.split() if 'p_obs' in vars[t]])
        map.fit(method='fmin_powell', verbose=1, tol=.01)


        map = mc.MAP(vars)
        map.fit(method='fmin_powell', verbose=1, tol=.01)
    except KeyboardInterrupt:
        pass

    m = mc.MCMC(vars)
    for k in 'i r f p pf'.split():
        for node in 'tau_alpha alpha gamma':
            if isinstance(vars[k].get(node), mc.Stochastic):
                m.use_step_method(mc.AdaptiveMetropolis, var[k][node])

    m.sample(11100, 10000, 10)

    for j, t in enumerate('i r f p pf'.split()):
        pl.subplot(2, 3, j+1)
        pl.plot(ages, vars[t]['mu_age'].stats()['mean'], 'k-', linewidth=2)
        pl.plot(ages, vars[t]['mu_age'].stats()['95% HPD interval'], 'k--')
       
    return m

if __name__ == '__main__':

    d = data.ModelData.from_gbd_json('tests/dismoditis.json')

    d.parameters['p']['level_value']['age_before'] = 1
    d.parameters['i']['level_value']['age_before'] = 25
    d.parameters['f']['level_value']['age_before'] = 25
    d.parameters['r']['level_value']['age_before'] = 100
    for t in 'irfp':
        d.parameters[t]['smoothness']['amount'] = 'Moderately'

    knots = pl.arange(30, 86, 5)
    ages = pl.arange(knots[0], knots[-1] + 1)
    for t in 'irfp':
        d.parameters[t]['parameter_age_mesh'] = knots

    #d.input_data['standard_error'] /= 10.

    # create model and priors for top level of hierarchy
    root = 'all'
    vars = consistent_model.consistent_model(d.input_data, d.parameters, d.hierarchy, root, ages=ages)

    pl.figure()
    plot_model_data(vars)
    #demo_model_fit(vars, 3, 3)
    m1 = fit_model(vars, ages)
    
    # generate estimates for super-region_5, male, 2005
    est_trace = {}
    priors = {}
    for t in 'irfp':
        est_trace[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'super-region_5', 'male', 2005, vars[t])
        priors[t] = est_trace[t].mean(0)

    for j, t in enumerate('irfp'):
        pl.subplot(2, 3, j+1)
        pl.plot(ages, priors[t], color='r', linewidth=1)

    # create model and priors for (asia_southeast, male, 2005), including estimate of
    # super-region_5 to borrow strength
    root = 'asia_southeast'
    subtree = nx.traversal.bfs_tree(d.hierarchy, root)
    relevant_rows = [i for i, r in d.input_data.T.iteritems() if r['area'] in subtree and r['year_end'] >= 1997 and r['sex'] in ['male', 'total']]
    vars = consistent_model.consistent_model(d.input_data.ix[relevant_rows], d.parameters, d.hierarchy, root, priors=priors, ages=ages)

    # fit initial conditions to data
    mc.MAP([vars['logit_C0'], vars['p']]).fit(tol=.01, verbose=1)

    pl.figure()
    plot_model_data(vars)
    for j, t in enumerate('irfp'):
        pl.subplot(2, 3, j+1)
        pl.plot(ages, priors[t], color='r', linewidth=1)

    #demo_model_fit(vars, 3, 3)
    m2 = fit_model(vars, ages)

    # generate estimates for THA, male, 2005
    posteriors = {}
    for t in 'irfp':
        posteriors[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'THA', 'male', 2005, vars[t])

    for j, t in enumerate('irfp'):
        pl.subplot(2, 3, j+1)
        pl.plot(ages, posteriors[t].mean(0), color='b', linewidth=1)

