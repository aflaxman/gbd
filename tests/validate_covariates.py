""" Validate covariate model

These tests are use randomized computation, so they might fail
occasionally due to stochastic variation
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import networkx as nx
import pandas

import data
import rate_model
import age_pattern
import covariate_model
import data_simulation
import data_model
import fit_model
import graphics

reload(covariate_model)

def validate_covariate_model_fe():
    # generate simulated data
    data_type = 'p'
    n = 1000
    sigma_true = .01
    a = pl.arange(0, 100, 1)
    pi_age_true = .05 * pl.ones_like(a)

    d = data.ModelData()
    d.input_data = data_simulation.simulated_age_intervals(data_type, n, a, pi_age_true, sigma_true)

    # add fixed effect to simulated data
    X = mc.rnormal(0., 1.**-2, size=(n,3))
    beta_true = [-.1, .1, .2]
    Y_true = pl.dot(X, beta_true)

    d.input_data['value'] = pl.maximum(0., d.input_data['value'] * pl.exp(Y_true))
    d.input_data['x_0'] = X[:,0]
    d.input_data['x_1'] = X[:,1]
    d.input_data['x_2'] = X[:,2]

    # adjust hierarchy and parameters
    d.hierarchy, d.output_template = data_simulation.small_output()
    d.parameters['p']['parameter_age_mesh'] = [0, 100]


    # create model and priors
    vars = data_model.data_model('fe_sim', d, 'p', 'all', 'total', 'all', None, None)


    # fit model
    m = fit_model.fit_data_model(vars, iter=2000, burn=1000, thin=10, tune_interval=100)

    # compare estimate to ground truth (skip endpoints, because they are extra hard to get right)
    print 'mean(beta) close to truth:', pl.allclose(m.beta.stats()['mean'], beta_true, atol=.05)
    lb, ub = m.beta.stats()['95% HPD interval'].T
    print 'probability coverage(beta) more than 50%:', pl.mean((lb <= beta_true) & (beta_true <= ub)) > .5

    print
    print 'truth:', pl.round_(beta_true, 2)
    print 'est:', m.beta.stats()['mean'].round(2)
    print 'ui:', m.beta.stats()['95% HPD interval'].round(2)

    graphics.plot_one_type(d, vars, {}, 'p')
    graphics.plot_one_ppc(vars, 'p')
    graphics.plot_convergence_diag(vars)

    return m, vars, d

def validate_covariate_model_re():
    # generate simulated data
    data_type = 'p'
    n = 100
    sigma_true = .01
    a = pl.arange(0, 100, 1)
    pi_age_true = .05 * pl.ones_like(a)

    d = data.ModelData()
    d.input_data = data_simulation.simulated_age_intervals(data_type, n, a, pi_age_true, sigma_true)
    d.input_data = d.input_data.drop('x_0 x_1 x_2'.split(), 1)

    # setup hierarchy
    hierarchy, d.output_template = data_simulation.small_output()

    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')
    hierarchy.add_edge('all', 'USA', weight=.1)
    hierarchy.add_edge('all', 'CAN', weight=.1)
    d.hierarchy = hierarchy

    # shift data differentially by area
    area_list = pl.array(['all', 'USA', 'CAN'])
    d.input_data['area'] = area_list[mc.rcategorical([.3, .3, .4], n)]

    sex_list = pl.array(['male', 'female', 'total'])
    sex = sex_list[mc.rcategorical([.3, .3, .4], n)]

    year = pl.array(mc.runiform(1990, 2010, n), dtype=int)
        
    alpha_true = dict(all=0., USA=.1, CAN=-.1)

    d.input_data['value'] = d.input_data['value'] * pl.exp([alpha_true[a] for a in d.input_data['area']])



    # adjust parameters
    d.parameters['p']['parameter_age_mesh'] = [0, 100]

    # create model and priors
    vars = data_model.data_model('re_sim', d, 'p', 'all', 'total', 'all', None, None)


    # fit model
    m = fit_model.fit_data_model(vars, iter=2000, burn=1000, thin=10, tune_interval=100)

    graphics.plot_one_type(d, vars, {}, 'p')
    pl.plot(covariate_model.predict_for(d.output_template, d.hierarchy, 'all', 'total', 'all', 'USA', 'male', 1990, vars).T, color='red', zorder=-100, linewidth=2, alpha=.1)
    pl.plot(covariate_model.predict_for(d.output_template, d.hierarchy, 'all', 'total', 'all', 'USA', 'male', 1990, vars).mean(0), color='red')
    pl.plot(covariate_model.predict_for(d.output_template, d.hierarchy, 'all', 'total', 'all', 'CAN', 'male', 1990, vars).T, color='blue', zorder=-100, linewidth=2, alpha=.1)
    pl.plot(covariate_model.predict_for(d.output_template, d.hierarchy, 'all', 'total', 'all', 'CAN', 'male', 1990, vars).mean(0), color='blue')


    #graphics.plot_one_ppc(vars, 'p')
    graphics.plot_convergence_diag(vars)

    # print results
    print 'alpha_true:', pl.round_(alpha_true.values(), 2)
    print 'alpha:', pl.round_([n.stats()['mean'] for n in m.alpha], 2)
    print 'sigma_alpha:', pl.round_([n.stats()['mean'] for n in m.sigma_alpha], 2)


    return m, vars, d

def test_covariate_model_dispersion():
    # simulate normal data
    n = 100

    hierarchy = nx.DiGraph()
    hierarchy.add_node('all')

    Z = mc.rcategorical([.5, 5.], n)
    zeta_true = -.2

    pi_true = .1
    ess = 10000
    eta_true = pl.log(50)
    delta_true = 50 + pl.exp(eta_true)

    p = mc.rnegative_binomial(pi_true*ess, delta_true*pl.exp(Z*zeta_true)) / float(ess)

    data = pandas.DataFrame(dict(value=p, z_0=Z))
    data['area'] = 'all'
    data['sex'] = 'total'
    data['year_start'] = 2000
    data['year_end'] = 2000

    # create model and priors
    vars = dict(mu=mc.Uninformative('mu_test', value=pi_true))
    vars.update(covariate_model.mean_covariate_model('test', vars['mu'], data, hierarchy, 'all'))
    vars.update(covariate_model.dispersion_covariate_model('test', data))
    vars.update(rate_model.neg_binom_model('test', vars['pi'], vars['delta'], p, ess))

    # fit model
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    m = mc.MCMC(vars)
    m.sample(20000, 10000, 10)

    # print summary results
    print 'mu:', m.mu.stats()
    print 'eta:', m.eta.stats()
    print 'zeta:', m.zeta.stats()

    # compare estimate to ground truth (skip endpoints, because they are extra hard to get right)
    #assert pl.allclose(m.zeta.stats()['mean'], zeta_true, rtol=.2)
    lb, ub = m.zeta.stats()['95% HPD interval'].T
    assert lb <= zeta_true <= ub

    lb, ub = m.eta.stats()['95% HPD interval'].T
    assert lb <= eta_true <= ub

if __name__ == '__main__':
    m, vars, model = validate_covariate_model_fe()
    m, vars, model = validate_covariate_model_re()
