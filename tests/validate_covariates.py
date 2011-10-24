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
reload(data)

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
    ## set simulation parameters
    data_type = 'p'
    N = 100
    delta_true = .5


    ## generate simulation data
    import dismod3
    import simplejson as json
    model = data.ModelData.from_gbd_jsons(json.loads(dismod3.disease_json.DiseaseJson().to_json()))
    model.parameters['p']['parameter_age_mesh'] = [0, 100]

    model.input_data = pandas.DataFrame(index=range(N))
    model.input_data['age_start'] = 0
    model.input_data['age_end'] = 0
    model.input_data['year_start'] = 2005.
    model.input_data['year_end'] = 2005.
    model.input_data['sex'] = 'total'
    model.input_data['data_type'] = data_type
    model.input_data['standard_error'] = pl.nan
    model.input_data['upper_ci'] = pl.nan
    model.input_data['lower_ci'] = pl.nan


    # 1. choose pi^true
    pi_true = mc.runiform(.0001, .1)
    a = pl.arange(0, 100, 1)
    pi_age_true = pi_true * pl.ones_like(a)


    # 2. choose sigma^true
    sigma_true = mc.runiform(.05, .5, 4)


    # 3. choose alpha^true
    alpha = dict(all=0.)
    sum_sr = 0.
    for sr in model.hierarchy['all']:
        sum_r = 0.
        for r in model.hierarchy[sr]:
            sum_c = 0.
            for c in model.hierarchy[r]:
                alpha[c] = mc.rnormal(0., sigma_true[3]**-2.)
                sum_c += alpha[c]
            alpha[c] -= sum_c

            alpha[r] = mc.rnormal(0., sigma_true[2]**-2.)
            sum_r += alpha[r]
        alpha[r] -= sum_r

        alpha[sr] = mc.rnormal(0., sigma_true[1]**-2.)
        sum_sr += alpha[sr]
    alpha[sr] -= sum_sr


    # 4. choose observed prevalence values
    model.input_data['effective_sample_size'] = mc.runiform(100, 10000, N)

    area_list = pl.array(nx.traversal.bfs_tree(model.hierarchy, 'australasia').nodes() + nx.traversal.bfs_tree(model.hierarchy, 'north_america_high_income').nodes())
    model.input_data['area'] = area_list[mc.rcategorical(pl.ones(len(area_list)) / float(len(area_list)), N)]

    model.input_data['value_true'] = pl.nan
    for i, a in model.input_data['area'].iteritems():
        model.input_data['value_true'][i] = pi_true * pl.exp(pl.sum([alpha[n] for n in nx.shortest_path(model.hierarchy, 'all', a)]))

    n = model.input_data['effective_sample_size']
    p = model.input_data['value_true']
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true*n*p) / n



    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('re_validation', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=5, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)

    pl.show()

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']

    model.re = pandas.DataFrame(index=[n for n in nx.traversal.dfs_preorder_nodes(model.hierarchy)])
    model.re['true'] = pandas.Series(dict(alpha))
    model.re['mu_pred'] = pandas.Series([n.stats()['mean'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.re['sigma_pred'] = pandas.Series([n.stats()['standard deviation'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.re['abs_err'] = model.re['true'] - model.re['mu_pred']
    model.re['rel_err'] = (model.re['true'] - model.re['mu_pred']) / model.re['true']
    model.re['coverage'] = (model.re['true'] >= model.re['mu_pred'] - 1.96*model.re['sigma_pred']) & (model.re['true'] <= model.re['mu_pred'] + 1.96*model.re['sigma_pred'])

    print 'Random Effects'
    print model.re.dropna()

    model.sigma = pandas.DataFrame(dict(true=sigma_true,
                                        mu_pred=[n.stats()['mean'] for n in model.vars['p']['sigma_alpha']]))
    print model.sigma

    return model


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
