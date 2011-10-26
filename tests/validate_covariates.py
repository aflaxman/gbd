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
    m = fit_model.fit_data_model(vars, iter=20000, burn=10000, thin=10, tune_interval=100)

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


def validate_covariate_model_re(N=500, delta_true=.15, pi_true=.01, sigma_true = [.1,.1,.1,.1,.1]):
    ## set simulation parameters
    import dismod3
    import simplejson as json
    model = data.ModelData.from_gbd_jsons(json.loads(dismod3.disease_json.DiseaseJson().to_json()))
    data_type = 'p'
    model.parameters[data_type]['parameter_age_mesh'] = [0, 100]
    area_list = []
    for sr in sorted(model.hierarchy.successors('all')):
        area_list.append(sr)
        for r in sorted(model.hierarchy.successors(sr)):
            area_list.append(r)
            area_list += sorted(model.hierarchy.successors(r))[:5]
    area_list = pl.array(area_list)


    ## generate simulation data
    model.input_data = pandas.DataFrame(index=range(N))
    model.input_data['age_start'] = 0
    model.input_data['age_end'] = 1
    model.input_data['year_start'] = 2005.
    model.input_data['year_end'] = 2005.
    model.input_data['sex'] = 'total'
    model.input_data['data_type'] = data_type
    model.input_data['standard_error'] = pl.nan
    model.input_data['upper_ci'] = pl.nan
    model.input_data['lower_ci'] = pl.nan

    # choose alpha^true
    alpha = dict(all=0.)
    sum_sr = 0.
    last_sr = -1
    for sr in model.hierarchy['all']:
        if sr not in area_list:
            continue

        sum_r = 0.
        last_r = -1
        for r in model.hierarchy[sr]:
            if r not in area_list:
                continue

            sum_c = 0.
            last_c = -1
            for c in model.hierarchy[r]:
                if c not in area_list:
                    continue

                alpha[c] = mc.rnormal(0., sigma_true[3]**-2.)
                sum_c += alpha[c]
                last_c = c
            if last_c >= 0:
                alpha[last_c] -= sum_c

            alpha[r] = mc.rnormal(0., sigma_true[2]**-2.)
            sum_r += alpha[r]
            last_r = r
        if last_r >= 0:
            alpha[last_r] -= sum_r

        alpha[sr] = mc.rnormal(0., sigma_true[1]**-2.)
        sum_sr += alpha[sr]
        last_sr = sr
    if last_sr >= 0:
        alpha[last_sr] -= sum_sr

    # choose observed prevalence values
    model.input_data['effective_sample_size'] = mc.runiform(100, 10000, N)

    model.input_data['area'] = area_list[mc.rcategorical(pl.ones(len(area_list)) / float(len(area_list)), N)]

    model.input_data['true'] = pl.nan
    for i, a in model.input_data['area'].iteritems():
        model.input_data['true'][i] = pi_true * pl.exp(pl.sum([alpha[n] for n in nx.shortest_path(model.hierarchy, 'all', a) if n in alpha]))

    n = model.input_data['effective_sample_size']
    p = model.input_data['true']
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true*n*p) / n



    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=20000, burn=10000, thin=10, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)

    pl.show()

    def metrics(df):
        df['abs_err'] = df['true'] - df['mu_pred']
        df['rel_err'] = (df['true'] - df['mu_pred']) / df['true']
        df['covered?'] = (df['true'] >= df['mu_pred'] - 1.96*df['sigma_pred']) & (df['true'] <= df['mu_pred'] + 1.96*df['sigma_pred'])

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    metrics(model.input_data)


    model.alpha = pandas.DataFrame(index=[n for n in nx.traversal.dfs_preorder_nodes(model.hierarchy)])
    model.alpha['true'] = pandas.Series(dict(alpha))
    model.alpha['mu_pred'] = pandas.Series([n.stats()['mean'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.alpha['sigma_pred'] = pandas.Series([n.stats()['standard deviation'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    metrics(model.alpha)

    print '\nalpha'
    print model.alpha.dropna()


    model.sigma = pandas.DataFrame(dict(true=sigma_true))
    model.sigma['mu_pred'] = [n.stats()['mean'] for n in model.vars['p']['sigma_alpha']]
    model.sigma['sigma_pred']=[n.stats()['standard deviation'] for n in model.vars['p']['sigma_alpha']]
    metrics(model.sigma)

    print 'sigma_alpha'
    print model.sigma

    
    model.results = dict(param=[], bias=[], mare=[], mae=[], pc=[])
    def add_to_results(df, name):
        model.results['param'].append(name)
        model.results['bias'].append(df['abs_err'].mean())
        model.results['mae'].append((pl.median(pl.absolute(df['abs_err'].dropna()))))
        model.results['mare'].append(pl.median(pl.absolute(df['rel_err'].dropna())))
        model.results['pc'].append(df['covered?'].mean())
    add_to_results(model.sigma[1:3], 'sigma')

    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    metrics(model.delta)

    print 'delta'
    print model.delta
    add_to_results(model.delta, 'delta')

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())
    print 'effect prediction MAE: %.3f, coverage: %.2f' % (pl.median(pl.absolute(model.alpha['abs_err'].dropna())),
                                                          model.alpha.dropna()['covered?'].mean())
    add_to_results(model.input_data, 'input_data')
    add_to_results(model.alpha, 'alpha')

    model.results = pandas.DataFrame(model.results)
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
    model = validate_covariate_model_re()
