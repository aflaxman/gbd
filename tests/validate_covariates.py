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

def initialize_input_data(input_data):
    input_data['age_start'] = 0
    input_data['age_end'] = 1
    input_data['year_start'] = 2005.
    input_data['year_end'] = 2005.
    input_data['sex'] = 'total'
    input_data['data_type'] = 'p'
    input_data['standard_error'] = pl.nan
    input_data['upper_ci'] = pl.nan
    input_data['lower_ci'] = pl.nan
    input_data['area'] = 'all'


def add_quality_metrics(df):
    df['abs_err'] = df['true'] - df['mu_pred']
    df['rel_err'] = (df['true'] - df['mu_pred']) / df['true']
    df['covered?'] = (df['true'] >= df['mu_pred'] - 1.96*df['sigma_pred']) & (df['true'] <= df['mu_pred'] + 1.96*df['sigma_pred'])


def add_to_results(model, name):
    df = getattr(model, name)
    model.results['param'].append(name)
    model.results['bias'].append(df['abs_err'].mean())
    model.results['mae'].append((pl.median(pl.absolute(df['abs_err'].dropna()))))
    model.results['mare'].append(pl.median(pl.absolute(df['rel_err'].dropna())))
    model.results['pc'].append(df['covered?'].mean())


def validate_covariate_model_fe(N=100, delta_true=3, pi_true=.01, beta_true=[.5, -.5, 0.], replicate=0):
    # set random seed for reproducibility
    mc.np.random.seed(1234567 + replicate)
    
    ## generate simulated data
    a = pl.arange(0, 100, 1)
    pi_age_true = pi_true * pl.ones_like(a)

    model = data.ModelData()
    model.parameters['p']['parameter_age_mesh'] = [0, 100]
    model.input_data = pandas.DataFrame(index=range(N))
    initialize_input_data(model.input_data)

    # add fixed effect to simulated data
    X = mc.rnormal(0., 1.**-2, size=(N,len(beta_true)))
    Y_true = pl.dot(X, beta_true)

    for i in range(len(beta_true)):
        model.input_data['x_%d'%i] = X[:,i]
    model.input_data['true'] = pi_true * pl.exp(Y_true)

    model.input_data['effective_sample_size'] = mc.runiform(100, 10000, N)

    n = model.input_data['effective_sample_size']
    p = model.input_data['true']
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true) / n


    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=5, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)

    pl.show()

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    add_quality_metrics(model.input_data)


    model.beta = pandas.DataFrame(index=model.vars['p']['X'].columns)
    model.beta['true'] = 0.
    for i in range(len(beta_true)):
        model.beta['true']['x_%d'%i] = beta_true[i]
    
    model.beta['mu_pred'] = [n.stats()['mean'] for n in model.vars['p']['beta']]
    model.beta['sigma_pred'] = [n.stats()['standard deviation'] for n in model.vars['p']['beta']]
    add_quality_metrics(model.beta)

    print '\nbeta'
    print model.beta
    
    model.results = dict(param=[], bias=[], mare=[], mae=[], pc=[])
    add_to_results(model, 'beta')

    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    add_quality_metrics(model.delta)

    print 'delta'
    print model.delta
    add_to_results(model, 'delta')

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())
    print 'effect prediction MAE: %.3f, coverage: %.2f' % (pl.median(pl.absolute(model.beta['abs_err'].dropna())),
                                                           model.beta.dropna()['covered?'].mean())
    add_to_results(model, 'input_data')
    add_to_results(model, 'beta')

    model.results = pandas.DataFrame(model.results)
    return model



def alpha_true_sim(model, area_list, sigma_true):
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

    return alpha


def validate_covariate_model_re(N=500, delta_true=.15, pi_true=.01, sigma_true = [.1,.1,.1,.1,.1], ess=1000):
    ## set simulation parameters
    import dismod3
    import simplejson as json
    model = data.ModelData.from_gbd_jsons(json.loads(dismod3.disease_json.DiseaseJson().to_json()))
    model.parameters['p']['parameter_age_mesh'] = [0, 100]
    model.parameters['p']['heterogeneity'] = 'Slightly'  # ensure heterogeneity is slightly

    area_list = []
    for sr in sorted(model.hierarchy.successors('all')):
        area_list.append(sr)
        for r in sorted(model.hierarchy.successors(sr)):
            area_list.append(r)
            area_list += sorted(model.hierarchy.successors(r))[:5]
    area_list = pl.array(area_list)


    ## generate simulation data
    model.input_data = pandas.DataFrame(index=range(N))
    initialize_input_data(model.input_data)

    alpha = alpha_true_sim(model, area_list, sigma_true)

    # choose observed prevalence values
    model.input_data['effective_sample_size'] = ess

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

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    add_quality_metrics(model.input_data)


    model.alpha = pandas.DataFrame(index=[n for n in nx.traversal.dfs_preorder_nodes(model.hierarchy)])
    model.alpha['true'] = pandas.Series(dict(alpha))
    model.alpha['mu_pred'] = pandas.Series([n.stats()['mean'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.alpha['sigma_pred'] = pandas.Series([n.stats()['standard deviation'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    add_quality_metrics(model.alpha)

    print '\nalpha'
    print model.alpha.dropna()


    model.sigma = pandas.DataFrame(dict(true=sigma_true))
    model.sigma['mu_pred'] = [n.stats()['mean'] for n in model.vars['p']['sigma_alpha']]
    model.sigma['sigma_pred']=[n.stats()['standard deviation'] for n in model.vars['p']['sigma_alpha']]
    add_quality_metrics(model.sigma)

    print 'sigma_alpha'
    print model.sigma

    
    model.results = dict(param=[], bias=[], mare=[], mae=[], pc=[])
    add_to_results(model, 'sigma')

    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    add_quality_metrics(model.delta)

    print 'delta'
    print model.delta
    add_to_results(model, 'delta')

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())
    print 'effect prediction MAE: %.3f, coverage: %.2f' % (pl.median(pl.absolute(model.alpha['abs_err'].dropna())),
                                                          model.alpha.dropna()['covered?'].mean())
    add_to_results(model, 'input_data')
    add_to_results(model, 'alpha')

    model.results = pandas.DataFrame(model.results)
    return model


def validate_covariate_model_dispersion(N=1000, delta_true=.15, pi_true=.01, zeta_true=[.5, -.5, 0.]):
    ## generate simulated data
    a = pl.arange(0, 100, 1)
    pi_age_true = pi_true * pl.ones_like(a)

    model = data.ModelData()
    model.parameters['p']['parameter_age_mesh'] = [0, 100]
    model.input_data = pandas.DataFrame(index=range(N))
    initialize_input_data(model.input_data)

    Z = mc.rbernoulli(.5, size=(N, len(zeta_true))) * 1.0
    delta = delta_true * pl.exp(pl.dot(Z, zeta_true))
    for i in range(len(zeta_true)):
        model.input_data['z_%d'%i] = Z[:,i]

    model.input_data['true'] = pi_true

    model.input_data['effective_sample_size'] = mc.runiform(100, 10000, N)

    n = model.input_data['effective_sample_size']
    p = model.input_data['true']
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta*n*p) / n


    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=5, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)

    pl.show()


    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    add_quality_metrics(model.input_data)


    model.zeta = pandas.DataFrame(index=model.vars['p']['Z'].columns)
    model.zeta['true'] = zeta_true
    
    model.zeta['mu_pred'] = model.vars['p']['zeta'].stats()['mean']
    model.zeta['sigma_pred'] = model.vars['p']['zeta'].stats()['standard deviation']
    add_quality_metrics(model.zeta)

    print '\nzeta'
    print model.zeta
    
    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    add_quality_metrics(model.delta)

    print 'delta'
    print model.delta

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())
    print 'effect prediction MAE: %.3f, coverage: %.2f' % (pl.median(pl.absolute(model.zeta['abs_err'].dropna())),
                                                           model.zeta.dropna()['covered?'].mean())


    model.results = dict(param=[], bias=[], mare=[], mae=[], pc=[])
    add_to_results(model, 'delta')
    add_to_results(model, 'input_data')
    add_to_results(model, 'zeta')
    model.results = pandas.DataFrame(model.results, columns='param bias mae mare pc'.split())

    return model


if __name__ == '__main__':
    fe_model = validate_covariate_model_fe()
    re_model = validate_covariate_model_re()
    gnb_model = validate_covariate_model_dispersion()
