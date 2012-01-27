""" Validate Age Integrating Random Effects Model
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import fit_model
import graphics
import data_simulation

reload(data_simulation)
reload(graphics)

pl.seterr('ignore')

def quadratic(a):
    return .0001 * (a * (100. - a) + 100.)

def validate_ai_re(N=500, delta_true=.15, sigma_true=[.1,.1,.1,.1,.1], pi_true=quadratic, smoothness='Moderately', heterogeneity='Slightly'):
    ## generate simulated data
    a = pl.arange(0, 101, 1)
    pi_age_true = pi_true(a)


    import dismod3
    import simplejson as json
    model = data.ModelData.from_gbd_jsons(json.loads(dismod3.disease_json.DiseaseJson().to_json()))
    gbd_hierarchy = model.hierarchy

    model = data_simulation.simple_model(N)
    model.hierarchy = gbd_hierarchy

    model.parameters['p']['parameter_age_mesh'] = range(0, 101, 10)
    model.parameters['p']['smoothness'] = dict(amount=smoothness)
    model.parameters['p']['heterogeneity'] = heterogeneity

    age_start = pl.array(mc.runiform(0, 100, size=N), dtype=int)
    age_end = pl.array(mc.runiform(age_start, 100, size=N), dtype=int)

    age_weights = pl.ones_like(a)
    sum_pi_wt = pl.cumsum(pi_age_true*age_weights)
    sum_wt = pl.cumsum(age_weights*1.)
    p = (sum_pi_wt[age_end] - sum_pi_wt[age_start]) / (sum_wt[age_end] - sum_wt[age_start])

    # correct cases where age_start == age_end
    i = age_start == age_end
    if pl.any(i):
        p[i] = pi_age_true[age_start[i]]

    model.input_data['age_start'] = age_start
    model.input_data['age_end'] = age_end
    model.input_data['effective_sample_size'] = mc.runiform(100, 10000, size=N)


    from validate_covariates import alpha_true_sim
    area_list = pl.array(['all', 'super-region_3', 'north_africa_middle_east', 'EGY', 'KWT', 'IRN', 'IRQ', 'JOR', 'SYR'])
    alpha = alpha_true_sim(model, area_list, sigma_true)
    print alpha

    model.input_data['true'] = pl.nan

    model.input_data['area'] = area_list[mc.rcategorical(pl.ones(len(area_list)) / float(len(area_list)), N)]
    
    for i, a in model.input_data['area'].iteritems():
        model.input_data['true'][i] = p[i] * pl.exp(pl.sum([alpha[n] for n in nx.shortest_path(model.hierarchy, 'all', a) if n in alpha]))
    p = model.input_data['true']

    n = model.input_data['effective_sample_size']
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true*n*p) / n

    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'north_africa_middle_east', 'total', 'all', None, None, None)
    #model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=1005, burn=500, thin=5, tune_interval=100)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=25, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)
    graphics.plot_one_type(model, model.vars['p'], {}, 'p')
    pl.plot(range(101), pi_age_true, 'r:', label='Truth')
    pl.legend(fancybox=True, shadow=True, loc='upper left')

    pl.show()

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    data_simulation.add_quality_metrics(model.input_data)

    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    data_simulation.add_quality_metrics(model.delta)

    model.alpha = pandas.DataFrame(index=[n for n in nx.traversal.dfs_preorder_nodes(model.hierarchy)])
    model.alpha['true'] = pandas.Series(dict(alpha))
    model.alpha['mu_pred'] = pandas.Series([n.stats()['mean'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.alpha['sigma_pred'] = pandas.Series([n.stats()['standard deviation'] for n in model.vars['p']['alpha']], index=model.vars['p']['U'].columns)
    model.alpha = model.alpha.dropna()
    data_simulation.add_quality_metrics(model.alpha)

    model.sigma = pandas.DataFrame(dict(true=sigma_true))
    model.sigma['mu_pred'] = [n.stats()['mean'] for n in model.vars['p']['sigma_alpha']]
    model.sigma['sigma_pred']=[n.stats()['standard deviation'] for n in model.vars['p']['sigma_alpha']]
    data_simulation.add_quality_metrics(model.sigma)

    print 'delta'
    print model.delta

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())

    model.mu = pandas.DataFrame(dict(true=pi_age_true,
                                     mu_pred=model.vars['p']['mu_age'].stats()['mean'],
                                     sigma_pred=model.vars['p']['mu_age'].stats()['standard deviation']))
    data_simulation.add_quality_metrics(model.mu)

    data_simulation.initialize_results(model)
    data_simulation.add_to_results(model, 'delta')
    data_simulation.add_to_results(model, 'mu')
    data_simulation.add_to_results(model, 'input_data')
    data_simulation.add_to_results(model, 'alpha')
    data_simulation.add_to_results(model, 'sigma')
    data_simulation.finalize_results(model)

    print model.results

    return model


if __name__ == '__main__':
    model = validate_age_integrating_model_sim()
