""" Validate Consistent Model
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import networkx as nx
import pymc as mc

import data
import consistent_model
import fit_model
import graphics
import data_simulation
import pandas

reload(consistent_model)
reload(fit_model)
reload(data_simulation)
reload(graphics)

pl.seterr('ignore')

def quadratic(a):
    return 1e-6 * (a * (100. - a) + 100.)

def constant(a):
    return .2 * pl.ones_like(a)

def validate_consistent_re(N=500, delta_true=.15, sigma_true=[.1,.1,.1,.1,.1], 
                           true=dict(i=quadratic, f=constant, r=constant)):
    types = pl.array(['i', 'r', 'f', 'p'])

    ## generate simulated data
    model = data_simulation.simple_model(N)
    model.input_data['effective_sample_size'] = 1.
    model.input_data['value'] = 0.
    # coarse knot spacing for fast testing
    for t in types:
        model.parameters[t]['parameter_age_mesh'] = range(0, 101, 20)

    sim = consistent_model.consistent_model(model, 'all', 'total', 'all', {})
    for t in 'irf':
        for i, k_i in enumerate(sim[t]['knots']):
            sim[t]['gamma'][i].value = pl.log(true[t](k_i))

    age_start = pl.array(mc.runiform(0, 100, size=N), dtype=int)
    age_end = pl.array(mc.runiform(age_start, 100, size=N), dtype=int)

    data_type = types[mc.rcategorical(pl.ones(len(types), dtype=float) / float(len(types)), size=N)]


    a = pl.arange(101)
    age_weights = pl.ones_like(a)
    sum_wt = pl.cumsum(age_weights)

    p = pl.zeros(N)
    for t in types:
        mu_t = sim[t]['mu_age'].value
        sum_mu_wt = pl.cumsum(mu_t*age_weights)
    
        p_t = (sum_mu_wt[age_end] - sum_mu_wt[age_start]) / (sum_wt[age_end] - sum_wt[age_start])

        # correct cases where age_start == age_end
        i = age_start == age_end
        if pl.any(i):
            p_t[i] = mu_t[age_start[i]]

        # copy part into p
        p[data_type==t] = p_t[data_type==t]


    # add covariate shifts
    import dismod3
    import simplejson as json
    gbd_model = data.ModelData.from_gbd_jsons(json.loads(dismod3.disease_json.DiseaseJson().to_json()))
    model.hierarchy = gbd_model.hierarchy

    from validate_covariates import alpha_true_sim
    area_list = pl.array(['all', 'super-region_3', 'north_africa_middle_east', 'EGY', 'KWT', 'IRN', 'IRQ', 'JOR', 'SYR'])
    alpha = {}
    for t in types:
        alpha[t] = alpha_true_sim(model, area_list, sigma_true)
    print json.dumps(alpha, indent=2)

    model.input_data['area'] = area_list[mc.rcategorical(pl.ones(len(area_list)) / float(len(area_list)), N)]
    
    for i, a in model.input_data['area'].iteritems():
        t = data_type[i]
        p[i] = p[i] * pl.exp(pl.sum([alpha[t][n] for n in nx.shortest_path(model.hierarchy, 'all', a) if n in alpha]))

    n = mc.runiform(100, 10000, size=N)

    model.input_data['data_type'] = data_type
    model.input_data['age_start'] = age_start
    model.input_data['age_end'] = age_end
    model.input_data['effective_sample_size'] = n
    model.input_data['true'] = p
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true) / n

    # coarse knot spacing for fast testing
    for t in types:
        model.parameters[t]['parameter_age_mesh'] = range(0, 101, 20)

    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars = consistent_model.consistent_model(model, 'all', 'total', 'all', {})
    #model.map, model.mcmc = fit_model.fit_consistent_model(model.vars, iter=101, burn=0, thin=1, tune_interval=100)
    model.map, model.mcmc = fit_model.fit_consistent_model(model.vars, iter=10000, burn=5000, thin=25, tune_interval=100)

    graphics.plot_convergence_diag(model.vars)

    graphics.plot_fit(model, model.vars, {}, {})
    for i, t in enumerate('i r f p rr pf'.split()):
        pl.subplot(2, 3, i+1)
        pl.plot(range(101), sim[t]['mu_age'].value, 'w-', label='Truth', linewidth=2)
        pl.plot(range(101), sim[t]['mu_age'].value, 'r-', label='Truth', linewidth=1)

    pl.show()

    model.input_data['mu_pred'] = 0.
    model.input_data['sigma_pred'] = 0.
    for t in types:
        model.input_data['mu_pred'][data_type==t] = model.vars[t]['p_pred'].stats()['mean']
        model.input_data['sigma_pred'][data_type==t] = model.vars[t]['p_pred'].stats()['standard deviation']
    data_simulation.add_quality_metrics(model.input_data)

    model.delta = pandas.DataFrame(dict(true=[delta_true for t in types if t != 'rr']))
    model.delta['mu_pred'] = [pl.exp(model.vars[t]['eta'].trace()).mean() for t in types if t != 'rr']
    model.delta['sigma_pred'] = [pl.exp(model.vars[t]['eta'].trace()).std() for t in types if t != 'rr']
    data_simulation.add_quality_metrics(model.delta)

    model.alpha = pandas.DataFrame()
    model.sigma = pandas.DataFrame()
    for t in types:
        alpha_t = pandas.DataFrame(index=[n for n in nx.traversal.dfs_preorder_nodes(model.hierarchy)])
        alpha_t['true'] = pandas.Series(dict(alpha[t]))
        alpha_t['mu_pred'] = pandas.Series([n.stats()['mean'] for n in model.vars[t]['alpha']], index=model.vars[t]['U'].columns)
        alpha_t['sigma_pred'] = pandas.Series([n.stats()['standard deviation'] for n in model.vars[t]['alpha']], index=model.vars[t]['U'].columns)
        alpha_t['type'] = t
        model.alpha = model.alpha.append(alpha_t.dropna(), ignore_index=True)

        sigma_t = pandas.DataFrame(dict(true=sigma_true))
        sigma_t['mu_pred'] = [n.stats()['mean'] for n in model.vars[t]['sigma_alpha']]
        sigma_t['sigma_pred'] = [n.stats()['standard deviation'] for n in model.vars[t]['sigma_alpha']]
        model.sigma = model.sigma.append(sigma_t.dropna(), ignore_index=True)

    data_simulation.add_quality_metrics(model.alpha)
    data_simulation.add_quality_metrics(model.sigma)


    print 'delta'
    print model.delta

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())

    model.mu = pandas.DataFrame()
    for t in types:
        model.mu = model.mu.append(pandas.DataFrame(dict(true=sim[t]['mu_age'].value,
                                                         mu_pred=model.vars[t]['mu_age'].stats()['mean'],
                                                         sigma_pred=model.vars[t]['mu_age'].stats()['standard deviation'])),
                                   ignore_index=True)
    data_simulation.add_quality_metrics(model.mu)
    print '\nparam prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.mu['abs_err'].mean(),
                                                                         pl.median(pl.absolute(model.mu['rel_err'].dropna())),
                                                                         model.mu['covered?'].mean())
    print


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
    model = validate_consistent_re()
