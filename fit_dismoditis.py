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
    for j, t in enumerate('irfp'):
        pl.subplot(2, 2, j+1)
        pl.title(t)
        age_start = vars[t]['data']['age_start']
        age_end = vars[t]['data']['age_end']
        p = vars[t]['data']['value']
        for a_0i, a_1i, p_i in zip(age_start, age_end, p):
            pl.plot([a_0i, a_1i], [p_i,p_i], 'ks-', mew=1, mec='w', ms=4)

def plot_model_params(vars):
    """ 2x2 tile plot of params for consistent model"""
    for j, t in enumerate('irfp'):
        pl.subplot(2, 2, j+1)
        pl.plot(vars[t]['mu_age'].value, color=pl.cm.spectral((i+.01)/n_iter))

def demo_model_fit(n_iter=12):
    """ fit model, showing the process
    don't try to run to completion, this is just for testing
    """

    for i  in range(n_iter):
        if i < 5:
            mc.MAP(vars).fit(method='fmin_powell', verbose=1, iterlim=1)
        else:
            m = mc.MCMC(vars)
            m.sample(1001)

        plot_model_params(vars)



d = data.ModelData.from_gbd_json('tests/dismoditis.json')


# create model and priors for top level of hierarchy
root = 'all'
vars = consistent_model.consistent_model(d.input_data, d.hierarchy, root)

pl.figure()
plot_model_data(vars)
demo_model_fit()


# generate estimates for super-region_5, male, 2005
est_trace = {}
priors = {}
for t in 'irfp':
    est_trace[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'super-region_5', 'male', 2005, vars[t])
    priors[t] = est_trace[t].mean(0)


# create model and priors for (asia_southeast, male, 2005), including estimate of
# super-region_5 to borrow strength
root = 'asia_southeast'
subtree = nx.traversal.bfs_tree(d.hierarchy, root)
relevant_rows = [i for i, r in d.input_data.T.iteritems() if r['area'] in subtree and r['year_end'] >= 1997 and r['sex'] in ['male', 'total']]
vars = consistent_model.consistent_model(d.input_data.ix[relevant_rows], d.hierarchy, root, priors)

pl.figure()
plot_model_data(vars)
demo_model_fit()

# generate estimates for THA, male, 2005
posterior = {}
for t in 'irfp':
    posterior[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'THA', 'male', 2005, vars[t])

for j, t in enumerate('irfp'):
    pl.subplot(2, 2, j+1)
    pl.plot(posterior[t].mean(0), color='r', linewidth=3)
