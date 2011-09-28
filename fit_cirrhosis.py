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
reload(data)

from fit_dismoditis import *

# load the model from disk
d = data.ModelData.from_gbd_json('/var/tmp/dismod_working/test/dm-19807/json/dm-19807.json')


# create model and priors for top level of hierarchy
root = 'all'
vars = consistent_model.consistent_model(d.input_data, d.parameters, d.hierarchy, root)

pl.figure()
plot_model_data(vars)
m1 = fit_model(vars)

# generate estimates for super-region_6, male, 2005
est_trace = {}
priors = {}
for t in 'irfp':
    est_trace[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'super-region_6', 'male', 2005, vars[t])
    priors[t] = est_trace[t].mean(0)

for j, t in enumerate('irfp'):
    pl.subplot(2, 2, j+1)
    pl.plot(ages, priors[t], color='r', linewidth=1)


# create model and priors for (latin_america_central, male, 2005), including estimate of
# super-region_5 to borrow strength
root = 'latin_america_central'
subtree = nx.traversal.bfs_tree(d.hierarchy, root)
relevant_rows = [i for i, r in d.input_data.T.iteritems() if r['area'] in subtree and r['year_end'] >= 1997 and r['sex'] in ['male', 'total']]
vars = consistent_model.consistent_model(d.input_data.ix[relevant_rows], d.parameters, d.hierarchy, root, priors=priors, ages=ages)

# fit initial conditions to data
mc.MAP([vars['logit_C0'], vars['p']]).fit(tol=.01, verbose=1)

pl.figure()
plot_model_data(vars)
for j, t in enumerate('irfp'):
    pl.subplot(2, 2, j+1)
    pl.plot(ages, priors[t], color='r', linewidth=1)


m2 = fit_model(vars)


# generate estimates for MEX, male, 2005
posteriors = {}
for t in 'irfp':
    posteriors[t] = data_model.predict_for(d.output_template, d.hierarchy, root, 'MEX', 'male', 2005, vars[t])

for j, t in enumerate('irfp'):
    pl.subplot(2, 2, j+1)
    pl.plot(ages, posteriors[t].mean(0), color='b', linewidth=1)

