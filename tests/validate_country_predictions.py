""" Fit a model with pf and i data and a narrow band on remission"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import data
import data_model
import fit_model
import covariate_model
import graphics

reload(data_model)
reload(covariate_model)

# load the model from disk, and adjust the data and parameters for this example
model = data.ModelData.from_gbd_json('/var/tmp/dismod_working/test/dm-19807/json/dm-19807.json')

model.parameters['p']['parameter_age_mesh'] = range(0,101,20)
model.parameters['p']['fixed_effects']['x_CODcorrected_Cirrhosis_ASDR'] = dict(dist='TruncatedNormal', mu=0, sigma=.01, lower=-1., upper=1.)
model.parameters['p']['fixed_effects']['x_IHME_alcohol_liters_pc_25July11'] = dict(dist='TruncatedNormal', mu=0, sigma=.01, lower=-1., upper=1.)
# create model for global prevalence
root_area = 'all'
t = 'p'
vars = data_model.data_model(t, model, t,
                             root_area='all', root_sex='total', root_year='all',
                             mu_age=None, mu_age_parent=None, sigma_age_parent=None)

m = fit_model.fit_data_model(vars, iter=4000, burn=2000, thin=20, tune_interval=100)


# generate estimates for all leaves of the hierarchy
est = {}
for n in model.hierarchy:
    if len(model.hierarchy.successors(n)) == 0:
        est[n] = pl.median(covariate_model.predict_for(model,
                                                       'all', 'total', 'all',
                                                       n, 'male', 2005, 0., vars, 0., 1.), axis=0)

graphics.plot_one_type(model, vars, {}, 'p')

r='europe_western'
N = len(model.hierarchy[r])
for i, n in enumerate(sorted(model.hierarchy[r], key=lambda n: est[n][-1])):
    col = pl.cm.Spectral((.1+i)/N)
    pl.semilogy(est[n], label=n, color=col, linewidth=2)
    pl.text(100,est[n][-1], n, color=col)

pl.axis([-5,110,1e-5,.02])
