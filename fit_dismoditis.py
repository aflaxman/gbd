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


d = data.ModelData.from_gbd_json('tests/dismoditis.json')

# create model and priors
vars = consistent_model.consistent_model(d.input_data, d.hierarchy, 'all')

# fit model
for i in pl.arange(0.01, .3, .1):
    mc.MAP(vars).fit(method='fmin_powell', verbose=1, iterlim=1)
    for j, t in enumerate('irfp'):
        pl.subplot(2, 2, j)
        pl.plot(vars[t]['mu_age'].value, color=pl.cm.spectral(i))

#m = mc.MCMC(vars)
#m.use_step_method(mc.AdaptiveMetropolis, [vars[k]['gamma_bar'] for k in 'irf'] + [vars[k]['gamma'] for k in 'irf'])
#m.sample(30000, 15000, 15)
