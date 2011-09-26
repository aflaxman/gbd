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

d = data.ModelData.from_gbd_json('tests/dismoditis.json')

# create model and priors
vars = consistent_model.consistent_model(d.input_data, d.hierarchy, 'all')

# plot the input data that has made it into the consistent model
plot_model_data(vars)

# fit model, showing the process
# don't try to run to completion, this is just for testing
n_iter = 8.
for i  in range(n_iter):
    if i < 3:
        mc.MAP(vars).fit(method='fmin_powell', verbose=1, iterlim=1)
    else:
        m = mc.MCMC(vars)
        m.sample(101)

    plot_model_params(vars)
