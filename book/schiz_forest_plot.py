""" Fit data with several rate models and generate forest plot"""

import sys
sys.path += ['..']

import pylab as pl
import pymc as mc
import simplejson as json

import dismod3
import book_graphics
reload(book_graphics)

vars = json.load(open('schiz_forest.json'))
cy = vars['cy']
r = pl.array(vars['r'])
n = pl.array(vars['n'])

### data only plot, for computational infrastructure appendix
book_graphics.forest_plot(r, n, data_labels=cy,
                          xmax=.0115,
                          subplot_params=dict(bottom=.1, right=.99, top=.95, left=.15),
                          figparams=book_graphics.quarter_page_params,
                          fname='ci-prev_meta_analysis-schiz_data.png')


### master graphic of data and models, for rate model section of stats chapter
book_graphics.forest_plot(r, n, data_labels=cy,
                          xmax=.0115,
                          model_keys='binomial beta-binomial poisson negative-binomial normal log-normal offset-log-normal'.split(),
                          results=vars['results'],
                          #subplot_params=dict(bottom=.1, right=.99, top=.95, left=.15),
                          #figparams=book_graphics.quarter_page_params,
                          fname='theory-rate_model-schiz_forest_plot.png')
