""" Fit moderate vision loss data set"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import covariate_model
import consistent_model
import fit_model
import graphics

reload(fit_model)
reload(consistent_model)
reload(data_model)
reload(data)

ages = pl.arange(101)
priors = {}

## load the model from disk or from web
import dismod3
import simplejson as json
id = 20525
dm = dismod3.load_disease_model(id)
model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))

## create priors for top level of hierarchy by fitting i, p, rr each individually
prior_models = {}
emp_priors = {}

prior_types = 'i p rr'.split()
for t in prior_types:
    print 'fitting', t
    vars = data_model.data_model('prior', model, t,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None)

    prior_models[t] = fit_model.fit_data_model(vars)
    
    emp_priors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                'all', 'total', 'all',
                                                'super-region_5', 'male', 2005, vars).mean(axis=0)
    
    graphics.all_plots_for(model, vars, emp_priors, t)


## create model and priors for (latin_america_central, male, 2005)
# including prediction for super-region as empirical prior
root_area = 'asia_southeast'

# select data that is about areas in this region, recent years, and sex of male or total only
subtree = nx.traversal.bfs_tree(model.hierarchy, root_area)

# TODO: add .iterrows method to Pandas (or is there an iteritems(axis=1)?)
relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                     if r['area'] in subtree \
                     and r['year_end'] >= 1997 \
                     and r['sex'] in ['male', 'total']]
model.input_data = model.input_data.ix[relevant_rows]

vars = consistent_model.consistent_model(model,
                                         root_area=root_area, root_sex='male', root_year=2005,
                                         priors=emp_priors)

# fit model to data
posterior_model = fit_model.fit_consistent_model(vars, 40000, 10000, 200)

# generate estimates for THA, 2005, male
predict_area = 'THA'
posteriors = {}
for t in 'i r f p rr pf'.split():
    posteriors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                root_area, 'male', 2005,
                                                predict_area, 'male', 2005, vars[t]).mean(axis=0)

graphics.all_plots(model, vars, emp_priors, posteriors)
