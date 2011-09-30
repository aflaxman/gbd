""" Fit a simulation data set"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data
import data_model
import consistent_model
import covariate_model
reload(consistent_model)
reload(data_model)
reload(data)
import fit_model
import graphics

ages = pl.arange(101)
priors = {}

# load the model from disk
model = data.ModelData.from_gbd_json('/var/tmp/dismod_working/test/dm-19807/json/dm-19807.json')

model.parameters['pf'] = {}
model.parameters['pf']['increasing'] = dict(age_start=80, age_end=100)
model.parameters['pf']['decreasing'] = dict(age_start=0, age_end=0)


## create priors for top level of hierarchy by fitting i, p, rr each individually
prior_models = {}
emp_priors = {}

prior_types = 'p pf'.split()
for t in prior_types:
    print 'fitting', t
    vars = data_model.data_model('prior', model, t,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None)

    prior_models[t] = fit_model.fit_data_model(vars)
    
    # generate estimates for super-region_6, 2005, male
    emp_priors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                'all', 'total', 'all',
                                                'super-region_6', 'male', 2005, vars).mean(axis=0)
    
    graphics.plot_one_type(model, vars, emp_priors, t)

# create model and priors for (latin_america_central, male, 2005), including estimate of
# super-region_6 to borrow strength
root_area = 'latin_america_central'
subtree = nx.traversal.bfs_tree(model.hierarchy, root_area)
relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                     if r['area'] in subtree \
                     and r['year_end'] >= 1997 \
                     and r['sex'] in ['male', 'total']]
model.input_data = model.input_data.ix[relevant_rows]

vars = consistent_model.consistent_model(model, root_area=root_area, root_sex='male', root_year=2005, priors=priors)

# fit consistent model at gbd region level
posterior_model = fit_model.fit_consistent_model(vars, 30000, 15000, 15)
import pdb; pdb.set_trace()

# generate estimates for MEX, male, 2005
predict_area = 'MEX'
posteriors = {}
for t in 'i r f p rr pf'.split():
    posteriors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                root_area, 'male', 2005,
                                                predict_area, 'male', 2005, vars[t])

graphics.plot_fit(model, vars, emp_priors, posteriors)
graphics.plot_effects(vars)
graphics.plot_one_ppc(vars['pf'], 'pf')
graphics.plot_convergence_diag(vars)
