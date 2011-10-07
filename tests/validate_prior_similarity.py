""" Demonstrate different levels of prior similarity weights

Use pancreatitis model (20945) to see the effect of prior similarity
on incidence for europe_eastern, male, 2005
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
import graphics
import dismod3
import simplejson as json
import fit_model

def validate_prior_similarity():
    dm = dismod3.load_disease_model(20945)
    dm.model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
    t = 'i'
    area, sex, year = 'europe_eastern', 'male', 2005

    for het in 'Slightly Moderately Very'.split():
        dm.model.parameters[t]['parameter_age_mesh'] = [0, 15, 20, 25, 35, 45, 55, 65, 75, 100]
        dm.model.parameters[t]['heterogeneity'] = het
        setup_regional_model(dm, area, sex, year)

        dm.vars = {}
        dm.vars[t] = data_model.data_model(t, dm.model, t,
                                           root_area=area, root_sex=sex, root_year=year,
                                           mu_age=None,
                                           mu_age_parent=dm.emp_priors[t, 'mu'],
                                           sigma_age_parent=dm.emp_priors[t, 'sigma'],
                                           rate_type=(t == 'rr') and 'log_normal' or 'neg_binom')

        fit_model.fit_data_model(dm.vars[t], iter=10050, burn=5000, thin=50, tune_interval=100)

        #2graphics.plot_one_effects(dm.vars[t], t, dm.model.hierarchy)
        #pl.title(het)

        graphics.plot_convergence_diag(dm.vars[t])
        pl.title(het)

        #graphics.plot_one_ppc(dm.vars[t], t)
        #pl.title(het)

        graphics.plot_one_type(dm.model, dm.vars[t], dm.emp_priors, t)
        pl.title(het)

    pl.show()
    return dm

def setup_regional_model(dm, predict_area, predict_sex, predict_year):
    model = dm.model
    # select data that is about areas in this region, recent years, and sex of male or total only
    subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
    relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                         if (r['area'] in subtree or r['area'] == 'all')\
                         and ((predict_year == 2005 and r['year_end'] >= 1997) or r['year_start'] <= 1997) \
                         and r['sex'] in [predict_sex, 'total']]
    model.input_data = model.input_data.ix[relevant_rows]

    # replace area 'all' with predict_area
    model.input_data['area'][model.input_data['area'] == 'all'] = predict_area

    ## load emp_priors dict from dm.params
    param_type = dict(i='incidence', p='prevalence', r='remission', f='excess-mortality')
    emp_priors = {}
    for t in 'irfp':
        key = dismod3.utils.gbd_key_for(param_type[t], model.hierarchy.predecessors(predict_area)[0], predict_year, predict_sex)
        mu = dm.get_mcmc('emp_prior_mean', key)
        sigma = dm.get_mcmc('emp_prior_std', key)
        if len(mu) == 101 and len(sigma) == 101:
            emp_priors[t, 'mu'] = mu
            emp_priors[t, 'sigma'] = sigma
    dm.emp_priors = emp_priors


if __name__ == '__main__':
    dm = validate_prior_similarity()
    
