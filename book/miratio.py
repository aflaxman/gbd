""" Estimate ratio of IHD deaths that are due to MI in Asia South and Asia East

TODO: make this all the code that is necessary::

    import dismod3

    # load model data
    dismod3.fetch_disease_model_if_necessary(29738, 'models/miratio/')
    model = dismod3.data.load('models/miratio/')

    # select portion to fit
    model.keep(areas=['super-region_5'], year_end=1997, sexes=['male', 'total'])

    # create and fit model
    model.vars += dismod3.ism.age_specific_rate(model, data_type='p',
                                                reference_area='super-region_5', reference_sex='male', reference_year=1990)
    dismod3.fit.fit_asr(model, iter=20000, burn=5000, thin=15, tune_interval=100)

    # use results of fit
    model.predict_for(data_type='p', area='asia_south', year=1990, sex='male')
    dismod3.graphics.summarize_fit(model)
    dismod3.upload_predictions(29738)
"""

predict_area='super-region_5'
predict_sex='female'
predict_year=1990


# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc

import dismod3
reload(dismod3)


## load MI Ratio model data
id = 29738
dir = dismod3.settings.JOB_WORKING_DIR % id
model = dismod3.data.ModelData.load(dir)  # TODO: move load function up one level, so the one line command is model = dismod3.data.load('models/miratio/')




## uncomment code to download model from web and save it locally
# TODO: make a dismod3 util that takes care of this, and also starts a git repository for the data
# TODO: make another util that merges in the covariates from the CODEm database
#dm = dismod3.load_disease_model(id)
#model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
#model.save(dir)



# examples of changing parameters programatically (for sensitivity
# analysis, etc)
# TODO: make dismod3 util functions to change all of these things
#model.parameters['p']['parameter_age_mesh'] = range(0, 101, 10)
#model.parameters['p']['smoothness'] = dict(amount='Slightly')


# select only the relevant data
# TODO: make this as simple as::
#
#     model.input_data = model.input_data_for(area='super-region_5', year_end=1997, sexes=['male', 'total'])


# filter out rows of data that are not relevant to the predictions
import networkx as nx
subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                     if (r['area'] in subtree or r['area'] == 'all')\
                     and ((predict_year >= 1997 and r['year_end'] >= 1997) or
                          (predict_year <= 1997 and r['year_start'] <= 1997)) \
                     and r['sex'] in [predict_sex, 'total']]


model.input_data = model.input_data.ix[relevant_rows]


# build a model of this data
# TODO: rethink name of "data_model", e.g. "rate_model", but then what should rate model be called?
# TODO: simplify method parameters, there is some redundancy here
t = 'p'
model.vars += dismod3.data_model.data_model(t, model, t,
                                    root_area=predict_area, root_sex=predict_sex, root_year=predict_year,
                                    mu_age=None, mu_age_parent=None, sigma_age_parent=None,
                                    rate_type='neg_binom')


# fit model of data with MCMC, using MAP for initial values
model.map, model.mcmc = dismod3.fit_model.fit_data_model(model.vars, iter=20000, burn=5000, thin=15, tune_interval=100)


# display results, and be sure to check MCMC convergence
dismod3.graphics.plot_one_type(model, model.vars, {}, t)
dismod3.graphics.plot_one_effects(model.vars, 'p', model.hierarchy)
model.vars.plot_acorr()

pl.show()
