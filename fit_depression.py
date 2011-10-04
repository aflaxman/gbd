""" Fit depression data"""

# matplotlib needs to use AGG on the cluster, because X is not
# installed there
import matplotlib
matplotlib.use("AGG")

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
reload(covariate_model)
reload(data_model)
reload(data)

import dismod3

def fit_depression(dm, predict_area, predict_sex, predict_year):
    priors = {}
    predict_year = int(predict_year)
    ## load the model from disk or from web
    import simplejson as json
    model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))


    ## fill all covariates missing in output template with zeros
    model.output_template = model.output_template.fillna(0)

    ## set parameters for i, r, f
    model.parameters['i']['smoothness'] = dict(age_start=0, age_end=100, amount='Slightly')
    #model.parameters['i']['decreasing'] = dict(age_start=50, age_end=100)
    model.parameters['i']['parameter_age_mesh'] = [0,10,20,30,40,50,100]
    model.parameters['r']['parameter_age_mesh'] = [0,100]
    model.parameters['f']['parameter_age_mesh'] = [0,100]

    model.parameters['rr']['level_bound'] = dict(lower=1.5, upper=10.)


    ## these parameters are only used for (inconsistent) empirical prior
    model.parameters['p']['smoothness'] = dict(age_start=0, age_end=100, amount='Moderately')
    model.parameters['p']['parameter_age_mesh'] = [0,10,20,30,40,50,100]

    model.parameters['rr']['smoothness'] = dict(amount='Very')
    model.parameters['rr']['parameter_age_mesh'] = [0,100]


    ## uncomment next line to fit with no covariates
    # model.input_data = model.input_data.drop([col for col in model.input_data.columns if col.startswith('x_')], axis=1)

    ## next block fills in missing covariates with zero
    for col in model.input_data.columns:
        if col.startswith('x_'):
            model.input_data[col] = model.input_data[col].fillna(0.)



    ## create priors for top level of hierarchy by fitting i, p, rr each individually
    prior_models = {}
    emp_priors = {}

    prior_types = 'p rr i'.split()
    for t in prior_types:
        print 'fitting', t
        vars = data_model.data_model('prior', model, t,
                                     root_area='all', root_sex='total', root_year='all',
                                     mu_age=None, mu_age_parent=None, rate_type=(t == 'rr') and 'log_normal' or 'neg_binom')

        prior_models[t] = fit_model.fit_data_model(vars, iter=101, burn=0, thin=1)

        emp_priors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                   'all', 'total', 'all',
                                                   model.hierarchy.successors(predict_area)[0], predict_sex, predict_year, vars).mean(axis=0)
        graphics.plot_one_type(model, vars, emp_priors, t)
        if t == 'p':
            pl.plot(dm.get_mcmc('mean', 'prevalence+europe_eastern+2005+male'), 'c', linewidth=3, alpha=.5)
            graphics.plot_one_effects(vars, t, model.hierarchy)
        graphics.plot_one_ppc(vars, t)


    ## create model and priors for region/sex/year
    # including prediction for super-region as empirical prior

    # select data that is about areas in this region, recent years, and sex of male or total only
    subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
    relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                         if r['area'] in subtree \
                         and ((predict_year == 2005 and r['year_end'] >= 1997) or r['year_start'] <= 2007) \
                         and r['sex'] in [predict_sex, 'total']]
    model.input_data = model.input_data.ix[relevant_rows]

    vars = consistent_model.consistent_model(model,
                                             root_area=predict_area, root_sex=predict_sex, root_year=predict_year,
                                             priors=emp_priors)

    ## fit model to data
    posterior_model = fit_model.fit_consistent_model(vars, 1005, 500, 5, 100)


    # generate estimates
    posteriors = {}
    for t in 'i r f p rr pf'.split():
        posteriors[t] = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                    predict_area, predict_sex, predict_year,
                                                    predict_area, predict_sex, predict_year, vars[t])

    graphics.plot_fit(model, vars, emp_priors, {})
    keys = []
    for i, (type, long_type) in enumerate([['i', 'incidence'],
                                           ['r', 'remission'],
                                           ['f', 'excess-mortality'],
                                           ['p', 'prevalence'],
                                           ['rr', 'relative-risk'],
                                           ['pf', 'prevalence_x_excess-mortality']]):
        key = '%s+%s+%s+%s' % (long_type, predict_area, predict_year, predict_sex)
        keys.append(key)
        # 
        pl.subplot(2,3,i+1)
        pl.plot(dm.get_mcmc('median', key), 'c')
        pl.plot(pl.median(posteriors[type], axis=0), 'b')
        #
        posteriors[type].sort(axis=0)
        dm.set_mcmc('mean', key, pl.mean(posteriors[type], axis=0))
        dm.set_mcmc('median', key, pl.median(posteriors[type], axis=0))
        dm.set_mcmc('lower_ui', key, posteriors[type][0,:])
        dm.set_mcmc('upper_ui', key, posteriors[type][-1,:])

    dm.save('dm-%d-posterior-%s-%s-%s.json' % (dm.id, predict_area, predict_sex, predict_year), keys_to_save=keys)



def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', default='male',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year', default='2005',
                      help='only estimate given year (valid settings ``1990``, ``2005``)')
    parser.add_option('-r', '--region', default='australasia',
                      help='only estimate given GBD Region')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')


    dm = dismod3.load_disease_model(id)
    fit_depression(dm, options.region, options.sex, options.year)
    return dm

if __name__ == '__main__':
    dm = main()

