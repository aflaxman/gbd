#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

import data_model
import covariate_model
import fit_model
import graphics

reload(data_model)
reload(covariate_model)
reload(fit_model)
reload(graphics)

def fit_emp_prior(id, param_type, fast_fit=False, generate_emp_priors=True):
    """ Fit empirical prior of specified type for specified model

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    param_type : str, one of incidence, prevalence, remission, excess-mortality, prevalence_x_excess-mortality
      The disease parameter to generate empirical priors for

    Example
    -------
    >>> import fit_emp_prior
    >>> fit_emp_prior.fit_emp_prior(2552, 'incidence')
    """

    dir = dismod3.settings.JOB_WORKING_DIR % id

    ## load the model from disk or from web
    import simplejson as json
    import data
    reload(data)

    dm = dismod3.load_disease_model(id)
    
    try:
        model = data.ModelData.load(dir)
        print 'loaded data from new format from %s' % dir
    except (IOError, AssertionError):
        model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
        #model.save(dir)
        print 'loaded data from json, saved in new format for next time in %s' % dir

    ## next block fills in missing covariates with zero
    for col in model.input_data.columns:
        if col.startswith('x_'):
            model.input_data[col] = model.input_data[col].fillna(0.)
    # also fill all covariates missing in output template with zeros
    model.output_template = model.output_template.fillna(0)

    # set all heterogeneity priors to Slightly for the global fit
    for t in model.parameters:
        if 'heterogeneity' in model.parameters[t]:
            model.parameters[t]['heterogeneity'] = 'Slightly'

    t = {'incidence': 'i', 'prevalence': 'p', 'remission': 'r', 'excess-mortality': 'f', 'prevalence_x_excess-mortality': 'pf'}[param_type]
    model.input_data = model.get_data(t)
    if len(model.input_data) == 0:
        print 'No data for type %s, exiting' % param_type
        return dm

    ### For testing:
    ## speed up computation by reducing number of knots
    ## model.parameters[t]['parameter_age_mesh'] = [0, 10, 20, 40, 60, 100]

    ## smooth Slightly, Moderately, or Very
    ## model.parameters[t]['smoothness'] = dict(age_start=0, age_end=100, amount='Very')

    ## speed up computation be reducing data size
    ## predict_area = 'super-region_0'
    ## predict_year=2005
    ## predict_sex='total'
    ## subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
    ## relevant_rows = [i for i, r in model.input_data.T.iteritems() \
    ##                      if (r['area'] in subtree or r['area'] == 'all')\
    ##                      and (r['year_end'] >= 1997) \
    ##                      and r['sex'] in [predict_sex, 'total']]
    ## model.input_data = model.input_data.ix[relevant_rows]

    # testing changes
    #model.input_data['effective_sample_size'] = pl.minimum(1.e3, model.input_data['effective_sample_size'])
    #missing_ess = pl.isnan(model.input_data['effective_sample_size'])
    #model.input_data['effective_sample_size'][missing_ess] = 1.
    #model.input_data['z_overdisperse'] = 1.
    #print model.describe(t)
    #model.input_data = model.input_data[model.input_data['area'].map(lambda x: x in nx.bfs_tree(model.hierarchy, 'super-region_5'))]
    #model.input_data = model.input_data = model.input_data.drop(['x_LDI_id_Updated_7July2011'], axis=1)
    #model.input_data = model.input_data.filter([model.input_data['x_nottroponinuse'] == 0.]
    #model.input_data = model.input_data[:100]

    ## speed up output by not making predictions for empirical priors
    #generate_emp_priors = False


    print 'fitting', t
    vars = data_model.data_model(t, model, t,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None, sigma_age_parent=None,
                                 rate_type=(t == 'rr') and 'log_normal' or 'neg_binom')
    dm.model = model
    dm.vars = vars

    if fast_fit:
        dm.map, dm.mcmc = fit_model.fit_data_model(vars, iter=101, burn=0, thin=1, tune_interval=100)
    else:
        dm.map, dm.mcmc = fit_model.fit_data_model(vars, iter=50000, burn=10000, thin=40, tune_interval=1000)

    stats = dm.vars['p_pred'].stats(batches=5)
    dm.vars['data']['mu_pred'] = stats['mean']
    dm.vars['data']['sigma_pred'] = stats['standard deviation']

    stats = dm.vars['pi'].stats(batches=5)
    dm.vars['data']['mc_error'] = stats['mc error']

    dm.vars['data']['residual'] = dm.vars['data']['value'] - dm.vars['data']['mu_pred']
    dm.vars['data']['abs_residual'] = pl.absolute(dm.vars['data']['residual'])

    graphics.plot_one_type(model, vars, {}, t)
    if generate_emp_priors:
        for a in [dismod3.utils.clean(a) for a in dismod3.settings.gbd_regions]:
            print 'generating empirical prior for %s' % a
            for s in dismod3.settings.gbd_sexes:
                for y in dismod3.settings.gbd_years:
                    key = dismod3.utils.gbd_key_for(param_type, a, y, s)
                    if t in model.parameters and 'level_bounds' in model.parameters[t]:
                        lower=model.parameters[t]['level_bounds']['lower']
                        upper=model.parameters[t]['level_bounds']['upper']
                    else:
                        lower=0
                        upper=pl.inf
                    emp_priors = covariate_model.predict_for(model,
                                                             'all', 'total', 'all',
                                                             a, dismod3.utils.clean(s), int(y),
                                                             1.,
                                                             vars, lower, upper)
                    n = len(emp_priors)
                    emp_priors.sort(axis=0)
                    
                    dm.set_mcmc('emp_prior_mean', key, emp_priors.mean(0))
                    dm.set_mcmc('emp_prior_median', key, pl.median(emp_priors, axis=0))

                    design_effect = 1
                    ## uncomment to calculate a "design effect" based on over-dispersion of negative binomial
                    #if 'eta' in vars:
                    #    mu_delta = pl.exp(vars['eta'].trace()).mean()
                    #    design_effect = pl.sqrt((1+mu_delta) / mu_delta)

                    dm.set_mcmc('emp_prior_std', key, emp_priors.std(0)*design_effect)

                    pl.plot(model.parameters['ages'], dm.get_mcmc('emp_prior_mean', key), color='grey', label=a, zorder=-10, alpha=.5)
    pl.savefig(dir + '/prior-%s.png'%param_type)

    store_effect_coefficients(dm, vars, param_type)

    #graphics.plot_one_ppc(vars, t)
    #pl.savefig(dir + '/prior-%s-ppc.png'%param_type)

    graphics.plot_convergence_diag(vars)
    pl.savefig(dir + '/prior-%s-convergence.png'%param_type)
    graphics.plot_trace(vars)
    pl.savefig(dir + '/prior-%s-trace.png'%param_type)
    
    graphics.plot_one_effects(vars, t, model.hierarchy)
    pl.savefig(dir + '/prior-%s-effects.png'%param_type)

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    try:
        dm.save('dm-%d-prior-%s.json' % (id, param_type))
    except IOError, e:
        print e

    return dm



def store_effect_coefficients(dm, vars, param_type):
    """ store effect coefficients"""
    # save the results in the param_hash
    prior_vals = {}
    if isinstance(vars.get('alpha'), mc.Node):
        stats = vars['alpha'].stats()
        stats = pandas.DataFrame(dict(mean=stats['mean'], std=stats['standard deviation']), index=vars['U'].columns)
    elif isinstance(vars.get('alpha'), list):
        res = []
        for n in vars['alpha']:
            if isinstance(n, mc.Node):
                res.append(n.trace())
            else:
                res.append([n for _ in vars['p_pred'].trace()])
        stats = pl.vstack(res)
        stats = pandas.DataFrame(dict(mean=stats.mean(1), std=stats.std(1)), index=vars['U'].columns)
    else:
        stats = pandas.DataFrame(dict(mean=[], std=[]))

    prior_vals['alpha'] = [sum([0] + [stats['mean'][n] for n in nx.shortest_path(dm.model.hierarchy, 'all', dismod3.utils.clean(a)) if n in stats['mean']]) for a in dismod3.settings.gbd_regions]
    prior_vals['sigma_alpha'] = [sum([0] + [stats['std'][n] for n in nx.shortest_path(dm.model.hierarchy, 'all', dismod3.utils.clean(a)) if n in stats['mean']]) for a in dismod3.settings.gbd_regions]

    index = []
    for level in ['Country_level', 'Study_level']:
        for cv in sorted(dm.params['covariates'][level]):
            if dm.params['covariates'][level][cv]['rate']['value']:

                # do some fiddly work to get the list of covariates in the correct order
                if 'X' in vars:
                    i_list = pl.where(vars['X'].columns == 'x_%s'%cv)[0]
                else:
                    i_list = []

                if len(i_list) == 0:
                    index.append(-1)
                else:
                    index.append(i_list[0])

    if 'X_shift' in vars:
        shift =  vars['X_shift'].__array__()
    else:
        shift = 0.

    if isinstance(vars.get('beta'), list):
        stats = pl.vstack((n.trace() for n in vars['beta'])).T + shift
    else:
        stats = pl.zeros((1, max([0]+index)+1))
    stats = pandas.DataFrame(dict(mean=stats.mean(0), std=stats.std(0)))
    stats = stats.append(pandas.DataFrame(dict(mean=[0.], std=[0.]), index=[-1]))

    prior_vals['beta'] = list(stats['mean'][index])
    prior_vals['sigma_beta'] = list(stats['std'][index])


    prior_vals['new_alpha'] = {}
    if 'alpha' in vars:
        for n, col in zip(vars['alpha'], vars['U'].columns):
            if isinstance(n, mc.Node):
                stats = n.stats()
                if stats:
                    #prior_vals['new_alpha'][col] = dict(dist='TruncatedNormal', mu=stats['mean'], sigma=stats['standard deviation'], lower=-5., upper=5.)
                    prior_vals['new_alpha'][col] = dict(dist='Constant', mu=stats['mean'], sigma=stats['standard deviation'])

        # uncomment below to save empirical prior on sigma_alpha, the dispersion of the random effects
        for n in vars['sigma_alpha']:
            stats = n.stats()
            prior_vals['new_alpha'][n.__name__] = dict(dist='TruncatedNormal', mu=stats['mean'], sigma=stats['standard deviation'], lower=.01, upper=.5)

    prior_vals['new_beta'] = {}
    if 'beta' in vars:
        for n, col in zip(vars['beta'], vars['X'].columns):
            stats = n.stats()
            if stats:
                #prior_vals['new_beta'][col] = dict(dist='normal', mu=stats['mean'], sigma=stats['standard deviation'], lower=-pl.inf, upper=pl.inf)
                prior_vals['new_beta'][col] = dict(dist='Constant', mu=stats['mean'])


    if 'x_sex' in prior_vals['new_beta']:
        prior_vals['alpha'] += [0., prior_vals['new_beta']['x_sex']['mu']]
        #prior_vals['sigma_alpha'] += [0., prior_vals['new_beta']['x_sex']['sigma']]
    else:
        prior_vals['alpha'] += [0., 0.]
        #prior_vals['sigma_alpha'] += [0., 0.]


    stats = pl.log(vars['mu_age'].trace())
    prior_vals['gamma'] = list(stats.mean(0))
    prior_vals['sigma_gamma'] = list(stats.std(0))

    if 'eta' in vars:
        stats = vars['eta'].trace()
        stats = pl.exp(stats) 
        prior_vals['delta'] = stats.mean()
        prior_vals['sigma_delta'] = stats.std()

    try:
        prior_vals['aic'] = dm.map.AIC
        prior_vals['bic'] = dm.map.BIC
    except AttributeError, e:
        print 'Saving AIC/BIC failed', e

    try:
        prior_vals['dic'] = dm.mcmc.dic
    except:
        print 'Saving DIC failed'

    dm.set_empirical_prior(param_type, prior_vals)



def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--type', default='prevalence',
                      help='only estimate given parameter type (valid settings ``incidence``, ``prevalence``, ``remission``, ``excess-mortality``, ``mortality, prevalence_x_excess-mortality``) (emp prior fit only)')
    parser.add_option('-f', '--fast', default='False',
                      help='fit faster for testing')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = fit_emp_prior(id, options.type, fast_fit=(options.fast.lower()=='true'), generate_emp_priors=(options.fast.lower()=='false'))
    return dm
      

if __name__ == '__main__':
    dm = main()
    pl.show()
