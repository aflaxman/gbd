#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import sys

import pylab as pl
import pymc as mc
import networkx as nx
import pandas

import consistent_model
import data_model
import covariate_model
import fit_model
import graphics

reload(consistent_model)
reload(covariate_model)
reload(data_model)
reload(fit_model)

import dismod3

iter=10000
burn=5000
thin=5

def inspect_vars(results, vars):
    for k in vars:
        if isinstance(vars[k], mc.Node):
            d = inspect_node(vars[k])
            results.update(d)
        elif isinstance(vars[k], dict):
            inspect_vars(results, vars[k])
        elif isinstance(vars[k], list):
            inspect_vars(results, dict(zip(range(len(vars[k])), vars[k])))
    results = pandas.Series(results).order(ascending=False, na_last=False)
    return results

def inspect_node(n):
    if isinstance(n, mc.Stochastic):
        return {n.__name__: n.logp}
        #print '%s: logp=%.2f, val=%s' % (n.__name__, n.logp, n.value.round(5))
        print '%65s: logp=%.2f' % (n.__name__, n.logp)
    #elif isinstance(n, mc.Deterministic):
    #    print '%s: val=%s' % (n.__name__, n.value.round(5))
    elif isinstance(n, mc.Potential):
        return {n.__name__: n.logp}
        print '%65s: logp=%.2f' % (n.__name__, n.logp)
    else:
        return {}

def fit_posterior(dm, region, sex, year, map_only=False, 
                  inconsistent_fit=False, params_to_fit=['p', 'r', 'i']):
    """ Fit posterior of specified region/sex/year for specified model

    Parameters
    ----------
    dm : DiseaseJson
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years

    map_only : sample 101 draws from posterior, don't try for convergence (fast for testing)
    inconsistent_fit : fit parameters  separately
    params_to_fit : list of params to fit, if not fitting all consistently

    Example
    -------
    >>> import fit_posterior
    >>> fit_posterior.fit_posterior(2552, 'asia_east', 'male', '2005')
    """
    dir = dismod3.settings.JOB_WORKING_DIR % dm.id

    ## load the model from disk or from web
    import simplejson as json
    import data
    reload(data)

    try:
        assert 0
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

    predict_area = dismod3.utils.clean(region)
    predict_sex = dismod3.utils.clean(sex)
    predict_year = int(year)

    ## load emp_priors dict from dm.params
    param_type = dict(i='incidence', p='prevalence', r='remission', f='excess-mortality', rr='relative-risk', pf='prevalence_x_excess-mortality', m_with='mortality')
    emp_priors = {}
    for t in 'i r pf p rr'.split():

        # uncomment below to not use empirical prior for rate with zero data
        # if pl.all(model.input_data['data_type'] != t):
        #     continue

        #key = dismod3.utils.gbd_key_for(param_type[t], model.hierarchy.predecessors(predict_area)[0], year, sex)
        key = dismod3.utils.gbd_key_for(param_type[t], predict_area, year, sex)
        mu = dm.get_mcmc('emp_prior_mean', key)
        sigma = dm.get_mcmc('emp_prior_std', key)
        
        if len(mu) == 101 and len(sigma) == 101:
            emp_priors[t, 'mu'] = mu
            emp_priors[t, 'sigma'] = sigma

        ## update model.parameters['random_effects'] if there is information in the disease model
        expert_priors = model.parameters[t].get('random_effects', {})
        model.parameters[t]['random_effects'] = dm.get_empirical_prior(param_type[t]).get('new_alpha', {})
        model.parameters[t]['random_effects'].update(expert_priors)

        ## update model.parameters['fixed_effects'] if there is information in the disease model
        expert_fe_priors = model.parameters[t].get('fixed_effects', {})
        model.parameters[t]['fixed_effects'] = dm.get_empirical_prior(param_type[t]).get('new_beta', {})

        ## uncomment next lines to drop non-significant priors on beta effects
        ##for effect in model.parameters[t]['fixed_effects']:
        ##    prior = model.parameters[t]['fixed_effects'][effect]
        ##    if 1.96*prior['sigma'] > abs(prior['mu']):
        ##        model.parameters[t]['fixed_effects'].pop(effect)
        
        model.parameters[t]['fixed_effects'].update(expert_fe_priors)

    ## for testing might want to discard empirical priors
    ## emp_priors = {}


    ## create model and priors for region/sex/year
    # including prediction for region as empirical prior

    # select data that is about areas in this region, recent years, and sex of male or total only
    assert predict_area in model.hierarchy, 'region %s not found in area hierarchy' % predict_area
    subtree = nx.traversal.bfs_tree(model.hierarchy, predict_area)
    relevant_rows = [i for i, r in model.input_data.T.iteritems() \
                         if (r['area'] in subtree or r['area'] == 'all')\
                         and ((predict_year >= 1997 and r['year_end'] >= 1997) or
                              (predict_year <= 1997 and r['year_start'] <= 1997)) \
                         and r['sex'] in [predict_sex, 'total']]
    model.input_data = model.input_data.ix[relevant_rows]

    # replace area 'all' with predict_area
    model.input_data['area'][model.input_data['area'] == 'all'] = predict_area

    # uncomment below to drop certain data types
    # for t in 'rr smr'.split():
    #     model.input_data = model.input_data[model.input_data['data_type'] != t]

    if inconsistent_fit:
        # generate fits for requested parameters inconsistently
        vars = {}
        for t in params_to_fit:
            vars[t] = data_model.data_model(t, model, t,
                                            root_area=predict_area, root_sex=predict_sex, root_year=predict_year,
                                            mu_age=None,
                                            mu_age_parent=emp_priors.get((t, 'mu')),
                                            sigma_age_parent=emp_priors.get((t, 'sigma')),
                                            rate_type=(t == 'rr') and 'log_normal' or 'neg_binom')
            if map_only:
                fit_model.fit_data_model(vars[t], iter=101, burn=0, thin=1, tune_interval=100)
            else:
                fit_model.fit_data_model(vars[t], iter=iter, burn=burn, thin=thin, tune_interval=100)

    else:
        vars = consistent_model.consistent_model(model,
                                                 root_area=predict_area, root_sex=predict_sex, root_year=predict_year,
                                                 priors=emp_priors)

        ## fit model to data
        if map_only:
            dm.map, dm.mcmc = fit_model.fit_consistent_model(vars, 105, 0, 1, 100)
        else:
            dm.map, dm.mcmc = fit_model.fit_consistent_model(vars, iter=iter, burn=burn, thin=thin, tune_interval=100)


    # generate estimates
    posteriors = {}
    for t in 'i r f p rr pf m_with X'.split():
        if t in vars:
            if t in model.parameters and 'level_bounds' in model.parameters[t]:
                lower=model.parameters[t]['level_bounds']['lower']
                upper=model.parameters[t]['level_bounds']['upper']
            else:
                lower=0
                upper=pl.inf
            posteriors[t] = covariate_model.predict_for(model,
                                                        predict_area, predict_sex, predict_year,
                                                        predict_area, predict_sex, predict_year,
                                                        .5, # TODO: inform with het prior
                                                        vars[t], lower, upper)
    try:
        graphics.plot_fit(model, vars, emp_priors, {})
        pl.savefig(dir + '/image/posterior-%s+%s+%s.png'%(predict_area, predict_sex, predict_year))
    except Exception, e:
        print 'Error generating output graphics'
        print e

    try:
        graphics.plot_convergence_diag(vars)
        pl.savefig(dir + '/image/posterior-%s+%s+%s-convergence.png'%(predict_area, predict_sex, predict_year))
    except Exception, e:
        print 'Error generating output graphics'
        print e

    dm.vars, dm.model, dm.emp_priors = vars, model, emp_priors
    for t in 'i r f p rr pf X m_with smr'.split():
        if t not in dm.vars:
            continue
        print 'saving tables for', t
        if 'data' in dm.vars[t] and 'p_pred' in dm.vars[t]:
            stats = dm.vars[t]['p_pred'].stats(batches=5)
            dm.vars[t]['data']['mu_pred'] = stats['mean']
            dm.vars[t]['data']['sigma_pred'] = stats['standard deviation']

            stats = dm.vars[t]['pi'].stats(batches=5)
            dm.vars[t]['data']['mc_error'] = stats['mc error']

            dm.vars[t]['data']['residual'] = dm.vars[t]['data']['value'] - dm.vars[t]['data']['mu_pred']
            dm.vars[t]['data']['abs_residual'] = pl.absolute(dm.vars[t]['data']['residual'])
            if 'delta' in dm.vars[t]:
                dm.vars[t]['data']['logp'] = [mc.negative_binomial_like(n*p_obs, n*p_pred, n*p_pred*d) for n, p_obs, p_pred, d \
                                                  in zip(dm.vars[t]['data']['effective_sample_size'], dm.vars[t]['data']['value'], dm.vars[t]['data']['mu_pred'], pl.atleast_1d(dm.vars[t]['delta'].stats()['mean']))]
            dm.vars[t]['data'].to_csv(dir + '/posterior/data-%s-%s+%s+%s.csv'%(t, predict_area, predict_sex, predict_year))
        if 'U' in dm.vars[t]:
            re = dm.vars[t]['U'].T
            columns = list(re.columns)
            mu = []
            sigma = []
            for n in dm.vars[t]['alpha']:
                if isinstance(n, mc.Node):
                    mu.append(n.stats()['mean'])
                    sigma.append(n.stats()['standard deviation'])
                else:
                    mu.append(n)
                    sigma.append(0.)
            re['mu_coeff'] = mu
            re['sigma_coeff'] = sigma

            re = re.reindex(columns=['mu_coeff', 'sigma_coeff'] + columns)
            re.to_csv(dir + '/posterior/re-%s-%s+%s+%s.csv'%(t, predict_area, predict_sex, predict_year))

        if 'X' in dm.vars[t]:
            fe = dm.vars[t]['X'].T
            columns = list(fe.columns)
            mu = []
            sigma = []
            for n in dm.vars[t]['beta']:
                if isinstance(n, mc.Node):
                    mu.append(n.stats()['mean'])
                    sigma.append(n.stats()['standard deviation'])
                else:
                    mu.append(n)
                    sigma.append(0)
            fe['mu_coeff'] = mu
            fe['sigma_coeff'] = sigma

            fe = fe.reindex(columns=['mu_coeff', 'sigma_coeff'] + columns)
            fe.to_csv(dir + '/posterior/fe-%s-%s+%s+%s.csv'%(t, predict_area, predict_sex, predict_year))
                                    

    save_country_level_posterior(dm, model, vars, predict_area, predict_sex, predict_year, ['incidence', 'prevalence', 'remission'])

    keys = []
    for i, (type, long_type) in enumerate([['i', 'incidence'],
                                           ['r', 'remission'],
                                           ['f', 'excess-mortality'],
                                           ['p', 'prevalence'],
                                           ['rr', 'relative-risk'],
                                           ['pf', 'prevalence_x_excess-mortality'],
                                           ['m_with', 'mortality'],
                                           ['X', 'duration']]):
        if type in posteriors:
            n = len(posteriors[type])
            if n > 0:
                key = '%s+%s+%s+%s' % (long_type, predict_area, predict_year, predict_sex)
                keys.append(key)

                posteriors[type].sort(axis=0)
                dm.set_mcmc('mean', key, pl.mean(posteriors[type], axis=0))
                dm.set_mcmc('median', key, pl.median(posteriors[type], axis=0))
                dm.set_mcmc('lower_ui', key, posteriors[type][.025*n,:])
                dm.set_mcmc('upper_ui', key, posteriors[type][.975*n,:])

                vars = dm.vars[type]
                effects = {}
                effects['alpha'] = {}
                effects['sigma_alpha'] = {}
                if 'alpha' in vars:
                    for n, col in zip(vars['alpha'], vars['U'].columns):
                        if isinstance(n, mc.Node):
                            stats = n.stats()
                            if stats:
                                effects['alpha'][col] = dict(mu=stats['mean'], sigma=stats['standard deviation'])
                    for n in vars['sigma_alpha']:
                        stats = n.stats()
                        effects['sigma_alpha'][n.__name__] = dict(mu=stats['mean'], sigma=stats['standard deviation'])

                effects['beta'] = {}
                if 'beta' in vars:
                    for n, col in zip(vars['beta'], vars['X'].columns):
                        if isinstance(n, mc.Node):
                            stats = n.stats()
                            if stats:
                                effects['beta'][col] = dict(mu=stats['mean'], sigma=stats['standard deviation'])
                                
                dm.set_key_by_type('effects', key, effects)

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    try:
        dm.save('dm-%d-posterior-%s-%s-%s.json' % (dm.id, predict_area, predict_sex, predict_year), keys_to_save=keys)
    except IOError, e:
        print 'WARNING: could not save file'
        print e
        
    return dm

def save_country_level_posterior(dm, model, vars, region, sex, year, rate_type_list):
    """ Save country level posterior in a csv file, and put the file in the 
    directory job_working_directory/posterior/country_level_posterior_dm-'id'
    """
    import csv
    
    # job working directory
    job_wd = dismod3.settings.JOB_WORKING_DIR % dm.id

    # directory to save the file
    dir = job_wd + '/posterior/'

    for rate_type in rate_type_list:
        try:
            # make an output file
            filename = 'dm-%s-%s-%s-%s-%s.csv' % (str(dm.id), rate_type, region, sex, year)
            # open a file to write
            f_file = open(dir + filename, 'w')

            # get csv file writer
            csv_f = csv.writer(f_file)
            print('writing csv file %s' % (dir + filename))

            # write header
            csv_f.writerow(['Iso3', 'Population', 'Rate type', 'Age'] + ['Draw%d'%i for i in range(1000)])

            t = {'incidence': 'i', 'prevalence': 'p', 'remission': 'r', 'excess-mortality': 'f',
                 'duration': 'X'}[rate_type]

            if t in vars:

                # loop over countries and rate_types
                for a in model.hierarchy[region]:
                    if t in model.parameters and 'level_bounds' in model.parameters[t]:
                        lower=model.parameters[t]['level_bounds']['lower']
                        upper=model.parameters[t]['level_bounds']['upper']
                    else:
                        lower=0
                        upper=pl.inf

                    posterior = covariate_model.predict_for(model,
                                                            region, sex, year,
                                                            a, sex, year,
                                                            .5, # TODO: inform with het prior
                                                            vars[t],
                                                            lower, upper)

                    # write a row
                    pop = dismod3.neg_binom_model.population_by_age[(a, str(year), sex)]
                    ages = model.parameters['ages']
                    for i, age in enumerate(ages):
                        csv_f.writerow([a, pop[age],
                                        rate_type, str(age)] +
                                       list(posterior[:,i]))

        except IOError, e:
            print 'WARNING: could not save country level output for %s' % rate_type
            print e
def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-s', '--sex', default='male',
                      help='only estimate given sex (valid settings ``male``, ``female``, ``all``)')
    parser.add_option('-y', '--year', default='2005',
                      help='only estimate given year (valid settings ``1990``, ``2005``, ``2010``)')
    parser.add_option('-r', '--region', default='australasia',
                      help='only estimate given GBD Region')
    parser.add_option('-f', '--fast', default='False',
                      help='use MAP only')
    parser.add_option('-i', '--inconsistent', default='False',
                      help='use inconsistent model for posteriors')
    parser.add_option('-t', '--types', default='pir',
                      help='with rate types to fit (only used if inconsistent=true)')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')


    dm = dismod3.load_disease_model(id)
    dm = fit_posterior(dm, options.region, options.sex, options.year,
                       options.fast.lower() == 'true',
                       options.inconsistent.lower() == 'true',
                       options.types)
    
    return dm

if __name__ == '__main__':
    dm = main()

