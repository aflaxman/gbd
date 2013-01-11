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

import ism
import fit
import covariate_model
import fit_model
import graphics

reload(fit_model)


def fit_world(id, fast_fit=False, zero_re=True, alt_prior=False, global_heterogeneity='Slightly'):
    """ Fit consistent for all data in world

    Parameters
    ----------
    id : int
      The model id number for the job to fit

    Example
    -------
    >>> import fit_world
    >>> dm = fit_world.dismod3.load_disease_model(1234)
    >>> fit_world.fit_world(dm)
    """

    dir = dismod3.settings.JOB_WORKING_DIR % id

    ## load the model from disk or from web
    import simplejson as json
    import data
    reload(data)

    try:
        model = data.ModelData.load(dir)
        print 'loaded data from new format from %s' % dir
        dm = dismod3.load_disease_model(id)
    except (IOError, AssertionError):
        dm = dismod3.load_disease_model(id)
        model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))
        try:
            model.save(dir)
            print 'loaded data from json, saved in new format for next time in %s' % dir
        except IOError:
            print 'loaded data from json, failed to save in new format'


    ## next block fills in missing covariates with zero
    for col in model.input_data.columns:
        if col.startswith('x_'):
            model.input_data[col] = model.input_data[col].fillna(0.)
    # also fill all covariates missing in output template with zeros
    model.output_template = model.output_template.fillna(0)

    # set all heterogeneity priors to Slightly for the global fit
    for t in model.parameters:
        if 'heterogeneity' in model.parameters[t]:
            model.parameters[t]['heterogeneity'] = global_heterogeneity

    ### For testing:
    ## speed up computation by reducing number of knots
    ## for t in 'irf':
    ##     model.parameters[t]['parameter_age_mesh'] = [0, 100]
    model.vars += dismod3.ism.consistent(model,
                                         reference_area='all',
                                         reference_sex='total',
                                         reference_year='all',
                                         priors={},
                                         zero_re=zero_re)

    ## fit model to data
    if fast_fit:
        dm.map, dm.mcmc = dismod3.fit.fit_consistent(model, 105, 0, 1, 100)
    else:
        dm.map, dm.mcmc = dismod3.fit.fit_consistent(model, iter=50000, burn=10000, thin=40, tune_interval=1000, verbose=True)

    dm.model = model

    # borrow strength to inform sigma_alpha between rate types post-hoc
    types_with_re = ['rr', 'f', 'i', 'm', 'smr', 'p', 'r', 'pf', 'm_with', 'X']
    ## first calculate sigma_alpha_bar from posterior draws from each alpha
    alpha_vals = []
    for type in types_with_re:
        if 'alpha' in model.vars[type]:
            for alpha_i in model.vars[type]['alpha']:
                alpha_vals += [a for a in alpha_i.trace() if a != 0]  # remove zeros because areas with no siblings are included for convenience but are pinned to zero
    ## then blend sigma_alpha_i and sigma_alpha_bar for each sigma_alpha_i
    if len(alpha_vals) > 0:
        sigma_alpha_bar = pl.std(alpha_vals)
        for type in types_with_re:
            if 'sigma_alpha' in model.vars[type]:
                for sigma_alpha_i in model.vars[type]['sigma_alpha']:
                    cur_val = sigma_alpha_i.trace()
                    sigma_alpha_i.trace._trace[0] = (cur_val + sigma_alpha_bar) * pl.ones_like(sigma_alpha_i.trace._trace[0])


    for t in 'p i r f rr pf m_with'.split():
        param_type = dict(i='incidence', r='remission', f='excess-mortality', p='prevalence', rr='relative-risk', pf='prevalence_x_excess-mortality', m_with='mortality')[t]
        #graphics.plot_one_type(model, model.vars[t], {}, t)
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
                                                             model.parameters.get(t, {}),
                                                             'all', 'total', 'all',
                                                             a, dismod3.utils.clean(s), int(y),
                                                             alt_prior,
                                                             model.vars[t], lower, upper)
                    dm.set_mcmc('emp_prior_mean', key, emp_priors.mean(0))
                    if 'eta' in model.vars[t]:
                        N,A = emp_priors.shape  # N samples, for A age groups
                        delta_trace = pl.transpose([pl.exp(model.vars[t]['eta'].trace()) for _ in range(A)])  # shape delta matrix to match prediction matrix
                        emp_prior_std = pl.sqrt(emp_priors.var(0) + (emp_priors**2 / delta_trace).mean(0))
                    else:
                        emp_prior_std = emp_priors.std(0)
                    dm.set_mcmc('emp_prior_std', key, emp_prior_std)


        from fit_emp_prior import store_effect_coefficients
        store_effect_coefficients(dm, model.vars[t], param_type)

    
        if 'p_pred' in model.vars[t]:
            graphics.plot_one_ppc(model, t)
            pl.savefig(dir + '/prior-%s-ppc.png'%param_type)

        if 'p_pred' in model.vars[t] or 'lb' in model.vars[t]:
            graphics.plot_one_effects(model, t)
            pl.savefig(dir + '/prior-%s-effects.png'%param_type)


    for t in 'i r f p rr pf X m_with smr'.split():
        fname = dir + '/empirical_priors/data-%s.csv'%t
        print 'saving tables for', t, 'to', fname
        if 'data' in model.vars[t] and 'p_pred' in model.vars[t]:
            stats = model.vars[t]['p_pred'].stats(batches=5)
            model.vars[t]['data']['mu_pred'] = stats['mean']
            model.vars[t]['data']['sigma_pred'] = stats['standard deviation']

            stats = model.vars[t]['pi'].stats(batches=5)
            model.vars[t]['data']['mc_error'] = stats['mc error']

            model.vars[t]['data']['residual'] = model.vars[t]['data']['value'] - model.vars[t]['data']['mu_pred']
            model.vars[t]['data']['abs_residual'] = pl.absolute(model.vars[t]['data']['residual'])
            #if 'delta' in model.vars[t]:
            #    model.vars[t]['data']['logp'] = [mc.negative_binomial_like(n*p_obs, n*p_pred, n*p_pred*d) for n, p_obs, p_pred, d \
            #                                  in zip(model.vars[t]['data']['effective_sample_size'],
            #                                         model.vars[t]['data']['value'],
            #                                         model.vars[t]['data']['mu_pred'],
            #                                         pl.atleast_1d(model.vars[t]['delta'].stats()['mean']))]
            model.vars[t]['data'].to_csv(fname)


    graphics.plot_fit(model)
    pl.savefig(dir + '/prior.png')

    graphics.plot_acorr(model)
    pl.savefig(dir + '/prior-convergence.png')

    graphics.plot_trace(model)
    pl.savefig(dir + '/prior-trace.png')
    
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    try:
        dm.save('dm-%d-prior-%s.json' % (dm.id, 'all'))
    except IOError, e:
        print e

    return dm

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-f', '--fast', default='False',
                      help='use MAP only')
    parser.add_option('-z', '--zerore', default='true',
                      help='enforce zero constraint on random effects')
    parser.add_option('-a', '--altprior', default='false',
                      help='use alternative aggregation for empirical prior')
    parser.add_option('-g', '--globalheterogeneity', default='Slightly',
                      help='negative binomial heterogeneity for global estimate')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = fit_world(id, options.fast.lower() == 'true',
                   zero_re=options.zerore.lower() == 'true',
                   alt_prior=options.altprior.lower() == 'true',
                   global_heterogeneity=options.globalheterogeneity)
    return dm
      

if __name__ == '__main__':
    dm = main()
