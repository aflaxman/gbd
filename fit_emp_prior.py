#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3

import pylab as pl

import data_model
import covariate_model
import fit_model
import graphics

reload(data_model)
reload(covariate_model)
reload(fit_model)
reload(graphics)

def fit_emp_prior(id, param_type, map_only=False):
    """ Fit empirical prior of specified type for specified model

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    param_type : str, one of incidence, prevalence, remission, excess-mortality
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
    model = data.ModelData.from_gbd_jsons(json.loads(dm.to_json()))

    ## next block fills in missing covariates with zero
    for col in model.input_data.columns:
        if col.startswith('x_'):
            model.input_data[col] = model.input_data[col].fillna(0.)
    # also fill all covariates missing in output template with zeros
    model.output_template = model.output_template.fillna(0)

    t = {'incidence': 'i', 'prevalence': 'p', 'remission': 'r', 'excess-mortality': 'f'}[param_type]
    data = model.get_data(t)
    if len(data) == 0:
        print 'No data for type %s, exiting' % param_type
        return dm

    print 'fitting', t
    vars = data_model.data_model('prior', model, t,
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None, rate_type=(t == 'rr') and 'log_normal' or 'neg_binom')
    if map_only:
        fit_model.fit_data_model(vars, iter=101, burn=0, thin=1, tune_interval=100)
    else:
        fit_model.fit_data_model(vars, iter=10050, burn=5000, thin=50, tune_interval=100)


    graphics.plot_one_type(model, vars, {}, t)
    for a in model.hierarchy['all']:
        print 'generating empirical prior for %s' % a
        for s in dismod3.settings.gbd_sexes:
            for y in dismod3.settings.gbd_years:
                key = dismod3.utils.gbd_key_for(param_type, a, y, s)
                emp_priors = covariate_model.predict_for(model.output_template, model.hierarchy,
                                                         'all', 'total', 'all',
                                                         a, dismod3.utils.clean(s), int(y),
                                                         vars)
                n = len(emp_priors)
                emp_priors.sort(axis=0)
                dm.set_mcmc('emp_prior_mean', key, emp_priors.mean(0))
                dm.set_mcmc('emp_prior_std', key, emp_priors.std(0))
    
                pl.plot(model.parameters['ages'], dm.get_mcmc('emp_prior_mean', key), 'r-')

    pl.savefig(dir + '/prior-%s.png'%param_type)

    graphics.plot_one_ppc(vars, t)
    pl.savefig(dir + '/prior-%s-ppc.png'%param_type)

    graphics.plot_convergence_diag(vars)
    pl.savefig(dir + '/prior-%s-convergence.png'%param_type)
        
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    try:
        dm.save('dm-%d-prior-%s.json' % (id, param_type))
    except IOError, e:
        print e

    return dm

def main():
    import optparse

    usage = 'usage: %prog [options] disease_model_id'
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--type', default='prevalence',
                      help='only estimate given parameter type (valid settings ``incidence``, ``prevalence``, ``remission``, ``excess-mortality``) (emp prior fit only)')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('incorrect number of arguments')

    try:
        id = int(args[0])
    except ValueError:
        parser.error('disease_model_id must be an integer')

    dm = fit_emp_prior(id, options.type)
    return dm
      

if __name__ == '__main__':
    dm = main()
