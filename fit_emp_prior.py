#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3

def fit_emp_prior(id, param_type):
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
    #dismod3.log_job_status(id, 'empirical_priors', param_type, 'Running')

    # load disease model
    dm = dismod3.load_disease_model(id)
    #dm.data = []  # remove all data to speed up computation, for test

    import dismod3.neg_binom_model as model
    dir = dismod3.settings.JOB_WORKING_DIR % id
    model.fit_emp_prior(dm, param_type, '%s/empirical_priors/pickle/dm-%d-emp_prior-%s.pickle' % (dir, id, param_type))

    # generate empirical prior plots
    from pylab import subplot
    for sex in dismod3.settings.gbd_sexes:
        for year in dismod3.settings.gbd_years:
            keys = dismod3.utils.gbd_keys(region_list=['all'], year_list=[year], sex_list=[sex], type_list=[param_type])
            dismod3.tile_plot_disease_model(dm, keys, defaults={})
            dm.savefig('dm-%d-emp_prior-%s-%s-%s.png' % (id, param_type, sex, year))

    # TODO: put this in a separate script, which runs after all empirical priors are computed
    for effect in ['alpha', 'beta', 'gamma', 'delta']:
        dismod3.plotting.plot_empirical_prior_effects([dm], effect)
        dm.savefig('dm-%d-emp-prior-%s-%s.png' % (id, param_type, effect))

    # summarize fit quality graphically, as well as parameter posteriors
    k0 = keys[0]
    dm.vars = {k0: dm.vars}   # hack to make posterior predictions plot
    dismod3.plotting.plot_posterior_predicted_checks(dm, k0)
    dm.savefig('dm-%d-emp-prior-check-%s.png' % (dm.id, param_type))
    dm.vars = dm.vars[k0]   # undo hack to make posterior predictions plot
    
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    dm.save('dm-%d-prior-%s.json' % (id, param_type))
    #dismod3.log_job_status(id, 'empirical_priors', param_type, 'Completed')
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
