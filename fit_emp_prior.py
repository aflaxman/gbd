#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import dismod3
import Matplot

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
    # load disease model
    dm = dismod3.load_disease_model(id)

    dir = dismod3.settings.JOB_WORKING_DIR % id
    dismod3.neg_binom_model.fit_emp_prior(dm, param_type, dbname='%s/empirical_priors/pickle/dm-%d-emp_prior-%s.pickle' % (dir, id, param_type))

    # summarize fit quality graphically, as well as parameter posteriors
    # skip saving images if there is no relevant data
    if not hasattr(dm, 'vars'):
        return dm
    try:
        dm.vars = {param_type: dm.vars}  # dm.vars dict is a hack to make posterior predictions plot
        dismod3.plotting.plot_posterior_predicted_checks(dm, param_type)
        dm.savefig('dm-%d-emp-prior-check-%s.png' % (dm.id, param_type))
        dm.vars = dm.vars[param_type]   # undo hack to make posterior predictions plot
    except ValueError, e:
        print e
        
    import pymc as mc
    path = '%s/image/'%dir
    try:
        Matplot.plot(dm.vars['dispersion'], path=path)
        Matplot.plot(dm.vars['age_coeffs_mesh'], path=path)
        Matplot.plot(dm.vars['study_coeffs'], path=path)
        Matplot.plot(dm.vars['region_coeffs'], path=path)
    except e:
        print e
        
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    dm.save('dm-%d-prior-%s.json' % (id, param_type))
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
