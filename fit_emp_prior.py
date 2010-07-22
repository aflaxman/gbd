#!/usr/bin/python2.5
""" Generate empirical prior of specified parameter type

Expects the disase model json to be saved already.
"""

import simplejson as json
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
    dismod3.log_job_status(id, 'empirical_priors', param_type, 'Running')

    # load disease model
    dm = dismod3.load_disease_model(id)

    import dismod3.neg_binom_model as model
    model.fit_emp_prior(dm, param_type)

    # remove all keys that have not been changed by running this model
    keys = dismod3.utils.gbd_keys(region_list=dismod3.gbd_regions,
                    year_list=dismod3.gbd_years,
                    sex_list=dismod3.gbd_sexes)
    for k in dm.params.keys():
        if type(dm.params[k]) == dict:
            for j in dm.params[k].keys():
                if not j in keys:
                    dm.params[k].pop(j)

    dm.save('dm-%d-prior-%s.json' % (id, param_type))
    
    dismod3.log_job_status(id, 'empirical_priors', param_type, 'Completed')


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

    fit_emp_prior(id, options.type)
      

if __name__ == '__main__':
    main()
