#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year

Examples
--------
$ python fit_posterior.py 3828 -t incidence   # command-line example

>>> # ipython example
>>> from fit_emp_prior import *
>>> dm = dismod3.get_disease_model(3828)
>>> import dismod3.neg_binom_model as model
>>> model.fit_emp_prior(dm, 'incidence')
>>> dismod3.post_disease_model(dm)
"""

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key
from dismod3.plotting import GBDDataHash

def fit_emp_prior(id, param_type):
    dismod3.log_job_status(id, 'empirical_priors', param_type, 'Running')

    dm = dismod3.get_disease_model(id)

    import dismod3.neg_binom_model as model
    model.fit_emp_prior(dm, param_type)

    # remove all keys that have not been changed by running this model
    keys = gbd_keys(region_list=dismod3.gbd_regions,
                    year_list=dismod3.gbd_years,
                    sex_list=dismod3.gbd_sexes)
    for k in dm.params.keys():
        if type(dm.params[k]) == dict:
            for j in dm.params[k].keys():
                if not j in keys:
                    dm.params[k].pop(j)

    # post results to dismod_data_server
    # "dumb" error handling, in case post fails (try: except: sleep random time, try again, stop after 3 tries)
    from twill.errors import TwillAssertionError
    import random
    import time

    for ii in range(3):
        try:
            url = dismod3.post_disease_model(dm)
        except TwillAssertionError:
            pass

        time.sleep(random.random()*30)

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
