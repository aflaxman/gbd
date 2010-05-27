#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year

Examples
--------

$ python fit_posterior.py 3828 -r australasia -s male -y 2005 

>>> # ipython example
>>> from fit_posterior import *
>>> dm = dismod3.get_disease_model(3828)
>>> mort = dismod3.get_disease_model('all-cause_mortality')
>>> dm.data += mort.data
>>> import dismod3.gbd_disease_model as model
>>> model.fit(dm, method='map', keys=keys)
>>> model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)
>>> dismod3.post_disease_model(dm)
"""

import dismod3
from dismod3.utils import clean, gbd_keys, type_region_year_sex_from_key
from dismod3.plotting import GBDDataHash

def fit_posterior(id, region, sex, year):
    dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (region, sex, year), 'Running')

    dm = dismod3.get_disease_model(id)

    keys = gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    import dismod3.gbd_disease_model as model

    # get the all-cause mortality data, and merge it into the model
    mort = dismod3.get_disease_model('all-cause_mortality')
    dm.data += mort.data

    dm.params['estimate_type'] = 'fit individually'

    # fit the model
    ## first generate decent initial conditions
    model.fit(dm, method='map', keys=keys, verbose=1)
    ## then sample the posterior via MCMC
    model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=5, burn=5000, verbose=1)
    #model.fit(dm, method='mcmc', keys=keys, iter=10, thin=5, burn=50, verbose=1) # quick version, for debugging

    # remove all keys that have not been changed by running this model
    # this prevents overwriting estimates that are being generated simulatneously
    # by other nodes in a cluster
    for k in dm.params.keys():
        if type(dm.params[k]) == dict:
            for j in dm.params[k].keys():
                if not j in keys:
                    dm.params[k].pop(j)

    # post results to dismod_data_server
    #url = dismod3.post_disease_model(dm) # the good way, commented out for experimental error handling

    # "dumb" error handling, in case post fails (try: except: sleep random time, try again, stop after 3 tries)
    from twill.errors import TwillAssertionError
    import random

    for ii in range(3):
        try:
            url = dismod3.post_disease_model(dm)
        except TwillAssertionError:
            pass
        import time
        time.sleep(random.random()*30)

    # update job status file
    dismod3.log_job_status(id, 'posterior',
                           '%s--%s--%s' % (region, sex, year), 'Completed')


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

    fit_posterior(id, options.region, options.sex, options.year)
        

if __name__ == '__main__':
    main()
