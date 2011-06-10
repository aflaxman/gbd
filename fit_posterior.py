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
>>> model.fit(dm, method='mcmc', keys=keys, iter=10000, thin=5, burn=5000, verbose=1)
>>> dismod3.post_disease_model(dm)
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

from dismod3.neg_binom_model import countries_for
import dismod3.neg_binom_model as nbm
import numpy as np

import dismod3

def fit_posterior(id, region, sex, year):
    """ Fit posterior of specified region/sex/year for specified model

    Parameters
    ----------
    id : int
      The model id number for the job to fit
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years

    Example
    -------
    >>> import fit_posterior
    >>> fit_posterior.fit_posterior(2552, 'asia_east', 'male', '2005')
    """
    #print 'updating job status on server'
    #dismod3.log_job_status(id, 'posterior', '%s--%s--%s' % (region, sex, year), 'Running')

    dm = dismod3.load_disease_model(id)
    #dm.data = []  # for testing, remove all data
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    # fit the model
    dir = dismod3.settings.JOB_WORKING_DIR % id
    import dismod3.gbd_disease_model as model
    model.fit(dm, method='map', keys=keys, verbose=1)     ## first generate decent initial conditions
    ## then sample the posterior via MCMC
    model.fit(dm, method='mcmc', keys=keys, iter=100000, thin=50, burn=50000, verbose=1, dbname='/dev/null')  # dbname was '%s/posterior/pickle/dm-%d-posterior-%s-%s-%s.pickle' % (dir, id, region, sex, year)
    #model.fit(dm, method='mcmc', keys=keys, iter=1000, thin=1, burn=500, verbose=1, dbname='/dev/null') # fast, for rapid development

    # generate plots of results
    dismod3.tile_plot_disease_model(dm, keys, defaults={})
    dm.savefig('dm-%d-posterior-%s.png' % (id, '+'.join(['all', region, sex, year])))  # TODO: refactor naming into its own function (disease_json.save_image perhaps)
    for param_type in dismod3.settings.output_data_types:
        keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex], type_list=[param_type])
        dismod3.tile_plot_disease_model(dm, keys, defaults={})
        dm.savefig('dm-%d-posterior-%s-%s-%s-%s.png' % (id, dismod3.utils.clean(param_type), region, sex, year))   # TODO: refactor naming into its own function


    # summarize fit quality graphically, as well as parameter posteriors
    for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
        if dm.vars[k].get('data'):
            dismod3.plotting.plot_posterior_predicted_checks(dm, k)
            dm.savefig('dm-%d-check-%s.png' % (dm.id, k))


    # make a rate_type_list
    rate_type_list = ['incidence', 'prevalence', 'remission', 'excess-mortality', 'mortality', 'duration']
    # save country level posterior
    save_country_level_posterior(dm, region, year, sex, rate_type_list)

    # update job status file
    #print 'updating job status on server'
    #dismod3.log_job_status(id, 'posterior',
    #                       '%s--%s--%s' % (region, sex, year), 'Completed')
    
    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.save('dm-%d-posterior-%s-%s-%s.json' % (id, region, sex, year), keys_to_save=keys)

    return dm

def save_country_level_posterior(dm, region, year, sex, rate_type_list):
    """ Save country level posterior in a csv file, and put the file in the 
    directory job_working_directory/posterior/country_level_posterior_dm-'id'
    
    Parameters:
    -----------
      dm : DiseaseJson object
        disease model
      region : str
      year : str
        1990 or 2005
      sex : str
        male or female
      rate_type_list : list
        list of rate types
    """
    import csv, os
    
    import dismod3.gbd_disease_model as model
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    #dm.vars = model.setup(dm, keys)

    # get covariate dict from dm
    covariates_dict = dm.get_covariates()
    derived_covariate = dm.get_derived_covariate_values()
    
    # job working directory
    job_wd = dismod3.settings.JOB_WORKING_DIR % dm.id

    # directory to save the file
    dir = job_wd + '/posterior/'
    
    #import pymc as mc
    #picklename = 'pickle/dm-%s-posterior-%s-%s-%s.pickle' % (str(dm.id), region, sex, year)
    #model_trace = mc.database.pickle.load(dir + picklename)

    # make an output file
    filename = 'dm-%s-%s-%s-%s.csv' % (str(dm.id), region, sex, year)
    # open a file to write
    f_file = open(dir + filename, 'w')

    # get csv file writer
    csv_f = csv.writer(f_file)
    #csv_f = csv.writer(f_file, dialect=csv.excel_tab)
    print('writing csv file %s' % filename)

    # write header
    csv_f.writerow(['Iso3', 'Rate type', 'Age', 'Value', 'Lower UI', 'Upper UI'])

    # loop over countries and rate_types
    for iso3 in countries_for[region]:
        for rate_type in rate_type_list:
            # make a key
            key = '%s+%s+%s+%s' % (rate_type, region, year, dismod3.utils.clean(sex))

            # modify rate type names
            if rate_type == 'mortality':
                rate_type = 'm_with'

            # get dm.vars by the key
            model_vars = dm.vars[key]
            if rate_type == 'duration':
                # make a value_list of 0s for ages
                value_list = np.zeros((dismod3.MAX_AGE, sample_size))

                # calculate value list for ages
                for i, value_trace in enumerate(model_vars['rate_stoch'].trace()):
                    value_list[:, i] = value_trace
            else:
                # get coeffs from dm.vars
                alpha=model_vars['region_coeffs']
                beta=model_vars['study_coeffs']
                #gamma_trace = model_trace.__getattribute__('age_coeffs_%s+%s+%s+%s' % (rate_type, region, year, dismod3.utils.clean(sex))).gettrace()
                gamma_trace = model_vars['age_coeffs'].trace()

                # get sample size
                sample_size = len(gamma_trace)

                # make a value_list of 0s for ages
                value_list = np.zeros((dismod3.MAX_AGE, sample_size))

                # calculate value list for ages
                for i, gamma in enumerate(gamma_trace):
                    value_trace = nbm.predict_country_rate(key, iso3, alpha, beta, gamma,
                                                           covariates_dict, derived_covariate,
                                                           model_vars['bounds_func'],
                                                           range(101))

                    value_list[:, i] = value_trace
            if rate_type == 'prevalence':
                print key, iso3, nbm.country_covariates(key, iso3, covariates_dict, derived_covariate)[1], np.sort(value_list, axis=1)[5, .5*sample_size]
                                
            # write a row
            for age in range(dismod3.MAX_AGE):
                csv_f.writerow([iso3, rate_type, str(age)] + list(np.sort(value_list, axis=1)[age, [.5*sample_size, .025*sample_size, .975*sample_size]]))

    # close the file
    f_file.close()

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

    import time
    import random
    time.sleep(random.random()*30)  # sleep random interval before start to distribute load
    dm = fit_posterior(id, options.region, options.sex, options.year)
    return dm

if __name__ == '__main__':
    dm = main()

