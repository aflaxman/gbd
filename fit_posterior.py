#!/usr/bin/python2.5
""" Generate a posterior estimate for a specific region, sex, and year
"""

# matplotlib backend setup
import matplotlib
matplotlib.use("AGG") 

import sys

import pylab as pl
import pymc as mc
import Matplot

import dismod3

def fit_posterior(dm, region, sex, year):
    """ Fit posterior of specified region/sex/year for specified model

    Parameters
    ----------
    dm : DiseaseJson
    region : str
      From dismod3.settings.gbd_regions, but clean()-ed
    sex : str, from dismod3.settings.gbd_sexes
    year : str, from dismod3.settings.gbd_years

    Example
    -------
    >>> import fit_posterior
    >>> fit_posterior.fit_posterior(2552, 'asia_east', 'male', '2005')
    """
    
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])

    print 'initializing model vars... ',
    dm.calc_effective_sample_size(dm.data)
    dm.vars = dismod3.gbd_disease_model.setup(dm, keys)
    print 'initialization completed'


    # fit the model
    dir = dismod3.settings.JOB_WORKING_DIR % dm.id

    ## first generate decent initial conditions
    verbose=1
    def map_fit(stoch_names):
        print '\nfitting', ' '.join(stoch_names)
        key = dismod3.utils.gbd_key_for('%s', region, year, sex)
        map = mc.MAP([dm.vars[key%type][subkey] for type in stoch_names for subkey in dm.vars[key%type] if not subkey.startswith('log_dispersion')])

        try:
            map.fit(method='fmin_powell', verbose=verbose)
        except KeyboardInterrupt:
            print 'User halted optimization routine before optimal value found'

        for type in stoch_names:
            for subkey in ['age_coeffs_mesh', 'dispersion']:
                if subkey in dm.vars[key%type]:
                    print key%type, subkey, pl.atleast_1d(dm.vars[key%type][subkey].value).round(0)

        return map

    print 'initializing MAP object... ',
    map_fit(['incidence'])
    map_fit(['remission'])
    map_fit('excess-mortality mortality relative-risk smr duration'.split())
    map_fit('remission excess-mortality mortality relative-risk smr duration'.split())
    map_fit(['incidence', 'bins', 'prevalence'])
    map_fit('incidence excess-mortality mortality relative-risk smr duration bins prevalence'.split())
    dm.map = map_fit('incidence remission excess-mortality mortality relative-risk smr duration bins prevalence'.split())
    print 'initialization completed'

    ## then sample the posterior via MCMC
    mc.warnings.warn = sys.stdout.write    # make pymc warnings go to stdout

    dm.mcmc = mc.MCMC(dm.vars)
    age_stochs = []
    for k in keys:
        #for subkey in 'log_dispersion age_coeffs_mesh'.split():
        #    if subkey in dm.vars[k] and isinstance(dm.vars[k][subkey], mc.Stochastic):
        #        dm.mcmc.use_step_method(mc.NoStepper, dm.vars[k][subkey])
        if 'dispersion_step_sd' in dm.vars[k]:
            dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                    proposal_sd=dm.vars[k]['dispersion_step_sd'])
        if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
            age_stochs.append(dm.vars[k]['age_coeffs_mesh'])
            dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                    cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)
    key = dismod3.utils.gbd_key_for('%s', region, year, sex)
    #dm.mcmc.use_step_method(mc.AdaptiveMetropolis, [dm.vars[key%type]['age_coeffs_mesh'] for type in 'incidence remission'.split()])
    #dm.mcmc.use_step_method(mc.AdaptiveMetropolis, [dm.vars[key%type]['age_coeffs_mesh'] for type in 'remission excess-mortality'.split()])
    #dm.mcmc.use_step_method(mc.AdaptiveMetropolis, [dm.vars[key%type]['age_coeffs_mesh'] for type in 'excess-mortality incidence'.split()])
    try:
        #dm.mcmc.sample(101, verbose=verbose)
        dm.mcmc.sample(iter=20000, burn=10000, thin=10, verbose=verbose)
    except KeyboardInterrupt:
        # if user cancels with cntl-c, save current values for "warm-start"
        pass

    # save results
    for k in keys:
        t,r,y,s = dismod3.utils.type_region_year_sex_from_key(k)

        if t in ['incidence', 'prevalence', 'remission', 'excess-mortality', 'mortality', 'prevalence_x_excess-mortality']:
            dismod3.neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])

        elif t in ['relative-risk', 'duration', 'incidence_x_duration']:
            dismod3.normal_model.store_mcmc_fit(dm, k, dm.vars[k])


    # generate plots of results
    for k in keys:
        t,r,y,s = dismod3.utils.type_region_year_sex_from_key(k)
        if t in ['incidence', 'prevalence', 'remission', 'excess-mortality']:
            for s in 'dispersion age_coeffs_mesh study_coeffs region_coeffs'.split():
                if s in dm.vars[k] and isinstance(dm.vars[k][s], mc.Node):
                    try:
                        Matplot.plot(dm.vars[k][s], path='%s/image/%s/%s/'%(dir, t, s))
                    except Exception, e:
                        print e

    dismod3.plotting.tile_plot_disease_model(dm, keys, defaults={})
    dm.savefig('dm-%d-posterior-%s.png' % (dm.id, dismod3.utils.gbd_key_for('all', region, year, sex)))

    # summarize fit quality graphically, as well as parameter posteriors
    for k in dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex]):
        if dm.vars[k].get('data'):
            dismod3.plotting.plot_posterior_predicted_checks(dm, k)
            dm.savefig('dm-%d-check-%s.png' % (dm.id, k))

    # save results (do this last, because it removes things from the disease model that plotting function, etc, might need
    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dm.save('dm-%d-posterior-%s-%s-%s.json' % (dm.id, region, sex, year), keys_to_save=keys)



from dismod3.neg_binom_model import countries_for

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
                value_list = pl.zeros((dismod3.settings.MAX_AGE, sample_size))

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
                value_list = pl.zeros((dismod3.settings.MAX_AGE, sample_size))

                # calculate value list for ages
                for i, gamma in enumerate(gamma_trace):
                    value_trace = neg_binom_model.predict_country_rate(key, iso3, alpha, beta, gamma,
                                                           covariates_dict, derived_covariate,
                                                           model_vars['bounds_func'],
                                                           range(dismod3.settings.MAX_AGE))

                    value_list[:, i] = value_trace
            if rate_type == 'prevalence':
                print key, iso3, neg_binom_model.country_covariates(key, iso3, covariates_dict, derived_covariate)[1], pl.sort(value_list, axis=1)[5, .5*sample_size]
                                
            # write a row
            for age in range(dismod3.settings.MAX_AGE):
                csv_f.writerow([iso3, rate_type, str(age)] + list(pl.sort(value_list, axis=1)[age, [.5*sample_size, .025*sample_size, .975*sample_size]]))

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


    dm = dismod3.load_disease_model(id)
    fit_posterior(dm, options.region, options.sex, options.year)
    return dm

if __name__ == '__main__':
    dm = main()

