import pylab as pl
import pymc as mc

import dismod3

import generic_disease_model as submodel
import neg_binom_model

def map_fit(dm, keys, map_method, stoch_names, verbose=1):
    """ Use MAP estimation from PyMC to get point estimates for the
    stochastics in dm.vars with keys that contain a substring matching
    something in the stoch_names list.

    Notes
    -----
    skip stochs that contain the string 'dispersion', since that was
    giving me trouble
    """
    vars = []
    for k in keys:
        for s in stoch_names:
            if k.find(s) != -1 and k.find('dispersion') == -1:
                vars.append(dm.vars[k])
    mc.MAP(vars).fit(method=map_method, iterlim=500, tol=.01, verbose=verbose)


def fit(dm, method='map', keys=None, iter=50000, burn=25000, thin=1, verbose=1,
        dbname='model.pickle'):
    """ Generate an estimate of the generic disease model parameters
    using maximum a posteriori liklihood (MAP) or Markov-chain Monte
    Carlo (MCMC)

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    method : string, optional
      the parameter estimation method, either 'map' or 'mcmc'

    keys : list, optional
      a list of gbd keys for the parameters to fit; it can speed up
      computation to holding some parameters constant while allowing
      others to vary

    iter : int, optional
    burn : int, optional
    thin : int, optional
      parameters for the MCMC, which control how long it takes, and
      how accurate it is
    """
    if not keys:
        keys = dismod3.utils.gbd_keys()

    if not hasattr(dm, 'vars'):
        print 'initializing model vars... ',
        dm.calc_effective_sample_size(dm.data)
        dm.vars = setup(dm, keys)
        print 'finished'

    if method == 'map':
        print 'initializing MAP object... ',
        map_method = 'fmin_powell'
        #map_method = 'fmin_l_bfgs_b'
        # consider removing dispersion stochs from these fits
        # TODO: refactor the submodel map fit into a separate method
        map_fit(dm, keys, map_method, ['incidence'])
        map_fit(dm, keys, map_method, ['remission'])
        map_fit(dm, keys, map_method, ['incidence', 'bins'])
        map_fit(dm, keys, map_method, ['excess-mortality', 'mortality', 'relative-risk', 'bins'])
        map_fit(dm, keys, map_method, ['incidence', 'bins', 'prevalence'])
        map_fit(dm, keys, map_method, ['excess-mortality', 'mortality', 'relative-risk', 'bins', 'prevalence'])
        dm.map = mc.MAP(dm.vars)
        print 'finished'

        try:
            #dm.map.fit(method=map_method, iterlim=500, tol=.001, verbose=verbose)
            map_fit(dm, keys, map_method, ['incidence', 'remission', 'excess-mortality', 'mortality', 'relative-risk', 'bins', 'prevalence'])
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        
        for k in keys:
            try:
                val = dm.vars[k]['rate_stoch'].value
                dm.set_map(k, val)
            except KeyError:
                pass

    if method == 'norm_approx':
        dm.na = mc.NormApprox(dm.vars, eps=.0001)

        try:
            dm.na.fit(method='fmin_powell', iterlim=500, tol=.00001, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

        for k in keys:
            if dm.vars[k].has_key('rate_stoch'):
                dm.set_map(k, dm.vars[k]['rate_stoch'].value)

        try:
            dm.na.sample(1000, verbose=verbose)
            for k in keys:
                # TODO: rename 'rate_stoch' to something more appropriate
                if dm.vars[k].has_key('rate_stoch'):
                    neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass

                        
    elif method == 'mcmc':
        # make pymc warnings go to stdout
        import sys
        mc.warnings.warn = sys.stdout.write
        
        dm.mcmc = mc.MCMC(dm.vars, db='pickle', dbname=dbname)
        for k in keys:
            if 'dispersion_step_sd' in dm.vars[k]:
                dm.mcmc.use_step_method(mc.Metropolis, dm.vars[k]['log_dispersion'],
                                        proposal_sd=dm.vars[k]['dispersion_step_sd'])
            if 'age_coeffs_mesh_step_cov' in dm.vars[k]:
                dm.mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars[k]['age_coeffs_mesh'],
                                        cov=dm.vars[k]['age_coeffs_mesh_step_cov'], verbose=0)

        try:
            dm.mcmc.sample(iter=iter, thin=thin, burn=burn, verbose=verbose)
        except KeyboardInterrupt:
            # if user cancels with cntl-c, save current values for "warm-start"
            pass
        dm.mcmc.db.commit()

        for k in keys:
            t,r,y,s = dismod3.utils.type_region_year_sex_from_key(k)
            
            if t in ['incidence', 'prevalence', 'remission', 'excess-mortality', 'mortality', 'prevalence_x_excess-mortality']:
                import neg_binom_model
                neg_binom_model.store_mcmc_fit(dm, k, dm.vars[k])
            elif t in ['relative-risk', 'duration', 'incidence_x_duration']:
                import normal_model
                normal_model.store_mcmc_fit(dm, k, dm.vars[k])


def setup(dm, keys):
    """ Generate the PyMC variables for a multi-region/year/sex generic
    disease model.

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
    
    Results
    -------
    vars : dict of PyMC stochs
      returns a dictionary of all the relevant PyMC objects for the
      multi-region/year/sex generic disease model.
    """
    
    vars = {}

    # for each region-year-sex triple among the keys
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for('%s', r, y, s)
                if not key%'prevalence' in keys:
                    continue

                dm.set_units(key%'prevalence', '(per person)')
                dm.set_units(key%'duration', '(years)')
                for t in 'incidence', 'remission', 'excess-mortality':
                    dm.set_units(key%t, '(per person-year)')
                    #dm.get_initial_estimate(key%t, [d for d in dm.data if relevant_to(d, t, r, y, s)])

                data = [d for d in dm.data if relevant_to(d, 'all', r, y, s)]
                # Also include data that is not sex specific (i.e. sex == 'total') here
                data += [d for d in dm.data if relevant_to(d, 'all', r, y, 'total')]  # FIXME: this will double include data with sex == 'all'
                
                sub_vars = submodel.setup(dm, key, data)
                vars.update(sub_vars)
    
    return vars

from neg_binom_model import countries_for

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
                value_list = pl.zeros((dismod3.MAX_AGE, sample_size))

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
                value_list = pl.zeros((dismod3.MAX_AGE, sample_size))

                # calculate value list for ages
                for i, gamma in enumerate(gamma_trace):
                    value_trace = neg_binom_model.predict_country_rate(key, iso3, alpha, beta, gamma,
                                                           covariates_dict, derived_covariate,
                                                           model_vars['bounds_func'],
                                                           range(101))

                    value_list[:, i] = value_trace
            if rate_type == 'prevalence':
                print key, iso3, neg_binom_model.country_covariates(key, iso3, covariates_dict, derived_covariate)[1], pl.sort(value_list, axis=1)[5, .5*sample_size]
                                
            # write a row
            for age in range(dismod3.MAX_AGE):
                csv_f.writerow([iso3, rate_type, str(age)] + list(pl.sort(value_list, axis=1)[age, [.5*sample_size, .025*sample_size, .975*sample_size]]))

    # close the file
    f_file.close()


def zip_country_level_posterior_files(id):
    """  Zip country level posterior files in the directory of the
    job_working_directory/posterior/country_level_posterior_dm-'id', 
    and then remove the directory containing the files

    Parameters
    ----------
    id : int
      The model id number
    """
    # job working directory
    job_wd = dismod3.settings.JOB_WORKING_DIR % id

    # directory containing the csv files
    directory = 'country_level_posterior_dm-' + str(id)

    try:
        import os
        import zipfile
        # open an archive for writing
        a = zipfile.ZipFile(directory + '.zip', 'w', zipfile.ZIP_DEFLATED)

        # put files into the archive
        for f in os.listdir(directory):
            print "archiving file %s" % f
            a.write(directory + '/' + f)

        # close the archive
        a.close()

    except Exception,e:
        print e
