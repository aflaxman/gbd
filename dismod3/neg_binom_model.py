""" Negative Binomial Model for a generic epidemological parameter

The Negative Binomial Model represents data for age-specific rates
according to the following formula::

    Y_i ~ NegativeBinomial(\mu_i N_i, \delta_i) / N_i
    \mu_i ~ \sum _{a = a_{i0}} ^{a_{i1}} w(a) \mu_{i, a}
    \log \mu_{i, a} = \alpha_{r_i} + \alpha_{s_i} + \alpha_{y_i} + \gamma_a + \beta^T X_i

Here Y_i, \N_i, a_{i0}, a_{i1}, r_i, s_i, y_i, and X_i are the
value, effective sample size, age range, region, sex, year,
and study-level covariates corresponding to a single age-range value
from a single study.  \alpha, \beta, \gamma, and \delta are parameters
(fixed effects and over-dispersion) that will be estimated from the
data.
"""
import sys

import pylab as pl
import pymc as mc

import dismod3
from dismod3.utils import debug, clean

def fit_emp_prior(dm, param_type, iter=100000, thin=50, burn=50000, dbname='/dev/null'):
    """ Generate an empirical prior distribution for a single disease parameter

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      The object containing all the data, (hyper)-priors, and additional
      information (like input and output age-mesh).

    param_type : str, one of 'incidence', 'prevalence', 'remission', 'excess-mortality'
      The disease parameter to work with

    Notes
    -----
    The results of this fit are stored in the disease model's params
    hash for use when fitting multiple paramter types together

    Example
    -------
    $ python2.5 gbd_fit.py 231 -t incidence
    """
    data = [d for d in dm.data if \
                d['data_type'] == '%s data' % param_type \
                and d.get('ignore') != -1]

    dm.clear_empirical_prior()

    dm.calc_effective_sample_size(data)
    dm.fit_initial_estimate(param_type, data)
    dm.vars = setup(dm, param_type, data)
    # don't do anything if there is no data for this parameter type
    if not dm.vars['data']:
        return

    debug('i: %s' % ', '.join(['%.2f' % x for x in dm.vars['rate_stoch'].value[::10]]))
    sys.stdout.flush()
    
    # fit the model
    def map_fit(stoch_names):
        print '\nfitting', ' '.join(stoch_names)
        map = mc.MAP([dm.vars[key] for key in stoch_names] + [dm.vars['observed_counts'], dm.vars['rate_potential'], dm.vars['priors']])
        try:
            map.fit(method='fmin_powell', verbose=verbose)
        except KeyboardInterrupt:
            debug('User halted optimization routine before optimal value found')
        for key in stoch_names:
            print key, dm.vars[key].value.round(2)
        sys.stdout.flush()

    def mcmc_fit(stoch_names):
        print '\nfitting', ' '.join(stoch_names)
        mcmc = mc.MCMC([dm.vars[key] for key in stoch_names] + [dm.vars['observed_counts'], dm.vars['rate_potential'], dm.vars['priors']])
        mcmc.use_step_method(mc.Metropolis, dm.vars['log_dispersion'],
                             proposal_sd=dm.vars['dispersion_step_sd'])
        # TODO: make a wrapper function for handling this adaptive metropolis setup
        stoch_list = [dm.vars['study_coeffs'], dm.vars['region_coeffs'], dm.vars['age_coeffs_mesh']]
        d1 = len(dm.vars['study_coeffs'].value)
        d2 = len(dm.vars['region_coeffs_step_cov'])
        d3 = len(dm.vars['age_coeffs_mesh_step_cov'])
        C = pl.eye(d1+d2+d3)
        C[d1:(d1+d2), d1:(d1+d2)] = dm.vars['region_coeffs_step_cov']
        C[(d1+d2):(d1+d2+d3), (d1+d2):(d1+d2+d3)] = dm.vars['age_coeffs_mesh_step_cov']
        C *= .01
        mcmc.use_step_method(mc.AdaptiveMetropolis, stoch_list, cov=C)

        # more step methods
        mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars['study_coeffs'])
        mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars['region_coeffs'], cov=dm.vars['region_coeffs_step_cov'])
        mcmc.use_step_method(mc.AdaptiveMetropolis, dm.vars['age_coeffs_mesh'], cov=dm.vars['age_coeffs_mesh_step_cov'])

        try:
            mcmc.sample(iter=10000, burn=5000, thin=5, verbose=verbose)
        except KeyboardInterrupt:
            debug('User halted optimization routine before optimal value found')
        sys.stdout.flush()

        # reset stoch values to sample mean
        for key in stoch_names:
            mean = dm.vars[key].stats()['mean']
            if isinstance(dm.vars[key], mc.Stochastic):
                dm.vars[key].value = mean
            print key, mean.round(2)

    verbose = 1
    stoch_names = 'region_coeffs age_coeffs_mesh study_coeffs'.split()
    ## start by optimizing parameters separately
    for key in stoch_names:
        map_fit([key])
    ## then fit them all together
    map_fit(stoch_names)
    # now find the over-dispersion parameter that matches these values
    map_fit(['log_dispersion'])

    # make pymc warnings go to stdout
    mc.warnings.warn = sys.stdout.write
    mcmc_fit(['log_dispersion', 'dispersion', 'study_coeffs', 'region_coeffs',
              'age_coeffs_mesh', 'age_coeffs',
              'predicted_rates', 'expected_rates', 'rate_stoch'])

    alpha = dm.vars['region_coeffs'].stats()['mean']
    beta = dm.vars['study_coeffs'].stats()['mean']
    gamma_mesh = dm.vars['age_coeffs_mesh'].stats()['mean']

    debug('a: %s' % ', '.join(['%.2f' % x for x in alpha]))
    debug('b: %s' % ', '.join(['%.2f' % x for x in pl.atleast_1d(beta)]))
    debug('g: %s' % ', '.join(['%.2f' % x for x in gamma_mesh]))
    debug('d: %.2f' % dm.vars['dispersion'].stats()['mean'])

    covariates_dict = dm.get_covariates()
    derived_covariate = dm.get_derived_covariate_values()
    X = covariates(data[0], covariates_dict)
    debug('p: %s' % ', '.join(['%.2f' % x for x in predict_rate(X, alpha, beta, gamma_mesh, dm.vars['bounds_func'], dm.get_param_age_mesh())]))

    # save the results in the param_hash
    prior_vals = dict(
        alpha=list(dm.vars['region_coeffs'].stats()['mean']),
        beta=list(pl.atleast_1d(dm.vars['study_coeffs'].stats()['mean'])),
        gamma=list(dm.vars['age_coeffs'].stats()['mean']),
        delta=float(dm.vars['dispersion'].stats()['mean']))

    prior_vals.update(
        sigma_alpha=list(dm.vars['region_coeffs'].stats()['standard deviation']),
        sigma_beta=list(pl.atleast_1d(dm.vars['study_coeffs'].stats()['standard deviation'])),
        sigma_gamma=list(dm.vars['age_coeffs'].stats()['standard deviation']),
        sigma_delta=float(dm.vars['dispersion'].stats()['standard deviation']))
    dm.set_empirical_prior(param_type, prior_vals)


    dispersion = prior_vals['delta']
    median_sample_size = pl.median([values_from(dm, d)[3] for d in dm.vars['data']] + [1000])
    debug('median effective sample size: %.1f' % median_sample_size)

    param_mesh = dm.get_param_age_mesh()
    age_mesh = dm.get_estimate_age_mesh()

    trace = zip(dm.vars['region_coeffs'].trace(), dm.vars['study_coeffs'].trace(), dm.vars['age_coeffs'].trace())[::5]
    
    for r in dismod3.gbd_regions:
        debug('predicting rates for %s' % r)
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.utils.gbd_key_for(param_type, r, y, s)
                rate_trace = []
                for a, b, g in trace:
                    rate_trace.append(predict_region_rate(key,
                                                          alpha=a,
                                                          beta=b,
                                                          gamma=g,
                                                          covariates_dict=covariates_dict,
                                                          derived_covariate=derived_covariate,
                                                          bounds_func=dm.vars['bounds_func'],
                                                          ages=dm.get_estimate_age_mesh()))
                mu = dismod3.utils.interpolate(param_mesh, pl.mean(rate_trace, axis=0)[param_mesh], age_mesh)
                dm.set_initial_value(key, mu)
                dm.set_mcmc('emp_prior_mean', key, mu)

                # similar to saving upper_ui and lower_ui in function store_mcmc_fit below
                rate_trace = pl.sort(rate_trace, axis=0)
                dm.set_mcmc('emp_prior_upper_ui', key, dismod3.utils.interpolate(param_mesh, rate_trace[.975 * len(rate_trace), :][param_mesh], age_mesh))
                dm.set_mcmc('emp_prior_lower_ui', key, dismod3.utils.interpolate(param_mesh, rate_trace[.025 * len(rate_trace), :][param_mesh], age_mesh))

def calc_rate_trace(dm, key, model_vars):

    covariates_dict = dm.get_covariates()
    derived_covariate = dm.get_derived_covariate_values()

    rate_trace = []

    if isinstance(model_vars['region_coeffs'], mc.Stochastic) and isinstance(model_vars['study_coeffs'], mc.Stochastic):
        for alpha, beta, gamma in zip(model_vars['region_coeffs'].trace(), model_vars['study_coeffs'].trace(), model_vars['age_coeffs'].trace()):
            mu = predict_region_rate(key, alpha, beta, gamma, covariates_dict, derived_covariate,
                                     model_vars['bounds_func'], dm.get_estimate_age_mesh())
            rate_trace.append(mu)
    else:
        alpha = model_vars['region_coeffs']
        beta = model_vars['study_coeffs']
        for gamma in model_vars['age_coeffs'].trace():
            mu = predict_region_rate(key, alpha, beta, gamma, covariates_dict, derived_covariate,
                                     model_vars['bounds_func'], dm.get_estimate_age_mesh())
            rate_trace.append(mu)

    return pl.array(rate_trace)

def store_mcmc_fit(dm, key, model_vars=None, rate_trace=None):
    """ Store the parameter estimates generated by an MCMC fit of the
    negative-binomial model in the disease_model object, keyed by key
    
    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)

    key : str

    model_vars : dict of PyMC stochastic or deterministic variable

    Results
    -------
    Save a regional estimate of the model prediction, with uncertainty
    """
    if rate_trace == None:
        rate_trace = calc_rate_trace(dm, key, model_vars)

    rate_trace = pl.sort(rate_trace, axis=0)
    rate = {}
    for x in [2.5, 50, 97.5]:
        rate[x] = rate_trace[x/100.*len(rate_trace), :]
    param_mesh = dm.get_param_age_mesh()
    age_mesh = dm.get_estimate_age_mesh()
    
    dm.set_mcmc('lower_ui', key, rate[2.5])
    dm.set_mcmc('median', key, rate[50])
    dm.set_mcmc('upper_ui', key, rate[97.5])
    dm.set_mcmc('mean', key, pl.mean(rate_trace,axis=0))

    if dm.vars[key].has_key('dispersion'):
        dm.set_mcmc('dispersion', key, dm.vars[key]['dispersion'].stats()['quantiles'].values())

def covariate_names(dm):
    covariate_list = []
    covariates_dict = dm.get_covariates()
    
    for level in ['Study_level', 'Country_level']:
        for k in sorted(covariates_dict[level]):
            if covariates_dict[level][k]['rate']['value'] == 1:
                covariate_list.append(k)
    if covariate_list == []:
        covariate_list.append('(none)')

    return covariate_list
                                
def covariates(d, covariates_dict):
    """ extract the covariates from a data point as a vector;

    Xa represents region-level covariates:
      Xa[0],...,Xa[21] = region indicators
      Xa[22] = .1*(year-1997)
      Xa[23] = .5 if sex == 'male', -.5 if sex == 'female'
    Xb represents study-level covariates, according to the covariates_dict
      
    """
    Xa = pl.zeros(len(dismod3.gbd_regions) + 2)
    for ii, r in enumerate(dismod3.gbd_regions):
        if clean(d['gbd_region']) == clean(r):
            Xa[ii] = 1.

    if d['year_start'] == 'all':
        Xa[ii+1] = 0.
    else:
        Xa[ii+1] = .1 * (.5 * (float(d['year_start']) + float(d['year_end'])) - 1997)
    
    if clean(d['sex']) == 'male':
        Xa[ii+2] = .5
    elif clean(d['sex']) == 'female':
        Xa[ii+2] = -.5
    else:
        Xa[ii+2] = 0.

    Xb = []
    for level in ['Study_level', 'Country_level']:
        for k in sorted(covariates_dict[level]):
            if covariates_dict[level][k]['rate']['value'] == 1:
                Xb.append(float(d.get(clean(k)) or 0.))
    #debug('%s-%s-%s-%s: Xb = %s' % (d['sex'], d['year_start'], d['gbd_region'], d.get('country_iso3_code', 'none'), str(Xb)))
    if Xb == []:
        Xb = [0.]
    return Xa, Xb


from dismod3.utils import clean
import csv
import settings
countries_for = dict(
    [[clean(x[0]), x[1:]] for x in csv.reader(open(settings.CSV_PATH + 'country_region.csv'))]
    )
population_by_age = dict(
    [[(d['Country Code'], d['Year'], d['Sex']),
      [max(.001,float(d['Age %d Population' % i])) for i in range(dismod3.settings.MAX_AGE)]] for d in csv.DictReader(open(settings.CSV_PATH + 'population.csv'))
     if len(d['Country Code']) == 3]
    )

def regional_population(key):
    """ calculate regional population for a gbd key"""
    t,r,y,s = dismod3.utils.type_region_year_sex_from_key(key)
    pop = pl.zeros(dismod3.settings.MAX_AGE)
    for c in countries_for[clean(r)]:
        if y == 'all' and s == 'all':
            for yy in dismod3.settings.gbd_years:
                for ss in dismod3.settings.gbd_sexes:
                    pop += population_by_age[(c, yy, dismod3.utils.clean(ss))]
        else:
            pop += population_by_age[(c, y, s)]
    return pop

def regional_average(derived_covariate, key, region, year, sex):
    """ handle region = iso3 code or region = clean(gbd_region)"""
    # TODO: make regional average weighted by population
    if key not in derived_covariate:
        debug('WARNING: derived covariate %s not found' % key)
        return 0.

    if region == 'world':
        return 0.

    cov_vals = [derived_covariate[key]['%s+%s+%s'%(iso3,year,sex)] for iso3 in countries_for[region]
                if derived_covariate[key].has_key('%s+%s+%s'%(iso3,year,sex))]
    return pl.mean(cov_vals)

# store computed covariate data for fast access later
# TODO: ensure that this hash table is cleared between runs of different models! otherwise it can break things.
covariate_hash = {}

def regional_covariates(key, covariates_dict, derived_covariate):
    """ form the covariates for a gbd key"""
    if not key in covariate_hash:
        try:
            t,r,y,s = dismod3.utils.type_region_year_sex_from_key(key)
        except KeyError:
            r = 'world'
            y = 1997
            s = 'total'

        d = {'gbd_region': r,
             'year_start': y,
             'year_end': y,
             'sex': s}
        for level in ['Study_level', 'Country_level']:
            for k in covariates_dict[level]:
                if k == 'none':
                    continue
                if covariates_dict[level][k]['rate']['value']:
                    d[clean(k)] = covariates_dict[level][k]['value']['value']
                    if level == 'Country_level':
                        d[clean(k)] = regional_average(derived_covariate, k, r, y, s)
                    else:
                        d[clean(k)] = float(d[clean(k)] or 0.)

        covariate_hash[key] = covariates(d, covariates_dict)
    
    return covariate_hash[key]

def country_covariates(key, iso3, covariates_dict, derived_covariate):
    """ form the covariates for a gbd key"""
    if not (key, iso3) in covariate_hash:
        t,r,y,s = dismod3.utils.type_region_year_sex_from_key(key)

        d = {'gbd_region': r,
             'year_start': y,
             'year_end': y,
             'sex': s}
        for level in ['Study_level', 'Country_level']:
            for k in covariates_dict[level]:
                if k == 'none':
                    continue
                if covariates_dict[level][k]['rate']['value']:
                    d[clean(k)] = covariates_dict[level][k]['value']['value']
                    if level == 'Country_level':
                        if k not in derived_covariate:
                            debug('WARNING: derived covariate %s not found' % key)
                            d[clean(k)] = 0.
                        elif not derived_covariate[k].has_key('%s+%s+%s'%(iso3,y,s)):
                            debug('WARNING: derived covariate %s not found for (%s, %s, %s)' % (k, iso3, y, s))
                            d[clean(k)] = 0.
                        else:
                            d[clean(k)] = derived_covariate[k].get('%s+%s+%s'%(iso3,y,s), 0.)
                    else:
                        d[clean(k)] = float(d[clean(k)] or 0.)

        covariate_hash[(key, iso3)] = covariates(d, covariates_dict)
    return covariate_hash[(key, iso3)]

def predict_rate(X, alpha, beta, gamma, bounds_func, ages):
    Xa, Xb = X
    return bounds_func(pl.exp(pl.dot(Xa, alpha) + pl.dot(Xb, beta) + gamma), ages)

def predict_country_rate(key, iso3, alpha, beta, gamma, covariates_dict, derived_covariate, bounds_func, ages):
    return predict_rate(country_covariates(key, iso3, covariates_dict, derived_covariate), alpha, beta, gamma, bounds_func, ages)

def predict_region_rate(key, alpha, beta, gamma, covariates_dict, derived_covariate, bounds_func, ages):
    t,r,y,s = dismod3.utils.type_region_year_sex_from_key(key)
    region_rate = pl.zeros(len(gamma))
    total_pop = pl.zeros(len(gamma))
    for iso3 in countries_for[r]:
        region_rate += predict_country_rate(key, iso3, alpha, beta, gamma, covariates_dict, derived_covariate, bounds_func, ages) * population_by_age.get((iso3,y,s), 1.)
        total_pop += population_by_age.get((iso3, y, s), 1.)
    return region_rate / total_pop
    
def setup(dm, key, data_list=[], rate_stoch=None, emp_prior={}, lower_bound_data=[]):
    """ Generate the PyMC variables for a negative-binomial model of
    a single rate function

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      the object containing all the data, priors, and additional
      information (like input and output age-mesh)
      
    key : str
      the name of the key for everything about this model (priors,
      initial values, estimations)

    data_list : list of data dicts
      the observed data to use in the negative binomial liklihood function

    rate_stoch : pymc.Stochastic, optional
      a PyMC stochastic (or deterministic) object, with
      len(rate_stoch.value) == len(dm.get_estimation_age_mesh()).
      This is used to link rate stochs into a larger model,
      for example.

    emp_prior : dict, optional
      the empirical prior dictionary, retrieved from the disease model
      if appropriate by::

          >>> t, r, y, s = dismod3.utils.type_region_year_sex_from_key(key)
          >>> emp_prior = dm.get_empirical_prior(t)

    Results
    -------
    vars : dict
      Return a dictionary of all the relevant PyMC objects for the
      rate model.  vars['rate_stoch'] is of particular
      relevance; this is what is used to link the rate model
      into more complicated models, like the generic disease model.
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    param_mesh = dm.get_param_age_mesh()

    if pl.any(pl.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    # calculate effective sample size for all data and lower bound data
    dm.calc_effective_sample_size(data_list)
    dm.calc_effective_sample_size(lower_bound_data)

    # generate regional covariates
    covariate_dict = dm.get_covariates()
    derived_covariate = dm.get_derived_covariate_values()
    X_region, X_study = regional_covariates(key, covariate_dict, derived_covariate)

    # use confidence prior from prior_str
    mu_delta = 100.
    sigma_delta = 10.
    from dismod3.settings import PRIOR_SEP_STR
    for line in dm.get_priors(key).split(PRIOR_SEP_STR):
        prior = line.strip().split()
        if len(prior) == 0:
            continue
        if prior[0] == 'heterogeneity':
            mu_delta = float(prior[1])
            sigma_delta = float(prior[2])

    # use the empirical prior mean if it is available
    if len(set(emp_prior.keys()) & set(['alpha', 'beta', 'gamma'])) == 3:
        mu_alpha = pl.array(emp_prior['alpha'])
        sigma_alpha = pl.array(emp_prior['sigma_alpha'])
        alpha = pl.array(emp_prior['alpha']) # TODO: make this stochastic
        vars.update(region_coeffs=alpha)

        beta = pl.array(emp_prior['beta']) # TODO: make this stochastic
        sigma_beta = pl.array(emp_prior['sigma_beta'])
        vars.update(study_coeffs=beta)

        mu_gamma = pl.array(emp_prior['gamma'])
        sigma_gamma = pl.array(emp_prior['sigma_gamma'])

        if 'delta' in emp_prior:
            mu_delta = emp_prior['delta']
            if 'sigma_delta' in emp_prior:
                sigma_delta = emp_prior['sigma_delta']
    else:
        import dismod3.regional_similarity_matrices as similarity_matrices
        n = len(X_region)
        mu_alpha = pl.zeros(n)
        sigma_alpha = .025  # TODO: make this a hyperparameter, with a traditional prior, like inverse gamma
        C_alpha = similarity_matrices.regions_nested_in_superregions(n, sigma_alpha)

        # use alternative region effect covariance structure if requested
        region_prior_key = 'region_effects'
        if region_prior_key in dm.params:
            if dm.params[region_prior_key] == 'uninformative':
                C_alpha = similarity_matrices.uninformative(n, sigma_alpha)

        region_prior_key = 'region_effect_%s'%key.split(dismod3.settings.KEY_DELIM_CHAR)[0]
        if region_prior_key in dm.params:
            if dm.params[region_prior_key] == 'uninformative':
                C_alpha = similarity_matrices.regions_nested_in_superregions(n, dm.params[region_prior_key]['std'])

        # add informative prior for sex effect if requested
        sex_prior_key = 'sex_effect_%s'%key.split(dismod3.settings.KEY_DELIM_CHAR)[0]
        if sex_prior_key in dm.params:
            print 'adjusting prior on sex effect coefficient for %s' % key
            mu_alpha[n-1] = pl.log(dm.params[sex_prior_key]['mean'])
            sigma_sex = (pl.log(dm.params[sex_prior_key]['upper_ci']) - pl.log(dm.params[sex_prior_key]['lower_ci'])) / (2*1.96)
            C_alpha[n-1, n-1]= sigma_sex**2.

        # add informative prior for time effect if requested
        time_prior_key = 'time_effect_%s'%key.split(dismod3.settings.KEY_DELIM_CHAR)[0]  # HACK: sometimes key is just parameter type, sometimes it is type+region+year+sex
        if time_prior_key in dm.params:
            print 'adjusting prior on time effect coefficient for %s' % key
            mu_alpha[n-2] = pl.log(dm.params[time_prior_key]['mean'])
            sigma_time = (pl.log(dm.params[time_prior_key]['upper_ci']) - pl.log(dm.params[time_prior_key]['lower_ci'])) / (2*1.96)
            C_alpha[n-2, n-2]= sigma_time**2.
        
        #C_alpha = similarity_matrices.all_related_equally(n, sigma_alpha)
        alpha = mc.MvNormalCov('region_coeffs_%s' % key, mu=mu_alpha,
                            C=C_alpha,
                            value=mu_alpha)
        vars.update(region_coeffs=alpha, region_coeffs_step_cov=.005*C_alpha)

        mu_beta = pl.zeros(len(X_study))
        sigma_beta = .05
        beta = mc.Normal('study_coeffs_%s' % key, mu=mu_beta, tau=sigma_beta**-2., value=mu_beta)
        vars.update(study_coeffs=beta)

        mu_gamma = 0.*pl.ones(len(est_mesh))
        sigma_gamma = 2.*pl.ones(len(est_mesh))

    if mu_delta != 0.:
        log_delta = mc.Normal('log_dispersion_%s' % key, mu=pl.log(mu_delta)/pl.log(10), tau=.25**-2, value=pl.log(mu_delta)/pl.log(10))
        delta = mc.Lambda('dispersion_%s' % key, lambda x=log_delta: .5 + 10.**x)
        
        vars.update(dispersion=delta, log_dispersion=log_delta, dispersion_step_sd=.1*log_delta.parents['tau']**-.5)
    else:
        delta = mc.Lambda('dispersion_%s' % key, lambda mu=mu_delta: 0)
        vars.update(dispersion=delta)

    if len(sigma_gamma) == 1:
        sigma_gamma = sigma_gamma[0]*pl.ones(len(est_mesh))

    # create varible for interpolated rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        # if the rate_stoch already exists, for example prevalence in the generic model,
        # we use it to back-calculate mu and eventually gamma
        mu = rate_stoch

        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(mu=mu, Xa=X_region, Xb=X_study, alpha=alpha, beta=beta):
            return pl.log(pl.maximum(dismod3.settings.NEARLY_ZERO, mu)) - pl.dot(alpha, Xa) - pl.dot(beta, Xb)

        @mc.potential(name='age_coeffs_potential_%s' % key)
        def gamma_potential(gamma=gamma, mu_gamma=mu_gamma, tau_gamma=1./sigma_gamma[param_mesh]**2, param_mesh=param_mesh):
            return mc.normal_like(gamma[param_mesh], mu_gamma[param_mesh], tau_gamma)

        vars.update(rate_stoch=mu, age_coeffs=gamma, age_coeffs_potential=gamma_potential)
    else:
        # if the rate_stoch does not yet exists, we make gamma a stoch, and use it to calculate mu
        # for computational efficiency, gamma is a linearly interpolated version of gamma_mesh
        initial_gamma = pl.log(dismod3.settings.NEARLY_ZERO + dm.get_initial_value(key))

        gamma_mesh = mc.Normal('age_coeffs_mesh_%s' % key, mu=mu_gamma[param_mesh], tau=sigma_gamma[param_mesh]**-2, value=initial_gamma[param_mesh])

        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(gamma_mesh=gamma_mesh, param_mesh=param_mesh, est_mesh=est_mesh):
            return dismod3.utils.interpolate(param_mesh, gamma_mesh, est_mesh)

        @mc.deterministic(name=key)
        def mu(Xa=X_region, Xb=X_study, alpha=alpha, beta=beta, gamma=gamma):
            return predict_rate([Xa, Xb], alpha, beta, gamma, lambda f, age: f, est_mesh)

        # Create a guess at the covariance matrix for MCMC proposals to update gamma_mesh
        from pymc.gp.cov_funs import matern
        a = pl.atleast_2d(param_mesh).T
        C = matern.euclidean(a, a, diff_degree = 2, amp = 1.**2, scale = 10.)

        vars.update(age_coeffs_mesh=gamma_mesh, age_coeffs=gamma, rate_stoch=mu, age_coeffs_mesh_step_cov=.005*pl.array(C))

        # adjust value of gamma_mesh based on priors, if necessary
        # TODO: implement more adjustments, currently only adjusted based on at_least priors
        for line in dm.get_priors(key).split(PRIOR_SEP_STR):
            prior = line.strip().split()
            if len(prior) == 0:
                continue
            if prior[0] == 'at_least':
                delta_gamma = pl.log(pl.maximum(mu.value, float(prior[1]))) - pl.log(mu.value)
                gamma_mesh.value = gamma_mesh.value + delta_gamma[param_mesh]

    # create potentials for priors
    dismod3.utils.generate_prior_potentials(vars, dm.get_priors(key), est_mesh)

    # create observed stochastics for data
    vars['data'] = []

    if mu_delta != 0.:  
        value = []
        N = []
        Xa = []
        Xb = []
        ai = []
        aw = []

        for d in data_list:
            try:
                age_indices, age_weights, Y_i, N_i = values_from(dm, d)
            except ValueError:
                debug('WARNING: could not calculate likelihood for data %d' % d['id'])
                continue

            value.append(Y_i*N_i)
            N.append(N_i)
            Xa.append(covariates(d, covariate_dict)[0])
            Xb.append(covariates(d, covariate_dict)[1])
            ai.append(age_indices)
            aw.append(age_weights)

            vars['data'].append(d)

        N = pl.array(N)
        Xa = pl.array(Xa)
        Xb = pl.array(Xb)
        value = pl.array(value)
        
        vars['effective_sample_size'] = list(N)
        
    if len(vars['data']) > 0:
        # TODO: consider using only a subset of the rates at each step of the fit to speed computation; say 100 of them
        k = 50000
        if len(vars['data']) < k:
            data_sample = range(len(vars['data']))
        else:
            import random
            @mc.deterministic(name='data_sample_%s' % key)
            def data_sample(n=len(vars['data']), k=k):
                return random.sample(range(n), k)

        @mc.deterministic(name='rate_%s' % key)
        def rates(S=data_sample,
                Xa=Xa, Xb=Xb,
                alpha=alpha, beta=beta, gamma=gamma,
                bounds_func=vars['bounds_func'],
                age_indices=ai,
                age_weights=aw):

            # calculate study-specific rate function
            shifts = pl.exp(pl.dot(Xa[S], alpha) + pl.dot(Xb[S], pl.atleast_1d(beta)))
            exp_gamma = pl.exp(gamma)
            mu = pl.zeros_like(shifts)
            for i,s in enumerate(S):
                mu[i] = pl.dot(age_weights[s], bounds_func(shifts[i] * exp_gamma[age_indices[s]], age_indices[s]))
                # TODO: evaluate speed increase and accuracy decrease of the following:
                #midpoint = age_indices[s][len(age_indices[s])/2]
                #mu[i] = bounds_func(shifts[i] * exp_gamma[midpoint], midpoint)
                # TODO: evaluate speed increase and accuracy decrease of the following: (to see speed increase, need to code this up using difference of running sums
                #mu[i] = pl.dot(pl.ones_like(age_weights[s]) / float(len(age_weights[s])),
                #               bounds_func(shifts[i] * exp_gamma[age_indices[s]], age_indices[s]))
            return mu
        vars['expected_rates'] = rates
        
        @mc.observed
        @mc.stochastic(name='data_%s' % key)
        def obs(value=value,
                S=data_sample,
                N=N,
                mu_i=rates,
                delta=delta):
            #zeta_i = .001
            #residual = pl.log(value[S] + zeta_i) - pl.log(mu_i*N[S] + zeta_i)
            #return mc.normal_like(residual, 0, 100. + delta)
            logp = mc.negative_binomial_like(value[S], N[S]*mu_i, delta)
            return logp

        vars['observed_counts'] = obs

        @mc.deterministic(name='predicted_data_%s' % key)
        def predictions(value=value,
                        N=N,
                        S=data_sample,
                        mu=rates,
                        delta=delta):
            r_S = mc.rnegative_binomial(N[S]*mu, delta)/N[S]
            r = pl.zeros(len(vars['data']))
            r[S] = r_S
            return r

        vars['predicted_rates'] = predictions
        debug('likelihood of %s contains %d rates' % (key, len(vars['data'])))

    # now do the same thing for the lower bound data
    # TODO: refactor to remove duplicated code
    vars['lower_bound_data'] = []
    value = []
    N = []
    Xa = []
    Xb = []
    ai = []
    aw = []
    for d in lower_bound_data:
        try:
            age_indices, age_weights, Y_i, N_i = values_from(dm, d)
        except ValueError:
            debug('WARNING: could not calculate likelihood for data %d' % d['id'])
            continue

        value.append(Y_i*N_i)
        N.append(N_i)
        Xa.append(covariates(d, covariate_dict)[0])
        Xb.append(covariates(d, covariate_dict)[1])
        ai.append(age_indices)
        aw.append(age_weights)

        vars['lower_bound_data'].append(d)

    N = pl.array(N)
    value = pl.array(value)

    if len(vars['lower_bound_data']) > 0:
        @mc.observed
        @mc.stochastic(name='lower_bound_data_%s' % key)
        def obs_lb(value=value, N=N,
                   Xa=Xa, Xb=Xb,
                   alpha=alpha, beta=beta, gamma=gamma,
                   bounds_func=vars['bounds_func'],
                   delta=delta,
                   age_indices=ai,
                   age_weights=aw):

            # calculate study-specific rate function
            shifts = pl.exp(pl.dot(Xa, alpha) + pl.dot(Xb, pl.atleast_1d(beta)))
            exp_gamma = pl.exp(gamma)
            mu_i = [pl.dot(weights, bounds_func(s_i * exp_gamma[ages], ages)) for s_i, ages, weights in zip(shifts, age_indices, age_weights)]  # TODO: try vectorizing this loop to increase speed
            rate_param = mu_i*N
            violated_bounds = pl.nonzero(rate_param < value)
            logp = mc.negative_binomial_like(value[violated_bounds], rate_param[violated_bounds], delta)
            return logp

        vars['observed_lower_bounds'] = obs_lb
        debug('likelihood of %s contains %d lowerbounds' % (key, len(vars['lower_bound_data'])))

    return vars


def values_from(dm, d):
    """ Extract the normalized values from a piece of data

    Parameters
    ----------
    dm : disease model

    d : data dict
    """
    est_mesh = dm.get_estimate_age_mesh()

    # get the index vector and weight vector for the age range
    age_indices = dismod3.utils.indices_for_range(est_mesh, d['age_start'], d['age_end'])
    age_weights = d.get('age_weights', pl.ones(len(age_indices))/len(age_indices))

    # ensure all rate data is valid
    Y_i = dm.value_per_1(d)
    if Y_i < 0:
        debug('WARNING: data %d < 0' % d['id'])
        raise ValueError

    N_i = max(1., d['effective_sample_size'])
    return age_indices, age_weights, Y_i, N_i
