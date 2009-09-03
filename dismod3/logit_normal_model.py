""" Logit-Normal Model for a generic epidemological parameter

The Logit-Normal Model represents data for age-specific risks/ratios
according to the following formula::

    logit(Y_i) ~ \sum _{a = a_{i0}} ^{a_{i1}} w(a) \mu_{i, a} + N(0, \sigma_i^2 + \sigma_d^2)
    \mu_{i, a} = \alpha_r + \gamma_a + \beta^T X_i

Here Y_i, \sigma_i, a_{i0}, a_{i1}, and X_i are the value, standard
error, age range, and covariates corresponding to a single age-range
value from a single study.  \alpha, \beta, \gamma, and \sigma_d are parameters
that will be estimated from the data.

"""

import numpy as np
import pymc as mc
import random

import dismod3
from dismod3.utils import debug, interpolate, rate_for_range, indices_for_range, generate_prior_potentials, gbd_regions, clean, type_region_year_sex_from_key
from dismod3.settings import NEARLY_ZERO, MISSING

# re-use the beta_binomial_model's store_mcmc_fit function
# (might need to override this in the future)
from dismod3.beta_binomial_model import store_mcmc_fit

def fit_emp_prior(dm, param_type):
    """ Generate an empirical prior distribution for a single disease parameter

    Parameters
    ----------
    dm : dismod3.DiseaseModel
      The object containing all the data, (hyper)-priors, and additional
      information (like input and output age-mesh).

    param_type : str, one of 'incidence', 'prevalence', 'remission', 'case-fatality'
      The disease parameter to work with

    Notes
    -----
    The results of this fit are stored in the disease model's params
    hash for use when fitting multiple paramter types together

    Example
    -------
    $ python2.5 gbd_fit.py 175 -t incidence -p 'zero 0 4, zero 41 100, smooth 25' # takes 7m to run
    """

    data = [d for d in dm.data if clean(d['data_type']).find(param_type) != -1]

    # don't do anything if there is no data for this parameter type
    if len(data) == 0:
        return
    
    dm.fit_initial_estimate(param_type, data)

    dm.vars = setup(dm, param_type, data)
    
    # fit the model
    dm.map = mc.MAP(dm.vars)
    try:
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.00001, verbose=1)
    except KeyboardInterrupt:
        print 'User halted optimization routine before optimal value found'
    
    # save the results in the param_hash
    dm.clear_empirical_prior()
    alpha = dm.vars['region_coeffs'].value
    beta = dm.vars['study_coeffs'].value
    gamma = dm.vars['age_coeffs'].value
    dispersion = dm.vars['dispersion'].value
    dm.set_empirical_prior(param_type, {'alpha': list(alpha),
                                        'beta': list(beta),
                                        'gamma': list(gamma),
                                        'dispersion': float(dispersion)})

    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for(param_type, r, y, s)
                logit_mu = predict_logit_risk(regional_covariates(key), alpha, beta, gamma)
                mu = mc.invlogit(logit_mu)
                dm.set_initial_value(key, mu)
                dm.set_mcmc('emp_prior_mean', key, mu)
                dm.set_mcmc('emp_prior_lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
                dm.set_mcmc('emp_prior_upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

    key = dismod3.gbd_key_for(param_type, 'world', 1997, 'total')
    logit_mu = predict_logit_risk(regional_covariates(key), alpha, beta, gamma)
    mu = mc.invlogit(logit_mu)
    dm.set_initial_value(key, mu)
    dm.set_mcmc('emp_prior_mean', key, mu)
    dm.set_mcmc('emp_prior_lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
    dm.set_mcmc('emp_prior_upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

def covariates(d):
    """ extract the covariates from a data point as a vector;
    X[0],...,X[21] = region indicators
    X[22] = year-1997
    X[23] = 1 if sex == 'male', -1 if sex == 'female'
    """
    Xa = np.zeros(len(gbd_regions) + 2)
    for ii, r in enumerate(gbd_regions):
        if clean(d['gbd_region']) == clean(r):
            Xa[ii] = 1.

    Xa[ii+1] = .5 * (float(d['year_start']) + float(d['year_end'])) - 1997

    if clean(d['sex']) == 'male':
        Xa[ii+2] = 1.
    elif clean(d['sex']) == 'female':
        Xa[ii+2] = -1.
    else:
        Xa[ii+2] = 0.

    Xb = np.zeros(5.)

    # TODO: instead of hard-coding this, store it in the disease model
    # (and let users set it through the web)
    if d.get('self_reported'):
        Xb[0] = 1.
        
    return Xa, Xb

def regional_covariates(key):
    """ form the covariates for a gbd key"""
    t,r,y,s = type_region_year_sex_from_key(key)

    d = {'gbd_region': r,
         'year_start': y,
         'year_end': y,
         'sex': s}
    return covariates(d)

def predict_risk(X, alpha, beta, gamma):
    return mc.invlogit(predict_logit_risk(X, alpha, beta, gamma))

def predict_logit_risk(X, alpha, beta, gamma):
    """ Calculate Y = gamma + X * beta"""
    Xa, Xb = X
    return np.dot(Xa, alpha) + np.dot(Xb, beta) + gamma

def setup(dm, key, data_list, rate_stoch=None, emp_prior={}):
    """ Generate the PyMC variables for a logit-normal model of
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
      the observed data to use in the logit-normal liklihood function

    rate_stoch : pymc.Stochastic, optional
      a PyMC stochastic (or deterministic) object, with
      len(rate_stoch.value) == len(dm.get_estimation_age_mesh()).
      This is used to link rate stochs into a larger model,
      for example.

    emp_prior : dict, optional
      the empirical prior dictionary, retrieved from the disease model
      if appropriate by::

          >>> t, r, y, s = type_region_year_sex_from_key(key)
          >>> emp_prior = dm.get_empirical_prior(t)

    Results
    -------
    vars : dict
      Return a dictionary of all the relevant PyMC objects for the
      rate model.  vars['rate_stoch'] is of particular
      relevance; this is what is used to link the rate model
      into more complicated models, like the generic disease model.

    Details
    -------
    The logit-normal model parameters are the following:
      * the mean age-specific rate function
      * systematic error in this mean
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    param_mesh = dm.get_param_age_mesh()

    if np.any(np.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    # for debugging
    #if key == 'prevalence+caribbean+2005+male':
    #    import pdb; pdb.set_trace()

    # use the empirical prior mean if it is available
    Xa, Xb = regional_covariates(key)

    if set(emp_prior.keys()) == set(['alpha', 'beta', 'gamma', 'dispersion']):
        mu_alpha = np.array(emp_prior['alpha'])
        mu_beta = np.array(emp_prior['beta'])
        mu_gamma = np.array(emp_prior['gamma'])
        sigma = emp_prior['dispersion']
        mu_dispersion = sigma
    else:
        mu_alpha = np.zeros(len(Xa))
        mu_beta = np.zeros(len(Xb))
        mu_gamma = -5.*np.ones(len(est_mesh))
        sigma = 10.
        mu_dispersion = .1
        
    # try using fully bayesian dispersion parameter
    mu_dispersion = .1

    # create varible for interpolated logit rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        vars['rate_stoch'] = rate_stoch

        @mc.deterministic(name='gamma_II(%s)' % key)
        def alpha_plus_gamma(rate_stoch=rate_stoch):
            return mc.logit(rate_stoch)
        vars['apg'] = alpha_plus_gamma

        @mc.deterministic(name='gamma(%s)' % key)
        def gamma(apg=alpha_plus_gamma, gamma=mu_gamma):
            return apg - gamma
        vars['age_coeffs'] = gamma

        # FIXME: if est_mesh doesn't start at 0, x[param_mesh] != interpolate(est_mesh, logit_mu, param_mesh)
        mu = predict_logit_risk([Xa, Xb], mu_alpha, mu_beta, mu_gamma)
        @mc.potential(name='interp_logit_empirical_prior_potential_%s' % key)
        def emp_prior_potential(mu=mu, Y=alpha_plus_gamma, tau=1./sigma**2, mesh=param_mesh):
            return mc.normal_like(Y[mesh], mu[mesh], tau)
        vars['rate_emp_prior_potential'] = emp_prior_potential
        
    else:
        initial_value = dm.get_initial_value(key)
        mu_logit_rate = predict_logit_risk([Xa, Xb], mu_alpha, mu_beta, mu_gamma)

        logit_rate = mc.Normal('logit(%s)' % key,
                               mu=mu_logit_rate[param_mesh],
                               tau=1./sigma**2,
                               value=mc.logit(initial_value)[param_mesh])
        vars['logit_rate'] = logit_rate

        @mc.deterministic(name='interp_logit(%s)' % key)
        def gamma(logit_rate=logit_rate):
            return interpolate(param_mesh, logit_rate, est_mesh)
        vars['age_coeffs'] = gamma

        @mc.deterministic(name=key)
        def rate_stoch(gamma=gamma):
            return mc.invlogit(gamma)
        vars['rate_stoch'] = rate_stoch


    region_coeffs = mc.Normal('region_coeffs_%s' % key, mu=mu_alpha, tau=1/.1**2)
    vars['region_coeffs'] = region_coeffs

    study_coeffs = mc.Normal('study_coeffs_%s' % key, mu=mu_beta, tau=1/.1**2)
    vars['study_coeffs'] = study_coeffs

    log_dispersion = mc.Uninformative('log(dispersion_%s)' % key, value=np.log(mu_dispersion))

    @mc.deterministic(name='dispersion_%s' % key)
    def dispersion(log_dispersion=log_dispersion):
        return np.exp(log_dispersion)

    @mc.potential(name='dispersion_potential_%s' % key)
    def dispersion_potential(dispersion=dispersion, alpha=10., beta=10. / mu_dispersion):
        return mc.gamma_like(dispersion, alpha, beta)

    vars['log_dispersion'] = log_dispersion
    vars['dispersion'] = dispersion
    vars['dispersion_potential'] = dispersion_potential


    # create potentials for priors
    vars['priors'] = generate_prior_potentials(dm.get_priors(key), est_mesh, rate_stoch)
    
    
    # create observed stochastics for data
    vars['data'] = data_list
    vars['observed_rates'] = []

    min_val = min([1.e-9] + [dm.value_per_1(d) for d in data_list if dm.value_per_1(d) > 0]) # TODO: assess validity of this minimum value
    max_se = max([.000001] + [dm.se_per_1(d) for d in data_list if dm.se_per_1(d) > 0])  # TODO: assess validity of this maximum std err

    for d in data_list:
        try:
            age_indices, age_weights, logit_val, logit_se = values_from(dm, d, min_val, max_se)
        except ValueError:
            continue

        @mc.observed
        @mc.stochastic(name='data_%d' % d['id'])
        def obs(value=logit_val,
                gamma=gamma,
                logit_se=logit_se,
                dispersion=dispersion,
                age_indices=age_indices,
                age_weights=age_weights,
                alpha=region_coeffs,
                beta=study_coeffs,
                X=covariates(d)):
            mean_val = rate_for_range(predict_logit_risk(X, alpha, beta, gamma), age_indices, age_weights)
            logit_disp = (1/mean_val + 1/(1-mean_val)) * dispersion
            return mc.normal_like(x=value, mu=mean_val, tau=1. / (logit_disp**2 + logit_se**2))
            
        vars['observed_rates'].append(obs)
        
    return vars


def values_from(dm, d, min_val=1.e-5, max_se=.1):
    """ Extract the normalized values from a piece of data

    Parameters
    ----------
    dm : disease model

    d : data dict

    min_val : float, optional
      the value to use instead of zero, since logit cannot model true zero

    max_se : float, optional
      the standard error to use for data with missing or zero standard error
    """
    est_mesh = dm.get_estimate_age_mesh()

    # get the index vector and weight vector for the age range
    age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
    age_weights = d.get('age_weights', np.ones(len(age_indices)))

    # ensure all rate data is valid
    d_val = dm.value_per_1(d)
    if d_val < 0 or d_val > 1:
        debug('WARNING: data %d not in range (0,1)' % d['id'])
        raise ValueError
    elif d_val == 0.:
        d_val = min_val / 10.  # TODO: determine if this is an acceptible way to deal with zero
    elif d_val == 1.:
        d_val = 1. - min_val / 10.

    logit_val = mc.logit(d_val)

    d_se = dm.se_per_1(d)
    if d_se == MISSING:
        d_se = max_se * 10. #TODO: determine if this is an acceptible way to deal with missing
    elif d_se == 0.:
        d_se = max_se

    logit_se = (1/d_val + 1/(1-d_val)) * d_se

    return age_indices, age_weights, logit_val, logit_se