""" Logit-Normal Model for a generic epidemological parameter

The Logit-Normal Model represents data for age-specific rates
according to the following formula::

    logit(Y_i) ~ \sum _{a = a_{i0}} ^{a_{i1}} w(a) \mu_{i, a} + N(0, \sigma_i^2 + \sigma^2)
    \mu_{i, a} = \alpha_{r_i} + \alpha_{s_i} + \alpha_{y_i} + \gamma_a + \beta^T X_i

Here Y_i, \sigma_i, a_{i0}, a_{i1}, r_i, s_i, y_i, and X_i are the
value, standard error in logit-space, age range, region, sex, year,
and study-level covariates corresponding to a single age-range value
from a single study.  \alpha, \beta, \gamma, and \sigma are parameters
(fixed effects and random effect) that will be estimated from the
data.
"""

import numpy as np
import pymc as mc

import dismod3
from dismod3.utils import debug, interpolate, rate_for_range, indices_for_range, generate_prior_potentials, gbd_regions, clean, type_region_year_sex_from_key
from dismod3.settings import MISSING

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
    prior_vals = dict(
        alpha=list(dm.vars['region_coeffs'].value),
        beta=list(dm.vars['study_coeffs'].value),
        gamma=list(dm.vars['age_coeffs'].value),
        sigma=float(dm.vars['dispersion'].value))
    dm.set_empirical_prior(param_type, prior_vals)

    dispersion = prior_vals['sigma']
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for(param_type, r, y, s)
                logit_mu = predict_logit_rate(regional_covariates(key), **prior_vals)
                mu = mc.invlogit(logit_mu)
                dm.set_initial_value(key, mu)
                dm.set_mcmc('emp_prior_mean', key, mu)
                dm.set_mcmc('emp_prior_lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
                dm.set_mcmc('emp_prior_upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

    key = dismod3.gbd_key_for(param_type, 'world', 1997, 'total')
    logit_mu = predict_logit_rate(regional_covariates(key), **prior_vals)
    mu = mc.invlogit(logit_mu)
    dm.set_initial_value(key, mu)
    dm.set_mcmc('emp_prior_mean', key, mu)
    dm.set_mcmc('emp_prior_lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
    dm.set_mcmc('emp_prior_upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

def covariates(d):
    """ extract the covariates from a data point as a vector;

    Xa represents region-level covariates:
      Xa[0],...,Xa[21] = region indicators
      Xa[22] = year-1997
      Xa[23] = 1 if sex == 'male', -1 if sex == 'female'
    Xb represented study-level covariates:
      Xb[0] = self-reported
      Xb[1] = threshold (integer)
    """
    Xa = np.zeros(len(gbd_regions) + 2)
    for ii, r in enumerate(gbd_regions):
        if clean(d['gbd_region']) == clean(r):
            Xa[ii] = 1.

    Xa[ii+1] = .1 * .5 * (float(d['year_start']) + float(d['year_end'])) - 1997

    if clean(d['sex']) == 'male':
        Xa[ii+2] = .5
    elif clean(d['sex']) == 'female':
        Xa[ii+2] = -.5
    else:
        Xa[ii+2] = 0.

    Xb = np.zeros(5.)

    # TODO: instead of hard-coding this, store it in the disease model
    # (and let users set it through the web)
    if clean(d.get('self_reported', '')) == 'true':
        Xb[0] = 1.
    if d.has_key('threshold'):
        Xb[0] = float(d['threshold'])
        
    return Xa, Xb

def regional_covariates(key):
    """ form the covariates for a gbd key"""
    t,r,y,s = type_region_year_sex_from_key(key)

    d = {'gbd_region': r,
         'year_start': y,
         'year_end': y,
         'sex': s}
    return covariates(d)

def predict_logit_rate(X, alpha, beta, gamma, sigma=0):
    """ Calculate logit(Y) = gamma + X * beta (sigma is unused, included for hacky convenience)"""
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
    """
    vars = {}
    est_mesh = dm.get_estimate_age_mesh()
    param_mesh = dm.get_param_age_mesh()

    if np.any(np.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    # for debugging
    #if key == 'incidence+asia_southeast+1990+female':
    #    import pdb; pdb.set_trace()

    # generate regional covariates
    X_region, X_study = regional_covariates(key)

    # use the empirical prior mean if it is available
    if set(emp_prior.keys()) == set(['alpha', 'beta', 'gamma', 'sigma']):
        mu_alpha = np.array(emp_prior['alpha'])
        sigma_alpha = .01

        beta = np.array(emp_prior['beta'])

        mu_gamma = np.array(emp_prior['gamma'])
        sigma_gamma = 1. #emp_prior['sigma']

        mu_sigma = .1
        conf_sigma = 10.
    else:
        mu_alpha = np.zeros(len(X_region))
        sigma_alpha = .1

        mu_beta = np.zeros(len(X_study))
        sigma_beta = .01
        beta = mc.Normal('study_coeffs_%s' % key, mu=mu_beta, tau=1/sigma_beta**2, value=mu_beta)
        vars.update(study_coeffs=beta)

        mu_gamma = -5.*np.ones(len(est_mesh))
        sigma_gamma = 1.

        mu_sigma = .1
        conf_sigma = 10.

    alpha = mc.Normal('region_coeffs_%s' % key, mu=mu_alpha, tau=1/sigma_alpha**2, value=mu_alpha)
    vars.update(region_coeffs=alpha)


    log_sigma = mc.Uninformative('log(dispersion_%s)' % key, value=np.log(mu_sigma))
    @mc.deterministic(name='dispersion_%s' % key)
    def sigma(log_sigma=log_sigma):
        return np.exp(log_sigma)
    # TODO: replace this potential in the generate_prior_potentials function if confidence is set
    @mc.potential(name='dispersion_potential_%s' % key)
    def sigma_potential(sigma=sigma, alpha=conf_sigma, beta=conf_sigma/mu_sigma):
        return mc.gamma_like(sigma, alpha, beta)
    vars.update(log_dispersion=log_sigma,
                dispersion=sigma,
                dispersion_potential=sigma_potential)


    # create varible for interpolated logit rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        # if the rate_stoch already exists, for example prevalence in the generic model,
        # we use it to back-calculate mu and eventually gamma
        @mc.deterministic(name='logit_%s' % key)
        def mu(invlogit_mu=rate_stoch):
            return mc.logit(invlogit_mu)

        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(mu=mu, Xa=X_region, Xb=X_study, alpha=alpha, beta=beta):
            return mu - np.dot(alpha, Xa) - np.dot(beta, Xb)

        @mc.potential(name='age_coeffs_potential_%s' % key)
        def gamma_potential(gamma=gamma, mu_gamma=mu_gamma, tau_gamma=1./sigma_gamma**2, param_mesh=param_mesh):
            return mc.normal_like(gamma[param_mesh], mu_gamma[param_mesh], tau_gamma)

        vars.update(rate_stoch=rate_stoch, logit_rate_stoch=mu, age_coeffs=gamma, age_coeffs_potential=gamma_potential)
        
    else:
        # if the rate_stoch does not yet exists, we make gamma a stoch, and use it to calculate mu
        # for computational efficiency, gamma is a linearly interpolated version of gamma_mesh
        initial_gamma = mu_gamma
        gamma_mesh = mc.Normal('age_coeffs_mesh_%s' % key, mu=mu_gamma[param_mesh], tau=1/sigma_gamma**2, value=initial_gamma[param_mesh])
        
        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(gamma_mesh=gamma_mesh, param_mesh=param_mesh, est_mesh=est_mesh):
            return interpolate(param_mesh, gamma_mesh, est_mesh)

        @mc.deterministic(name='logit_%s' % key)
        def mu(Xa=X_region, Xb=X_study, alpha=alpha, beta=beta, gamma=gamma):
            return np.dot(alpha, Xa) + np.dot(beta, Xb) + gamma

        @mc.deterministic(name=key)
        def rate_stoch(mu=mu):
            return mc.invlogit(mu)

        vars.update(age_coeffs_mesh=gamma_mesh, age_coeffs=gamma, logit_rate_stoch=mu, rate_stoch=rate_stoch)


    # create potentials for priors
    vars['priors'] = generate_prior_potentials(dm.get_priors(key), est_mesh, rate_stoch)
    
    
    # create observed stochastics for data
    vars['data'] = data_list
    vars['observed_rates'] = []

    min_val = min([1.e-9] + [dm.value_per_1(d) for d in data_list if dm.value_per_1(d) > 0]) # TODO: assess validity of this minimum value
    max_se = max([.000001] + [dm.se_per_1(d) for d in data_list if dm.se_per_1(d) > 0])  # TODO: assess validity of this maximum std err

    #import pdb; pdb.set_trace()
    for d in data_list:
        try:
            age_indices, age_weights, logit_val, logit_se = values_from(dm, d, min_val, max_se)
        except ValueError:
            continue

        @mc.observed
        @mc.stochastic(name='data_%d' % d['id'])
        def obs(value=logit_val, logit_se=logit_se,
                X=covariates(d),
                alpha=alpha, beta=beta, gamma=gamma, sigma=sigma,
                age_indices=age_indices,
                age_weights=age_weights):

            # calculate study-specific rate function
            mu = predict_logit_rate(X, alpha, beta, gamma)
            mu_i = rate_for_range(mu, age_indices, age_weights)
            
            tau_i = 1. / (sigma**2 + logit_se**2)
            logp = mc.normal_like(x=value, mu=mu_i, tau=tau_i)
            return logp
            
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
        d_se = max_se #TODO: determine if this is an acceptible way to deal with missing
    elif d_se == 0.:
        d_se = max_se

    logit_se = (1/d_val + 1/(1-d_val)) * d_se

    return age_indices, age_weights, logit_val, logit_se
