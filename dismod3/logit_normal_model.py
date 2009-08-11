""" Logit-Normal Model for a generic epidemological parameter

The Logit-Normal Model represents data for age-specific risks/ratios
according to the following formula::

    logit(Y_i) ~ \sum _{a = a_{i0}} ^{a_{i1}} w(a) \mu_{i, a} + N(0, \sigma_i^2 + \sigma_d^2)
    \mu_{i, a} = \gamma_a + \beta^T X_i

Here Y_i, \sigma_i, a_{i0}, a_{i1}, and X_i are the value, standard
error, age range, and covariates corresponding to a single age-range
value from a single study.  \beta, \gamma, and \sigma_d are parameters
that will be estimated from the data.

"""

import numpy as np
import pymc as mc
import random

import dismod3
from dismod3.utils import debug, interpolate, rate_for_range, indices_for_range, generate_prior_potentials, gbd_regions, clean
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

    data =  [d for d in dm.data if clean(d['data_type']).find(param_type) != -1]
    # use a random subset of the data if there is a lot of it,
    # to speed things up
    if len(data) > 25:
        dm.fit_initial_estimate(param_type, random.sample(data,25))
    else:
        dm.fit_initial_estimate(param_type, data)

    dm.vars = setup(dm, param_type, data)
    
    # fit the model
    dm.map = mc.MAP(dm.vars)
    try:
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.0001, verbose=1)
    except KeyboardInterrupt:
        print 'User halted optimization routine before optimal value found'

    #print 'coefficient values: ', dm.vars['coefficients'].value
    #import pdb; pdb.set_trace()
    
    # save the results in the param_hash
    dm.clear_empirical_prior()
    beta = dm.vars['coefficients'].value
    gamma = dm.vars['interp_logit_rate'].value
    dispersion = dm.vars['dispersion'].value
    dm.set_empirical_prior(param_type, {'beta': list(beta),
                                        'gamma': list(gamma),
                                        'dispersion': float(dispersion)})

    # TODO: make an easier-to-understand way to see the results of the
    # empirical prior fit
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for(param_type, r, y, s)
                logit_mu = predict_logit_risk(regional_covariates(r), beta, gamma)
                mu = mc.invlogit(logit_mu)
                dm.set_initial_value(key, mu)
                dm.set_map(key, mu)
                dm.set_mcmc('lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
                dm.set_mcmc('upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

    key = dismod3.gbd_key_for(param_type, 'world', 'total', 'total')
    logit_mu = predict_logit_risk(regional_covariates('world'), beta, gamma)
    mu = mc.invlogit(logit_mu)
    dm.set_initial_value(key, mu)
    dm.set_map(key, mu)
    dm.set_mcmc('lower_ui', key, mc.invlogit(logit_mu - 1.96*dispersion))
    dm.set_mcmc('upper_ui', key, mc.invlogit(logit_mu + 1.96*dispersion))

def covariates(d):
    """ extract the covariates from a data point as a vector"""
    X = np.zeros(len(gbd_regions))
    for ii, r in enumerate(gbd_regions):
        if clean(str(d['gbd_region'])) == clean(r):
            X[ii] = 1.
    return X

def regional_covariates(r='world'):
    """ form the covariates for a region r"""
    d = {'gbd_region': r}
    return covariates(d)

def predict_risk(X, beta, gamma):
    return mc.invlogit(predict_logit_risk(X, beta, gamma))

def predict_logit_risk(X, beta, gamma):
    return gamma + np.dot(X, beta)

def setup(dm, key, data_list, rate_stoch=None, emp_prior={}, r_cov=regional_covariates('world')):
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

    r_cov : dict, optional
      the covariates to use when predicting the rate_stoch (which is
      connected into the larger model, etc)

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
    # FIXME: names should be more informative!
    if emp_prior.has_key('gamma'):
        gamma = np.array(emp_prior['gamma'])
    else:
        gamma = -5.*np.ones(len(est_mesh))

    logit_mu = gamma
    mu = mc.invlogit(gamma)

    # use the empirical prior covariates if available
    if emp_prior.has_key('beta'):
        X = regional_covariates(r_cov)
        beta = np.array(emp_prior['beta'])
        logit_mu = predict_logit_risk(X, beta, gamma)
        mu = mc.invlogit(logit_mu)

    # use the empirical prior standard error if it is available
    if emp_prior.has_key('dispersion'):
        sigma_d = emp_prior['dispersion']
    else:
        sigma_d = 10.


    # create varible for interpolated logit rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        vars['rate_stoch'] = rate_stoch

        @mc.deterministic(name='interp_logit(%s)' % key)
        def interp_logit_rate(rate_stoch=rate_stoch):
            return mc.logit(rate_stoch)
        vars['interp_logit_rate'] = interp_logit_rate

        # FIXME: if est_mesh doesn't start at 0, x[param_mesh] != interpolate(est_mesh, logit_mu, param_mesh)
        @mc.potential(name='interp_logit_empirical_prior_potential_%s' % key)
        def emp_prior_potential(x=interp_logit_rate, mu=logit_mu, tau=1./sigma_d**2, mesh=param_mesh):
            return mc.normal_like(x[param_mesh], mu[param_mesh], tau)
        vars['rate_emp_prior_potential'] = emp_prior_potential
        
    else:
        # find the logit of the emp prior and initial values, which is
        # a little bit of work because initial values are sampled from
        # the est_mesh, but the logit_initial_values are needed on the
        # param_mesh
        initial_value = dm.get_initial_value(key)
        logit_initial_value = mc.logit(
            interpolate(est_mesh, initial_value, param_mesh))

        param_logit_mu = interpolate(est_mesh, logit_mu, param_mesh)
        logit_rate = mc.Normal('logit(%s)' % key,
                               mu=param_logit_mu,
                               tau=1./sigma_d**2,
                               value=logit_initial_value)
        vars['logit_rate'] = logit_rate

        @mc.deterministic(name='interp_logit(%s)' % key)
        def interp_logit_rate(logit_rate=logit_rate):
            return interpolate(param_mesh, logit_rate, est_mesh)
        vars['interp_logit_rate'] = interp_logit_rate

        @mc.deterministic(name=key)
        def rate_stoch(interp_logit_rate=interp_logit_rate):
            return mc.invlogit(interp_logit_rate)
        vars['rate_stoch'] = rate_stoch


    # create stochastic variable for dispersion/"random effect"
    if emp_prior.has_key('dispersion'):
        mu_dispersion = emp_prior['dispersion']
    else:
        mu_dispersion = .1
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


    # create covariate coefficient stoch
    mu_coefficients = emp_prior.get('coefficients', np.zeros(len(regional_covariates())))
    coefficients = mc.Normal('coefficients_%s' % key, mu_coefficients, 1.e2)
    vars['coefficients'] = coefficients
    
    
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
                logit_rate=interp_logit_rate,
                logit_se=logit_se,
                dispersion=dispersion,
                age_indices=age_indices,
                age_weights=age_weights,
                beta=coefficients,
                X=covariates(d)):
            mean_val = rate_for_range(predict_logit_risk(X, beta, logit_rate), age_indices, age_weights)
            return mc.normal_like(x=value, mu=mean_val, tau=1. / (dispersion**2 + logit_se**2))
            
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
