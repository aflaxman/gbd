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

import numpy as np
import pymc as mc
import sys

import dismod3
from dismod3.utils import debug, interpolate, rate_for_range, indices_for_range, generate_prior_potentials, gbd_regions, clean, type_region_year_sex_from_key
from dismod3.settings import MISSING, NEARLY_ZERO

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

    data = [d for d in dm.data if clean(d['data_type']).find(param_type) != -1 and not d.get('ignore')]
    dm.calc_effective_sample_size(data)

    dm.clear_empirical_prior()
    dm.fit_initial_estimate(param_type, data)

    # don't do anything if there is no data for this parameter type
    if len(data) == 0:
        return

    dm.vars = setup(dm, param_type, data)
    print 'i', '%s' % ', '.join(['%.2f' % x for x in dm.get_initial_value(param_type)[::10]])
    sys.stdout.flush()
    
    # fit the model
    dm.map = mc.MAP(dm.vars)
    try:
        dm.map.fit(method='fmin_powell', iterlim=500, tol=.1, verbose=1)
    except KeyboardInterrupt:
        print 'User halted optimization routine before optimal value found'
    sys.stdout.flush()

    dm.mcmc = mc.MCMC(dm.vars)
    dm.mcmc.sample(10000)

    dm.vars['region_coeffs'].value = dm.vars['region_coeffs'].stats()['mean']
    dm.vars['study_coeffs'].value = dm.vars['study_coeffs'].stats()['mean']
    dm.vars['age_coeffs_mesh'].value = dm.vars['age_coeffs_mesh'].stats()['mean']
    dm.vars['log_dispersion'].value = dm.vars['log_dispersion'].stats()['mean']

    alpha = dm.vars['region_coeffs'].stats()['mean']
    beta = dm.vars['study_coeffs'].stats()['mean']
    gamma_mesh = dm.vars['age_coeffs_mesh'].stats()['mean']
    print 'a', '%s' % ', '.join(['%.2f' % x for x in alpha])
    print 'b', '%s' % ', '.join(['%.2f' % x for x in beta])
    print 'g', '%s' % ', '.join(['%.2f' % x for x in gamma_mesh])
    print 'd', '%.2f' % dm.vars['dispersion'].stats()['mean']
    print 'm', '%s' % ', '.join(['%.2f' % x for x in dm.vars['rate_stoch'].stats()['mean'][::10]])
    X = covariates(data[0])
    #print X
    print 'p', '%s' % ', '.join(['%.2f' % x for x in predict_rate(X, alpha, beta, gamma_mesh)])
    # save the results in the param_hash
    prior_vals = dict(
        alpha=list(dm.vars['region_coeffs'].stats()['mean']),
        beta=list(dm.vars['study_coeffs'].stats()['mean']),
        gamma=list(dm.vars['age_coeffs'].stats()['mean']),
        delta=float(dm.vars['dispersion'].stats()['mean']))

    prior_vals.update(
        sigma_alpha=list(dm.vars['region_coeffs'].stats()['standard deviation']),
        sigma_beta=list(dm.vars['study_coeffs'].stats()['standard deviation']),
        sigma_gamma=list(dm.vars['age_coeffs'].stats()['standard deviation']),
        sigma_delta=float(dm.vars['dispersion'].stats()['standard deviation']))
    dm.set_empirical_prior(param_type, prior_vals)

    dispersion = prior_vals['delta']
    for r in dismod3.gbd_regions:
        for y in dismod3.gbd_years:
            for s in dismod3.gbd_sexes:
                key = dismod3.gbd_key_for(param_type, r, y, s)
                mu = predict_rate(regional_covariates(key),
                                  alpha=prior_vals['alpha'],
                                  beta=prior_vals['beta'],
                                  gamma=prior_vals['gamma'])
                dm.set_initial_value(key, mu)
                dm.set_mcmc('emp_prior_mean', key, mu)

                emp_p = mc.NegativeBinomial('emp_p', mu*1000, dispersion)
                mc.MCMC([emp_p]).sample(25)
                dm.set_mcmc('emp_prior_lower_ui', key, emp_p.stats()['quantiles'][2.5]/1000)
                dm.set_mcmc('emp_prior_upper_ui', key, emp_p.stats()['quantiles'][97.5]/1000)

    key = dismod3.gbd_key_for(param_type, 'world', 1997, 'total')
    mu = predict_rate(regional_covariates(key),
                      alpha=prior_vals['alpha'],
                      beta=prior_vals['beta'],
                      gamma=prior_vals['gamma'])
    dm.set_initial_value(key, mu)
    dm.set_mcmc('emp_prior_mean', key, mu)
    emp_p = mc.NegativeBinomial('emp_p', mu*84000, dispersion)
    mc.MCMC([emp_p]).sample(25)
    dm.set_mcmc('emp_prior_lower_ui', key, emp_p.stats()['quantiles'][2.5]/84000)
    dm.set_mcmc('emp_prior_upper_ui', key, emp_p.stats()['quantiles'][97.5]/84000)

from logit_normal_model import covariates, regional_covariates

def predict_rate(X, alpha, beta, gamma):
    """ Calculate logit(Y) = gamma + X * beta"""
    Xa, Xb = X
    return np.exp(np.dot(Xa, alpha) + np.dot(Xb, beta) + gamma)

def setup(dm, key, data_list, rate_stoch=None, emp_prior={}):
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

    dm.calc_effective_sample_size(data_list)

    # for debugging
    #if key == 'incidence+asia_southeast+1990+female':

    # generate regional covariates
    X_region, X_study = regional_covariates(key)

    # use the empirical prior mean if it is available
    if len(set(emp_prior.keys()) & set(['alpha', 'beta', 'gamma', 'delta'])) == 4:
        mu_alpha = np.array(emp_prior['alpha'])
        sigma_alpha = max([.1] + emp_prior['sigma_alpha'])
        alpha = np.array(emp_prior['alpha'])

        beta = np.array(emp_prior['beta'])
        sigma_beta = max([.1] + emp_prior['sigma_beta'])

        mu_gamma = np.array(emp_prior['gamma'])
        sigma_gamma = max([.1] + emp_prior['sigma_gamma'])

        mu_delta = emp_prior['delta']
        sigma_delta = emp_prior['sigma_delta']

    else:
        mu_alpha = np.zeros(len(X_region))
        sigma_alpha = .5
        alpha = mc.Normal('region_coeffs_%s' % key, mu=mu_alpha, tau=sigma_alpha**-2., value=mu_alpha)
        vars.update(region_coeffs=alpha)

        mu_beta = np.zeros(len(X_study))
        sigma_beta = .1
        beta = mc.Normal('study_coeffs_%s' % key, mu=mu_beta, tau=sigma_beta**-2., value=mu_beta)
        vars.update(study_coeffs=beta)

        mu_gamma = -5.*np.ones(len(est_mesh))
        sigma_gamma = 5.

        mu_delta = 100.
        sigma_delta = 1.

    log_delta = mc.Uninformative('log(dispersion_%s)' % key, value=np.log(mu_delta - 1.))
    delta = mc.Lambda('dispersion_%s' % key, lambda x=log_delta: 1. + np.exp(x))
    @mc.potential(name='dispersion_potential_%s' % key)
    def delta_potential(delta=delta, mu=mu_delta, tau=sigma_delta**-2):
        return mc.normal_like(delta, mu, tau)
    vars.update(log_dispersion=log_delta,
                dispersion=delta,
                dispersion_potential=delta_potential)


    # create varible for interpolated rate;
    # also create variable for age-specific rate function, if it does not yet exist
    if rate_stoch:
        # if the rate_stoch already exists, for example prevalence in the generic model,
        # we use it to back-calculate mu and eventually gamma
        mu = rate_stoch

        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(mu=mu, Xa=X_region, Xb=X_study, alpha=alpha, beta=beta):
            return np.log(mu) - np.dot(alpha, Xa) - np.dot(beta, Xb)

        @mc.potential(name='age_coeffs_potential_%s' % key)
        def gamma_potential(gamma=gamma, mu_gamma=mu_gamma, tau_gamma=1./sigma_gamma**2, param_mesh=param_mesh):
            return mc.normal_like(gamma[param_mesh], mu_gamma[param_mesh], tau_gamma)

        vars.update(rate_stoch=mu, age_coeffs=gamma, age_coeffs_potential=gamma_potential)
        
    else:
        # if the rate_stoch does not yet exists, we make gamma a stoch, and use it to calculate mu
        # for computational efficiency, gamma is a linearly interpolated version of gamma_mesh
        initial_gamma = np.log(np.maximum(dm.get_initial_value(key), NEARLY_ZERO))
        gamma_mesh = mc.Normal('age_coeffs_mesh_%s' % key, mu=mu_gamma[param_mesh], tau=sigma_gamma**-2, value=initial_gamma[param_mesh])
        
        @mc.deterministic(name='age_coeffs_%s' % key)
        def gamma(gamma_mesh=gamma_mesh, param_mesh=param_mesh, est_mesh=est_mesh):
            return interpolate(param_mesh, gamma_mesh, est_mesh)

        @mc.deterministic(name=key)
        def mu(Xa=X_region, Xb=X_study, alpha=alpha, beta=beta, gamma=gamma):
            return predict_rate([Xa, Xb], alpha, beta, gamma)

        vars.update(age_coeffs_mesh=gamma_mesh, age_coeffs=gamma, rate_stoch=mu)


    # create potentials for priors
    vars['priors'] = generate_prior_potentials(dm.get_priors(key), est_mesh, mu)
    
    
    # create observed stochastics for data
    vars['data'] = data_list
    vars['observed_rates'] = []

    for d in data_list:
        try:
            age_indices, age_weights, Y_i, N_i = values_from(dm, d)
        except ValueError:
            print 'WARNING: could not calculate likelihood for data %d' % d['id']

        @mc.observed
        @mc.stochastic(name='data_%d' % d['id'])
        def obs(value=Y_i*N_i, N_i=N_i,
                X=covariates(d),
                alpha=alpha, beta=beta, gamma=gamma, delta=delta,
                age_indices=age_indices,
                age_weights=age_weights):

            # calculate study-specific rate function
            mu = predict_rate(X, alpha, beta, gamma)
            mu_i = rate_for_range(mu, age_indices, age_weights)
            logp = mc.negative_binomial_like(value, mu_i*N_i, delta)
            return logp
            
        vars['observed_rates'].append(obs)
        
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
    age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
    age_weights = d.get('age_weights', np.ones(len(age_indices))/len(age_indices))

    # ensure all rate data is valid
    Y_i = dm.value_per_1(d)
    # TODO: allow Y_i > 1, extract effective sample size appropriately in this case
    if Y_i < 0:
        debug('WARNING: data %d < 0' % d['id'])
        raise ValueError

    N_i = max(1000., d['effective_sample_size'])
    print N_i
    
    return age_indices, age_weights, Y_i, N_i
