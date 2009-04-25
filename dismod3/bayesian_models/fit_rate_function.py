import inspect

import pymc as mc

import probabilistic_utils
#import beta_binomial_rate as rate_model
import urbanicity_covariate_rate as rate_model

MAP_PARAMS = {
    'bfgs': [ 500, 'fmin_l_bfgs_b'],
    'powell': [ 500, 'fmin_powell'],
    'most accurate': [ 1500, 'fmin_powell'],
    'fast': [ 50, 'fmin_powell'],
    'testing fast': [ 1, 'fmin' ],
    }

MCMC_PARAMS = {
    'over accurate': [500, 200, 10000],
    'most accurate': [100, 50, 5000],
    'fast': [500, 10, 5000],
    'testing fast': [5, 5, 5],
    }

def initialize_model(asrf):
    # store the rate model code in the asrf for future reference
    asrf.fit['rate_model'] = inspect.getsource(rate_model)
    asrf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = probabilistic_utils.normal_approx(asrf)

    # define the model variables
    rate_model.setup_rate_model(asrf)

def map_fit(asrf, speed='most accurate'):
    """
    The Maximum A Posteriori (MAP) fit of the model is a point
    estimate of the model parameters, which is found using numerical
    optimization to attempt to maximize the posterior liklihood.
    Since this is a local optimization method, it might not find the
    global optimum.
    """
    initialize_model(asrf)
    
    map = mc.MAP(asrf.vars)
    iterlim, method = MAP_PARAMS[speed]
    print "searching for maximum likelihood point estimate (%s)" % method
    map.fit(verbose=1, iterlim=iterlim, method=method)

    probabilistic_utils.save_map(asrf)

def mcmc_fit(asrf, speed='most accurate'):
    """
    The Markov Chain Monte Carlo (MCMC) fit of the model works by
    making successive draws of the model parameters from the posterior
    distribution.  This provides confidence intervals, and should be
    more robust against local maxima in the posterior liklihood.  But
    the question is, did the chain run for long enough to mix?
    """
    map_fit(asrf, speed)
    
    print "drawing samples from posterior distribution (MCMC) (speed: %s)" % speed
    mcmc = mc.MCMC(asrf.vars)
    if asrf.vars.has_key('beta_binom_stochs'):
        mcmc.use_step_method(mc.AdaptiveMetropolis, asrf.vars['beta_binom_stochs'], verbose=0)
    if asrf.vars.has_key('logit(Erf_%d)'%asrf.id):
        mcmc.use_step_method(mc.AdaptiveMetropolis, asrf.vars['logit(Erf_%d)'%asrf.id], verbose=0)

    trace_len, thin, burn = MCMC_PARAMS[speed]
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    probabilistic_utils.save_mcmc(asrf)


    
import twill.commands as twc
import simplejson as json
from dismod3.settings import *

def get_disease_model(disease_model_id):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_LOGIN_URL)
    twc.fv('1', 'username', DISMOD_USERNAME)
    twc.fv('1', 'password', DISMOD_PASSWORD)
    twc.submit()
    twc.url('accounts/profile')

    twc.go(DISMOD_DOWNLOAD_URL % disease_model_id)

    return json.loads(twc.show())

def post_disease_model(disease_model):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_LOGIN_URL)
    twc.fv('1', 'username', DISMOD_USERNAME)
    twc.fv('1', 'password', DISMOD_PASSWORD)
    twc.submit()
    twc.url('accounts/profile')

    twc.go(DISMOD_UPLOAD_URL)
    twc.fv('1', 'model_json', json.dumps(disease_model))
    twc.submit()

    return twc.browser.get_url()


def fit(disease_model, data_type='prevalence data'):
    """
    download a disease model with twill, and fit just the
    data corresponding to the specified data_type, and upload
    the results as a new model
    """
    dm = get_disease_model(disease_model)

    # filter out all data with type != data_type
    dm['data'] = [[i,d] for [i,d] in dm['data'] if d['data_type'] == data_type]

    # store the probabilistic model code for future reference
    dm['params']['bayesian_model'] = inspect.getsource(rate_model)
    dm['params']['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    fit_normal_approx(dm, data_type)
    
    # define the model variables
    # rate_model.setup_rate_model(dm['data'])
    
    return dm


def fit_normal_approx(dm, data_type):
    """
    This 'normal approximation' estimate for an age-specific dataset
    is formed by using each datum to produce an estimate of the
    function value at a single age, and then saying that the logit of
    the true rate function is a gaussian process and these
    single age estimates are observations of this gaussian process.

    This is less valid and less accurate than using MCMC or MAP, but
    it is much faster.  It is used to generate an initial value for
    the maximum-liklihood estimate.
    """
    from pymc import gp
    from probabilistic_utils import uninformative_prior_gp, NEARLY_ZERO, MAX_AGE
    MISSING = -99
    

    param_hash = dm['params']
    data_list = [d for i,d in dm['data'] if d['data_type'] == data_type]

    M,C = uninformative_prior_gp()

    age = []
    val = []
    V = []
    for d in data_list:
        scale = float(d['units'].split()[-1])

        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE

        if d['standard_error'] == 0.:
            d['standard_error'] = .001

        age.append(.5 * (d['age_start'] + d['age_end']))
        val.append(d['value'] / scale + .00001)
        V.append((d['standard_error'] / scale) ** 2.)

    if len(data_list) > 0:
        gp.observe(M, C, age, mc.logit(val), V)

    # use prior to set rate near zero as requested
    near_zero = min(1., val)**2
    if near_zero == 1.:
        near_zero = 1e-9
        
    for prior_str in param_hash.get('priors', '').split('\n'):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])

            gp.observe(M, C, range(age_start, age_end+1), mc.logit(near_zero), [0.])
        
    x = param_hash['out_age_mesh']
    normal_approx_vals = mc.invlogit(M(x))
    param_hash['normal_approx'] = list(normal_approx_vals)


