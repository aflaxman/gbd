import inspect

import numpy as np
import pymc as mc
from pymc import gp
from probabilistic_utils import uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

MISSING = -99

import probabilistic_utils
from probabilistic_utils import trim

import beta_binomial_rate as rate_model
#import urbanicity_covariate_rate as rate_model

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


