import inspect

import pymc as mc

import probabilistic_utils
import rate_beta_binomial as rate_model

def map_fit(asrf, speed='most accurate'):
    """
    The Maximum A Posteriori (MAP) fit of the model is a point
    estimate of the model parameters, which is found using numerical
    optimization to attempt to maximize the posterior liklihood.
    Since this is a local optimization method, it might not find the
    global optimum.
    """
    # store the rate model code in the asrf for future reference
    asrf.fit['rate_model'] = inspect.getsource(rate_model)
    asrf.fit['out_age_mesh'] = range(probabilistic_utils.MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = probabilistic_utils.normal_approx(asrf)

    # define the model variables
    rate_model.setup_rate_model(asrf)
    map = mc.MAP(asrf.vars)
    print "searching for maximum likelihood point estimate"
    if speed == 'most accurate':
        iterlim, method = 500, 'fmin_powell'
    elif speed == 'fast':
        iterlim, method = 25, 'fmin_powell'
    elif speed == 'testing fast':
        iterlim, method = 1, 'fmin'

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

    # TODO: make these part of the rate_model, since different models
    # will require more or less burn-in and thinning
    if speed == 'most accurate':
        trace_len, thin, burn = 1000, 100, 10000
    elif speed == 'fast':
        trace_len, thin, burn = 500, 10, 5000
    elif speed == 'testing fast':
        trace_len, thin, burn = 10, 1, 1

    mcmc = mc.MCMC(asrf.vars)
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=2)

    probabilistic_utils.save_mcmc(asrf)
