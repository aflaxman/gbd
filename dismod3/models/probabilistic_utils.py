# Copyright 2008 Abraham Flaxman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.

import numpy as np
import pylab as pl
import pymc as mc
import pymc.gp as gp
import simplejson as json

NEARLY_ZERO = 1.e-10
MAX_AGE = 101

def const_func(x, c):
    """
    useful function for defining a non-informative
    prior on a Gaussian process
    >>> const_func([1,2,3], 17.0)
    [17., 517., 17.]
    """
    return np.zeros(np.shape(x)) + c

def uninformative_prior_gp(c=-10.,  diff_degree=2., amp=100., scale=200.):
    """
    return mean and covariance objects for an uninformative prior on
    the age-specific rate
    """
    M = gp.Mean(const_func, c=c)
    C = gp.Covariance(gp.matern.euclidean, diff_degree=diff_degree,
                      amp=amp, scale=scale)

    return M,C

def gp_interpolate(in_mesh, values, out_mesh):
    """
    interpolate a set of values given at
    points on in_mesh to find values on
    out_mesh.
    """
    M,C = uninformative_prior_gp()
    gp.observe(M,C,in_mesh,values)
    return M(out_mesh)

def rate_for_range(raw_rate,a0,a1,pop_mesh):
    """
    calculate rate for a given age-range,
    using the age-specific population numbers
    given by entries in years t0-t1 of pop_table,
    for given country and sex
    """
    a = range(a0,a1+1)
    age_adjusted_rate = np.sum(raw_rate[a]*pop_mesh)/np.sum(pop_mesh)
    return age_adjusted_rate
        
def logit_rate_from_range(rate):
    """
    calculate age-specific rates and variances
    in logit space from a Rate model object
    """
    logit_mesh = np.arange(rate.age_start, rate.age_end+1)
    pop_vals = np.array(rate.population())

    n = (rate.numerator + 1.*NEARLY_ZERO) * pop_vals / np.sum(pop_vals)
    d = (rate.denominator + 2.*NEARLY_ZERO) * pop_vals / np.sum(pop_vals)
    
    logit_rate = mc.logit(np.minimum(n/d, 1.-NEARLY_ZERO))
    logit_V = ( logit_rate - mc.logit( n/d + (n/d)*(1.-n/d)/np.sqrt(d) ) )**2.

    # filter out the points where the denominator is very close to zero
    good_mesh = []
    good_rate = []
    good_V = []
    
    for ii in range(len(logit_mesh)):
        if n[ii] > 0. and n[ii] < d[ii] and d[ii] > .01:
            good_mesh.append(logit_mesh[ii])
            good_rate.append(logit_rate[ii])
            good_V.append(logit_V[ii])
    
    return good_mesh, good_rate, good_V


def population_during(rate):
    """
    calculate the age-specific population counts
    for years {t0,t0+1,...,t1} of pop_table for
    the specified country and sex
    """
    import numpy as np
    from dismod3.models import Population

    a = range(rate.age_start,rate.age_end+1)
    total = np.zeros(len(a))

    relevant_populations = Population.objects.filter(country=rate.country, sex=rate.sex,
                                                     year__gte=rate.epoch_start, year__lte=rate.epoch_end)
    if relevant_populations.count() == 0:
        print "WARNING: Population for %s not found, using World,%d-%d,%s population instead (rate_id=%d)" \
              % (rate.country, rate.epoch_start, rate.epoch_end, rate.sex, rate.pk)
        relevant_populations = Population.objects.filter(country="World", sex=rate.sex,
                                                         year__gte=rate.epoch_start, year__lte=rate.epoch_end)
        
    for population in relevant_populations:
        M,C = population.gaussian_process()
        total += M(a)

    return np.maximum(NEARLY_ZERO, total/(rate.epoch_end + 1. - rate.epoch_start))

def predict_rate_from_asrf(asrf, observed_rate, fit_type='mcmc_mean'):
    predicted_rate = np.array(asrf.fit[fit_type])
    return rate_for_range(predicted_rate, observed_rate.age_start, observed_rate.age_end, observed_rate.population())



#################### Code for generating a "normal approximation" fit of the data
def normal_approx(asrf):
    """
    This 'normal approximation' of the age-specific rate function is
    formed by using each rate to produce an estimate of the
    age-specific rate, and then saying that that logit of the true
    rate function is a gaussian process and these age-specific rates
    are observations of this gaussian process.

    This is less valid and less accurate than using mcmc or map on the
    vars produced by the model_rate_list method below, but maybe it
    will be faster.
    """
    M,C = uninformative_prior_gp()

    # start with rate near zero at age zero
    gp.observe(M, C, [0.], [-10.], [0.])
               
    for r in asrf.rates.all():
        mesh, obs, V = logit_rate_from_range(r)

        # make sure that there is something to observe
        if mesh == []:
            continue
        
        # uncomment the following line to make more inferences than
        # are valid from the data
        #gp.observe(M, C, mesh, obs, V)

        # uncomment the following 2 lines to make less inferences than
        # possible: it may be better to waste information than have
        # false confidence
        ii = len(mesh)/2
        gp.observe(M, C, [mesh[ii]], [obs[ii]], [V[ii]])

    x = asrf.fit['out_age_mesh']
    na_rate = mc.invlogit(M(x))
    asrf.fit['normal_approx'] = list(na_rate)
    asrf.save()

    return M, C

#################### Code for finding the MAP fit for an age-specific rate function
def map_fit(asrf, speed='most accurate'):
    """
    The Maximum A Posteriori (MAP) fit of the model is a point
    estimate of the model parameters, which is found using numerical
    optimization to attempt to maximize the posterior liklihood.
    Since this is a local optimization method, it might not find the
    global optimum.
    """
    import dismod3.bayesian_models.rate_single_binomial as rate_model

    # store the rate model code in the asrf for future reference
    import inspect
    asrf.fit['rate_model'] = inspect.getsource(rate_model)
    asrf.fit['out_age_mesh'] = range(MAX_AGE)

    # do normal approximation first, to generate a good starting point
    M,C = normal_approx(asrf)

    # define the model variables
    vars = rate_model.model_vars(asrf)
    map = mc.MAP(vars)
    print "searching for maximum likelihood point estimate"
    if speed == 'most accurate':
        iterlim, method = 50, 'fmin_powell'
    elif speed == 'fast':
        iterlim, method = 10, 'fmin_powell'
    elif speed == 'testing fast':
        iterlim, method = 1, 'fmin'

    map.fit(verbose=10, iterlim=iterlim, method=method)
    rate_model.save_map(vars, asrf)

    return vars

#################### Code for finding the MAP fit for an age-specific rate function
def mcmc_fit(asrf, speed='most accurate'):
    """
    The Markov Chain Monte Carlo (MCMC) fit of the model works by
    making successive draws of the model parameters from the posterior
    distribution.  This provides confidence intervals, and should be
    more robust against local maxima in the posterior liklihood.  But
    the question is, did the chain run for long enough to mix?
    """
    import dismod3.bayesian_models.rate_single_binomial as rate_model

    vars = map_fit(asrf, speed)

    print "drawing samples from posterior distribution (MCMC) (speed: %s)" % speed

    if speed == 'most accurate':
        trace_len, thin, burn = 1000, 100, 10000
    elif speed == 'fast':
        trace_len, thin, burn = 500, 10, 5000
    elif speed == 'testing fast':
        trace_len, thin, burn = 10, 1, 1

    mcmc = mc.MCMC(vars)
    mcmc.sample(trace_len*thin+burn, burn, thin, verbose=1)

    rate_model.save_mcmc(vars, asrf)
    return vars

def add_rate_stochs(vars, name, mesh, out_mesh, transform='logit', inv_transform=mc.invlogit):
    """
    generate stochastic random vars for the logit of the age-specific
    rate function called name, measured at points given by mesh, as
    well as its gaussian interpolated inverse logit (i.e. the actual
    rate function)

    save them in the variable dict vars
    """
    # logit_rates have uninformative priors
    #
    # for computational convenience, store values only
    # at mesh points
    transformed_rate = mc.Normal('%s(%s)' % (transform, name), mu = np.zeros(len(mesh)),
                                 tau = 1.e-12, value = -7.*np.ones(len(mesh)))
    # the rate function is obtained by "non-parametric regression"
    # using a Gaussian process with a nice covariance function to fill
    # in the mesh of logit_rate, and then looking at the inverse logit
    #
    # the interpolation is done in logit space to ensure that
    # the rate is always in the interval [0,1]
    @mc.deterministic(name=name)
    def rate(transformed_rate=transformed_rate):
        return inv_transform(gp_interpolate(mesh, transformed_rate, out_mesh))

    vars['%s(%s)' % (transform, name)] = transformed_rate
    vars[name] = rate

def observed_rates_stochs(rates, rate_gp):
    """
    for each rate on the rate_list, set up an observed stochastic variable
    which accounts for the probability of the observation given the true
    rate (and the age-specific population size during the years of the observation)
    
    model the rate observations as a binomial random variables, all independent,
    after conditioning on the rate function.
    """
    vars = []
    for r in rates:
        @mc.observed
        @mc.stochastic(name="rate_%d" % r.id)
        def d_stoc(value=(r.numerator,r.denominator,r.age_start,r.age_end),
                   rate=rate_gp,
                   pop_vals=r.population()):
            n,d,a0,a1 = value
            return mc.binomial_like(x=n, n=d,
                                    p=rate_for_range(rate, a0, a1, pop_vals))
        vars.append(d_stoc)

    return vars

