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
        gp.observe(M, C, mesh, obs, V)

        # uncomment the following 2 lines to make less inferences than
        # possible: it may be better to waste information than have
        # false confidence
        #ii = len(mesh)/2
        #gp.observe(M, C, [mesh[ii]], [obs[ii]], [V[ii]])

    x = asrf.fit['out_age_mesh']
    na_rate = mc.invlogit(M(x))
    asrf.fit['normal_approx'] = list(na_rate)
    asrf.save()

    return M, C

def trim(x, a, b):
    return np.maximum(a, np.minimum(b, x))

INV_TRANSFORM = {
    'logit': mc.invlogit
    }

TRANSFORM = {
    'logit': mc.logit
    }

def add_stochs(rf, name, initial_value, transform='logit'):
    """
    generate stochastic random var, represented in a transformed space at
    points given by rf.fit['age_mesh'], and mapped back to the original space
    by a gaussian interpolated inverse transform, at points given by rf.fit['out_age_mesh']

    save them in rf.vars dictionary
    """
    mesh = rf.fit['age_mesh']
    out_mesh = rf.fit['out_age_mesh']

    inv_transform_func = INV_TRANSFORM[transform]
    transform_func = TRANSFORM[transform]
    # for computational convenience, store values only
    # at mesh points
    transformed_rate = mc.Normal('%s(%s)' % (transform, name), mu = np.zeros(len(mesh)),
                                 tau = 1.e-12, value = transform_func(initial_value[mesh]))
    # the rate function is obtained by "non-parametric regression"
    # using a Gaussian process with a nice covariance function to fill
    # in the mesh of logit_rate, and then looking at the inverse logit
    #
    # the interpolation is done in transformed space to ensure that
    # the rate is always in the image of the inverse transform
    @mc.deterministic(name=name)
    def rate(transformed_rate=transformed_rate):
        return inv_transform_func(gp_interpolate(mesh, transformed_rate, out_mesh))

    rf.vars['%s(%s)' % (transform, name)] = transformed_rate
    rf.vars[name] = rate

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


def save_map(asrf):
    asrf.fit['map'] = list(asrf.rate_stoch.value)
    asrf.save()

def save_mcmc(asrf):
    rate = asrf.rate_stoch.trace()
    trace_len = len(rate)
    
    sr = []
    for ii in asrf.fit['out_age_mesh']:
        sr.append(sorted(rate[:,ii]))
    asrf.fit['mcmc_lower_cl'] = [sr[ii][int(.025*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_median'] = [sr[ii][int(.5*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_upper_cl'] = [sr[ii][int(.975*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_mean'] = list(np.mean(rate, 0))

    asrf.save()
