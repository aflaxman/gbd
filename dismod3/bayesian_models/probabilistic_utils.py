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

def spline_interpolate(in_mesh, values, out_mesh):
    from scipy.interpolate import interp1d
    f = interp1d(in_mesh, values, kind='linear')
    return f(out_mesh)

# def gp_interpolate(in_mesh, values, out_mesh):
#     """
#     interpolate a set of values given at
#     points on in_mesh to find values on
#     out_mesh.
#     """
#     M,C = uninformative_prior_gp()
#     gp.observe(M,C,in_mesh,values)
#     return M(out_mesh)

def interpolate(in_mesh, values, out_mesh):
    """
    wrapper so that it is only necessary to
    make one change to try different interpolation
    methods
    """
    return spline_interpolate(in_mesh, values, out_mesh)

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

    # use prior to set rate near zero as requested
    for prior_str in asrf.fit.get('priors', '').split('\n'):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])

            gp.observe(M, C, range(age_start, age_end+1), [-10.], [0.])
               
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

def trim(x, a, b):
    return np.maximum(a, np.minimum(b, x))

def flatten(l):
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))
        else:
            out.append(item)
    return out

INV_TRANSFORM = {
    'logit': mc.invlogit,
    'log': np.exp,
    }

TRANSFORM = {
    'logit': mc.logit,
    'log': np.log,
    }

def add_stoch_to_rf_vars(rf, name, initial_value, transform='logit'):
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
    transformed_rate = mc.Normal('%s(%s)' % (transform, name), mu=np.zeros(len(mesh)),
                                 tau=1.e-2, value=transform_func(initial_value[mesh]),
                                 verbose=0)
    # the rate function is obtained by "non-parametric regression"
    # using a Gaussian process with a nice covariance function to fill
    # in the mesh of logit_rate, and then looking at the inverse logit
    @mc.deterministic(name=name)
    def rate(transformed_rate=transformed_rate):
        return interpolate(mesh, inv_transform_func(transformed_rate), out_mesh)

    rf.vars['%s(%s)' % (transform, name)] = transformed_rate
    rf.vars[name] = rate

def add_priors_to_rf_vars(rf):
    """
    include priors specified in rf.params['priors'] in the rf model
    """

    # TODO: refactor the rate function vars structure so that it is
    # more straight-forward where all of these priors are applied

    # rf.fit['priors'] shall have the following format
    # smooth <tau> <age_start> <age_end>
    # zero <age_start> <age_end>
    # confidence <mean> <tau>
    # increasing <age_start> <age_end>
    # decreasing <age_start> <age_end>
    # convex_up <age_start> <age_end>
    # convex_down <age_start> <age_end>
    # unimodal <age_start> <age_end>
    #
    # for example: 'smooth .1 \n zero 0 5 \n zero 95 100'


    def derivative_sign_prior(rf, prior, deriv, sign):
        age_start = int(prior[1])
        age_end = int(prior[2])
        @mc.potential(name='deriv_sign-%d-%d-%d-%d^%d' % (deriv, sign, age_start, age_end, rf.id))
        def deriv_sign_rate(f=rf.vars['Erf_%d'%rf.id],
                            age_start=age_start, age_end=age_end, tau=1000.,
                            deriv=deriv, sign=sign):
            df = np.diff(f[age_start:(age_end+1)], deriv)
            return -tau * np.dot(df**2, (sign * df < 0))
        return [deriv_sign_rate]

    rf.vars['prior hyper-params'] = []
    rf.vars['prior'] = []
    
    for prior_str in rf.fit.get('priors', '').split('\n'):
        prior = prior_str.split()
        if len(prior) == 0:
            continue
        if prior[0] == 'smooth':
            # tau_smooth_rate = mc.InverseGamma('smooth_rate_tau_%d'%rf.id, .01, .05, value=5.)
            tau_smooth_rate = float(prior[1])

            if len(prior) == 4:
                age_start = int(prior[2])
                age_end = int(prior[3])
            else:
                age_start = 0
                age_end = MAX_AGE
                
            rf.vars['prior hyper-params'] += [tau_smooth_rate]

            @mc.potential(name='smooth-%d-%d^%d'%(age_start, age_end, rf.id))
            def smooth_rate(f=rf.vars['Erf_%d'%rf.id], age_start=age_start, age_end=age_end, tau=tau_smooth_rate):
                return mc.normal_like(np.diff(np.log(np.maximum(NEARLY_ZERO, f[range(age_start, age_end)]))), 0.0, tau)
            rf.vars['prior'] += [smooth_rate]

        elif prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])
                               
            @mc.potential(name='zero-%d-%d^%d' % (age_start, age_end, rf.id))
            def zero_rate(f=rf.vars['Erf_%d'%rf.id], age_start=age_start, age_end=age_end, tau=1./(1e-4)**2):
                return mc.normal_like(f[range(age_start, age_end+1)], 0.0, tau)
            rf.vars['prior'] += [zero_rate]

        elif prior[0] == 'confidence':
            # prior only affects beta_binomial_rate model
            if not rf.vars.has_key('confidence'):
                continue

            mu = float(prior[1])
            tau = float(prior[2])

            @mc.potential(name='conf^%d'%rf.id)
            def confidence(f=rf.vars['confidence'], mu=mu, tau=tau):
                return mc.normal_like(f, mu, tau)
            rf.vars['prior'] += [confidence]

        elif prior[0] == 'increasing':
            rf.vars['prior'] += derivative_sign_prior(rf, prior, deriv=1, sign=1)
        elif prior[0] == 'decreasing':
            rf.vars['prior'] += derivative_sign_prior(rf, prior, deriv=1, sign=-1)
        elif prior[0] == 'convex_down':
            rf.vars['prior'] += derivative_sign_prior(rf, prior, deriv=2, sign=-1)
        elif prior[0] == 'convex_up':
            rf.vars['prior'] += derivative_sign_prior(rf, prior, deriv=2, sign=1)

        elif prior[0] == 'unimodal':
            age_start = int(prior[1])
            age_end = int(prior[2])
            @mc.potential(name='unimodal-%d-%d^%d' % (age_start, age_end, rf.id))
            def unimodal_rate(f=rf.vars['Erf_%d'%rf.id], age_start=age_start, age_end=age_end, tau=1000.):
                df = np.diff(f[age_start:(age_end + 1)])
                sign_changes = pl.find((df[:-1] > NEARLY_ZERO) & (df[1:] < -NEARLY_ZERO))
                sign = np.ones(age_end-age_start-1)
                if len(sign_changes) > 0:
                    change_age = sign_changes[len(sign_changes)/2]
                    sign[change_age:] = -1.
                return -tau*np.dot(np.abs(df[:-1]), (sign * df[:-1] < 0))
            rf.vars['prior'] += [unimodal_rate]

        else:
            raise KeyException, 'Unrecognized prior: %s' % prior_str

def save_map(asrf):
    asrf.fit['map'] = list(asrf.map_fit_stoch.value)
    asrf.save()

def save_mcmc(asrf):
    rate = asrf.mcmc_fit_stoch.trace()
    trace_len = len(rate)
    
    sr = []
    for ii in asrf.fit['out_age_mesh']:
        sr.append(sorted(rate[:,ii]))
    asrf.fit['mcmc_lower_cl'] = [sr[ii][int(.025*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_median'] = [sr[ii][int(.5*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_upper_cl'] = [sr[ii][int(.975*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_mean'] = list(np.mean(rate, 0))

    asrf.save()
