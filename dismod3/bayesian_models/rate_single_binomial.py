"""

This Python module uses a Bayesian probability framework, supported by
the PyMC package.

"""

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

import dismod3.models.probabilistic_utils as pu

#################### Code for generating a one-level bayesian model of the data
def add_rate_stochs(vars, name, mesh):
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
    logit_rate = mc.Normal('logit(%s)'%name, mu = np.zeros(len(mesh)),
                           tau = 1.e-12, value = -7.*np.ones(len(mesh)))
    # the rate function is obtained by "non-parametric regression"
    # using a Gaussian process with a nice covariance function to fill
    # in the mesh of logit_rate, and then looking at the inverse logit
    #
    # the interpolation is done in logit space to ensure that
    # the rate is always in the interval [0,1]
    @mc.deterministic(name=name)
    def rate(logit_rate=logit_rate):
        return mc.invlogit(pu.gp_interpolate(mesh,logit_rate))

    vars['logit(%s)'%name] = logit_rate
    vars[name] = rate

def observed_rates_stochs(asrf, rate_gp):
    """
    for each rate on the rate_list, set up an observed stochastic variable
    which accounts for the probability of the observation given the true
    rate (and the age-specific population size during the years of the observation)
    
    model the rate observations as a binomial random variables, all independent,
    after conditioning on the rate function.
    
    TODO:  maybe there is more information we want to include in this probability,
    like dependence between rates from the same study, or inverse-binomialness in the
    observation -- this is where the _hierarchical_ modelling can come in to play
    """
    vars = []
    for r in asrf.rates():
        pop_vals = pu.population_during(r)
        @mc.observed
        @mc.stochastic(name="rate_%d" % r.id)
        def d_stoc(value=(r.rate_numerator,r.rate_denominator,r.age_start,r.age_end),
                   rate=rate_gp,
                   pop_vals=pop_vals):
            n,d,a0,a1 = value

            # experimental code to reduce influence of large samples
            # (an unprincipled approach to addressing compositional bias)
            # my_rate=float(n)/float(d)
            # d = np.sqrt(d)*10
            # n = my_rate*d
            
            return mc.binomial_like(x=n, n=d,
                                    p=rate_for_range(rate, a0, a1, pop_vals))
        vars.append(d_stoc)

    return vars

def rate_model_vars(asrf):
    """
    
    Generate a list of the PyMC variables which model an age-specific
    rate function that could produce the rates associated with asrf.

    See comments in the code for exact details on the model.
    """
    vars = []
#    import pdb; pdb.set_trace()

    mesh = asrf.fit['age_mesh']

    # logit_rate has an uninformative prior
    #
    # for computational convenience, store values only
    # at mesh points
    logit_rate = mc.Normal('logit(age-specific rate function)',
                           mu    = np.zeros(len(mesh)),
                           tau   = 1.e-12,
                           value = -7.*np.ones(len(mesh)))

    # the rate function is obtained by "non-parametric regression"
    # using a Gaussian process with a nice covariance function to fill
    # in the mesh of logit_rate, and then looking at the inverse logit
    #
    # the interpolation is done in logit space to ensure that
    # the rate is always in the interval [0,1]
    @mc.deterministic
    def age_specific_rate_function(logit_rate=logit_rate):
        M,C = pu.uninformative_prior_gp()
        gp.observe(M,C,mesh,logit_rate)
        return mc.invlogit(M(range(pu.MAX_AGE)))

    vars += [logit_rate, age_specific_rate_function]

    # the following block creates an additional
    # prior on the smoothness of the rate function
    #
    # TODO:  include similar potential code for priors on monotonicity or unimodality
    @mc.potential
    def smoothing_logp(f=logit_rate,tau=1./.5**2):
        return mc.normal_like(np.diff(f), 0.0, tau)
#        return mc.normal_like(np.diff(f,2), 0.0, tau)

    @mc.potential
    def zero_logp(f=age_specific_rate_function, age_start=0, age_end=5):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, 1./.0001**2)
    vars += [smoothing_logp, zero_logp]

    # for each rate on the rate_list, set up an observed stochastic variable
    # which accounts for the probability of the observation given the true
    # rate (and the age-specific population size during the years of the observation)
    #
    # model the rate observations as a binomial random variables, all independent,
    # after conditioning on the rate function.
    #
    # TODO:  maybe there is more information we want to include in this probability,
    # like dependence between rates from the same study, or inverse-binomialness in the
    # observation -- this is where the _hierarchical_ modelling can come in to play
    for r in asrf.rates.all():
        pop_vals = pu.population_during(r)
        @mc.observed
        @mc.stochastic(name="rate_%d" % r.id)
        def d_stoc(value=(r.numerator,r.denominator,r.age_start,r.age_end),
                   rate=age_specific_rate_function,
                   pop_vals=pop_vals):
            n,d,a0,a1 = value
            n = np.minimum(n,d)
            return mc.binomial_like(x=n, n=d,
                                    p=pu.rate_for_range(rate, a0, a1, pop_vals))
        vars.append(d_stoc)

    return vars

