import numpy as np
import pymc as mc

from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import NEARLY_ZERO
from model_utils import *

MIN_CONFIDENCE = 1
MAX_CONFIDENCE = 100000

def initialize(dm, data_type='prevalence data'):
    dm.set_units(data_type, '(per person-year)')
    dm.fit_normal_approx(data_type)

def map_fit(dm, vars, data_type='prevalence data'):
    map = mc.MAP(vars)
    map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
    dm.set_map(data_type, vars['rate_stoch'].value)
    return map

def mcmc_fit(dm, vars, data_type='prevalence data'):
    mcmc = mc.MCMC(vars)
    mcmc.sample(iter=4000, burn=1000, thin=2, verbose=1)
    store_mcmc_fit(dm, vars['rate_stoch'], data_type)

def store_mcmc_fit(dm, rate_stoch, data_type):
    rate = rate_stoch.trace()
    trace_len = len(rate)
    age_len = len(dm.get_estimate_age_mesh())
    
    sr = []
    for ii in xrange(age_len):
        sr.append(sorted(rate[:,ii]))
    dm.set_mcmc('lower_ui', data_type, [sr[ii][int(.025*trace_len)] for ii in xrange(age_len)])
    dm.set_mcmc('median', data_type, [sr[ii][int(.5*trace_len)] for ii in xrange(age_len)])
    dm.set_mcmc('upper_ui', data_type, [sr[ii][int(.975*trace_len)] for ii in xrange(age_len)])
    dm.set_mcmc('mean', data_type, np.mean(rate, 0))

def setup(dm, data_type='prevalence data', rate_stoch=None):
    """
    Generate the PyMC variables for a beta binomial model of
    a single rate function, and return it as a dict
    """
    vars = {}
    rate_str = data_type.replace('data','')
    est_mesh = dm.get_estimate_age_mesh()
    if np.any(np.diff(est_mesh) != 1):
        raise ValueError, 'ERROR: Gaps in estimation age mesh must all equal 1'

    #############################################################################
    # set up age-specific rate function, if it does not yet exist
    #
    if not rate_stoch:
        initial_value = dm.get_initial_value(data_type)
        param_mesh = dm.get_param_age_mesh()

        # find the logit of the initial values, which is a little bit
        # of work because initial values are sampled from the est_mesh,
        # but the logit_initial_values are needed on the param_mesh
        logit_initial_value = mc.invlogit(
            probabilistic_utils.interpolate(est_mesh, initial_value, param_mesh))
        
        logit_rate = mc.Normal('logit(%s)' % rate_str,
                               mu=np.zeros(len(param_mesh)),
                               tau=1.e-2,
                               value=logit_initial_value,
                               verbose=0)
        vars['logit_rate'] = logit_rate

        @mc.deterministic(name=rate_str)
        def rate_stoch(logit_rate=logit_rate):
            return probabilistic_utils.interpolate(param_mesh, mc.invlogit(logit_rate), est_mesh)

    vars['rate_stoch'] = rate_stoch

    confidence = mc.Normal('conf_%s' % rate_str, mu=1000.0, tau=1./(3.)**2)
    
    @mc.deterministic(name='alpha_%s' % rate_str)
    def alpha(rate=rate_stoch, confidence=confidence):
        return rate * probabilistic_utils.trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    @mc.deterministic(name='beta_%s' % rate_str)
    def beta(rate=rate_stoch, confidence=confidence):
        return (1. - rate) * probabilistic_utils.trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    vars['additional rate params'] = [confidence, alpha, beta]

    ########################################################################
    # set up priors and observed data
    #
    prior_str = dm.get_priors(data_type)
    print 'setting up priors from:\n%s' % prior_str
    vars['priors'] = generate_prior_potentials(prior_str, est_mesh, rate_stoch, confidence)

    vars['logit_p_stochs'] = []
    vars['p_stochs'] = []
    vars['beta_potentials'] = []
    vars['observed_rates'] = []
    for d in dm.data:
        # set up observed stochs for all relevant data
        id = d['id']
        
        if d['data_type'] != data_type:
            continue
        if d['value'] == MISSING:
            print 'WARNING: data %d missing value' % id
            continue

        # ensure all rate data is valid
        d_val = dm.value_per_1(d)
        d_se = dm.se_per_1(d)
        
        if d_val < 0 or d_val > 1:
            print 'WARNING: data %d not in range [0,1]' % id
            continue

        if d['age_start'] < est_mesh[0] or d['age_end'] > est_mesh[-1]:
            raise ValueError, 'Data %d is outside of estimation range---([%d, %d] is not inside [%d, %d])' \
                % (d['id'], d['age_start'], d['age_end'], est_mesh[0], est_mesh[-1])
        age_indices = indices_for_range(est_mesh, d['age_start'], d['age_end'])
        age_weights = d['age_weights']
        # if the data has a standard error, model it as a realization
        # of a beta binomial r.v.
        if d_se > 0:
            logit_p = mc.Normal('logit(p_%d)' % id, 0., 1/(10.)**2,
                                value=mc.logit(d_val + NEARLY_ZERO),
                                verbose=0)
            p = mc.InvLogit('p_%d' % id, logit_p)

            @mc.potential(name='beta_potential_%d' % id)
            def potential_p(p=p,
                            alpha=alpha, beta=beta,
                            age_indices=age_indices,
                            age_weights=age_weights):
                a = probabilistic_utils.rate_for_range(alpha, age_indices, age_weights)
                b = probabilistic_utils.rate_for_range(beta, age_indices, age_weights)
                return mc.beta_like(probabilistic_utils.trim(p, NEARLY_ZERO, 1. - NEARLY_ZERO), a, b)

            denominator = max(100., d['value'] * (1 - d['value']) / d['standard_error']**2.)
            numerator = d['value'] * denominator
            obs = mc.Binomial('data_%d' % id, value=numerator, n=denominator, p=p, observed=True)

            vars['logit_p_stochs'].append(logit_p)
            vars['p_stochs'].append(p)
            vars['beta_potentials'].append(potential_p)
        else:
            # if the data is a point estimate with no uncertainty
            # recorded, model it as a realization of a beta
            @mc.observed
            @mc.stochastic(name='data_%d' % id)
            def obs(value=d_val,
                    alpha=alpha, beta=beta,
                    age_indices=age_indices,
                    age_weights=age_weights):
                a = probabilistic_utils.rate_for_range(alpha, age_indices, age_weights)
                b = probabilistic_utils.rate_for_range(beta, age_indices, age_weights)
                return mc.beta_like(probabilistic_utils.trim(value, NEARLY_ZERO, 1. - NEARLY_ZERO), a, b)
            
        vars['observed_rates'].append(obs)
        
    return vars
