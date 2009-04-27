import numpy as np
import pymc as mc

from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import NEARLY_ZERO
from model_utils import generate_prior_potentials

MIN_CONFIDENCE = 1
MAX_CONFIDENCE = 100000

def setup(dm, data_type):
    #############################################################################
    # set up age-specific rate function
    #
    
    initial_value = dm.get_initial_value(data_type)
    param_mesh = dm.get_param_age_mesh()
    est_mesh = dm.get_estimate_age_mesh()
    rate_str = data_type.replace('data','')

    logit_rate = mc.Normal('logit(%s)' % rate_str,
                           mu=np.zeros(len(param_mesh)),
                           tau=1.e-2,
                           value=mc.logit(initial_value[param_mesh]),
                           verbose=0)
    @mc.deterministic(name=rate_str)
    def rate(logit_rate=logit_rate):
        return probabilistic_utils.interpolate(param_mesh, mc.invlogit(logit_rate), est_mesh)
    
    confidence = mc.Normal('conf_%s' % rate_str, mu=1000.0, tau=1./(300.)**2)
    
    @mc.deterministic(name='alpha_%s' % rate_str)
    def alpha(rate=rate, confidence=confidence):
        return rate * probabilistic_utils.trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    @mc.deterministic(name='beta_%s' % rate_str)
    def beta(rate=rate, confidence=confidence):
        return (1. - rate) * probabilistic_utils.trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    ########################################################################
    # set up priors and observed data
    #
    prior_str = dm.get_priors(data_type)
    priors = generate_prior_potentials(prior_str, rate, confidence)

    logit_p_stochs = []
    p_stochs = []
    beta_potentials = []
    observed_rates = []
    for d in dm.data:
        if d['data_type'] != data_type:
            continue
        id = d['id']
        
        # ensure all rate data is valid
        # TODO: raise exceptions to have users fix an errors
        d.update(value=probabilistic_utils.trim(d['value'], NEARLY_ZERO, 1.-NEARLY_ZERO),
                 standard_error=max(d['standard_error'], 0.0001))

        logit_p = mc.Normal('logit(p_%d)' % id, 0., 1/(10.)**2,
                            value=mc.logit(d['value'] + NEARLY_ZERO),
                            verbose=0)

        p = mc.InvLogit('p_%d' % id, logit_p)

        @mc.potential(name='beta_potential_%d' % id)
        def potential_p(p=p,
                        alpha=alpha, beta=beta,
                        a0=d['age_start'], a1=d['age_end'],
                        age_weights=d['age_weights']):
            a = probabilistic_utils.rate_for_range(alpha, a0, a1, age_weights)
            b = probabilistic_utils.rate_for_range(beta, a0, a1, age_weights)
            return mc.beta_like(probabilistic_utils.trim(p, NEARLY_ZERO, 1. - NEARLY_ZERO), a, b)

        denominator = max(100., d['value'] * (1 - d['value']) / d['standard_error']**2.)
        numerator = d['value'] * denominator
        obs = mc.Binomial("data_%d" % id, value=numerator, n=denominator, p=p, observed=True)

        logit_p_stochs.append(logit_p)
        p_stochs.append(p)
        beta_potentials.append(potential_p)
        observed_rates.append(obs)
        
    return mc.Model(locals())
