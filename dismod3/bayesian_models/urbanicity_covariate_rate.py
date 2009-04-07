# model observed rates as binomial draws from a common rate with a beta distribution
# that includes a covariate for urban/rural (and estimates a value which is rural)

from probabilistic_utils import *

MIN_CONFIDENCE = 1
MAX_CONFIDENCE = 100000

def setup_rate_model(rf, rate_stoch=None):
    rf.vars = {}

    #############################################################################
    # set up the age-specific Beta stochastic variables
    #
         
    if rate_stoch:
        Erf = rate_stoch
        rf.vars['Erf_%d'%rf.id] = Erf
    else:
        add_stoch_to_rf_vars(rf, 'Erf_%d'%rf.id,
                             initial_value=trim(rf.fit['normal_approx'], NEARLY_ZERO, 1. - NEARLY_ZERO),
                             transform='logit')
        Erf = rf.vars['Erf_%d'%rf.id]

    confidence = mc.Normal('conf_%d'%rf.id, mu=1000.0, tau=1./(300.)**2)
    rf.vars['confidence'] = confidence
    
    @mc.deterministic(name='alpha_%d'%rf.id)
    def alpha(rate=Erf, confidence=confidence):
        return rate * trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    @mc.deterministic(name='beta_%d'%rf.id)
    def beta(rate=Erf, confidence=confidence):
        return (1. - rate) * trim(confidence, MIN_CONFIDENCE, MAX_CONFIDENCE)

    rf.vars['alpha'], rf.vars['beta'] = alpha, beta

    ##########################################################################
    # set up the stochs that will be saved and model the interactions
    # between this and the rest of a disease model
    
    rf.map_fit_stoch = Erf
    rf.mcmc_fit_stoch = Erf
    rf.rate_stoch = Erf

    ###################################################
    # set up stoch for urban covariate adjustment:
    #  p_i ~ Beta(alpha, beta) + gamma * urban_i
    gamma = mc.Normal('gamma', 0., 1./10.**2)
    rf.vars['gamma'] = gamma

    ########################################################################
    # set up stochs for the priors and observed data
    #
    
    add_priors_to_rf_vars(rf)

    rf.vars['beta_binom_stochs'] = []
    rf.vars['observed_rates'] = []
    for r in rf.rates.all():
        # ensure all rate data is valid
        # TODO: raise exceptions to have users fix an errors
        r.numerator = min(r.numerator, r.denominator)

        logit_p = mc.Normal('logit_p_%d' % r.id, 0., 1/(10.)**6,
                            value=mc.logit((1. + r.numerator)/(2. + r.denominator)),
                            verbose=0)

        @mc.deterministic(name='p_%d' % r.id)
        def p(logit_p=logit_p, gamma=gamma, urban=(r.params['Urbanicity'] == 'Urban')):
            return trim(mc.invlogit(logit_p + gamma * urban), NEARLY_ZERO, 1. - NEARLY_ZERO)

        @mc.potential(name='beta_potential_%d' % r.id)
        def potential_p(logit_p=logit_p,
                        alpha=alpha, beta=beta,
                        a0=r.age_start, a1=r.age_end,
                        pop_vals=r.population()):
            a = rate_for_range(alpha, a0, a1, pop_vals)
            b = rate_for_range(beta, a0, a1, pop_vals)
            return mc.beta_like(mc.invlogit(logit_p), a, b)

        obs = mc.Binomial("rate_%d" % r.id, value=r.numerator, n=r.denominator, p=p, observed=True)
        rf.vars['beta_binom_stochs'] += [logit_p]
        rf.vars['observed_rates'] += [p, potential_p, obs]

