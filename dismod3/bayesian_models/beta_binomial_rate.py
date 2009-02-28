# model observed rates as binomial draws from a common rate with a beta distribution

from probabilistic_utils import *

MIN_CONFIDENCE = 1
MAX_CONFIDENCE = 100000

def setup_rate_model(rf, rate_stoch=None):
    rf.vars = {}

    #############################################################################
    # set up the age-specific Beta stochastic variables
    #

    initial_value=trim(rf.fit['normal_approx'], NEARLY_ZERO, 1. - NEARLY_ZERO)
    add_stoch_to_rf_vars(rf, 'Erf_%d'%rf.id,
                         initial_value,
                         transform='logit')
    Erf = rf.vars['Erf_%d'%rf.id]

    confidence = mc.Normal('conf_%d'%rf.id, mu=0.0, tau=1./(5.)**2)
    rf.vars['confidence'] = confidence
    
    @mc.deterministic(name='alpha_%d'%rf.id)
    def alpha(rate=Erf, confidence=confidence):
        return rate * (MIN_CONFIDENCE + mc.invlogit(confidence)*MAX_CONFIDENCE)

    @mc.deterministic(name='beta_%d'%rf.id)
    def beta(rate=Erf, confidence=confidence):
        return (1. - rate) * (MIN_CONFIDENCE + mc.invlogit(confidence)*MAX_CONFIDENCE)

    rf.vars['alpha'], rf.vars['beta'] = alpha, beta

    ##########################################################################
    # set up the stochs that will be saved and model the interactions
    # between this and the rest of a disease model
    
    rf.map_fit_stoch = Erf
    rf.mcmc_fit_stoch = Erf
    rf.rate_stoch = Erf
         
    if rate_stoch:
        rf.rate_stoch = rate_stoch
        @mc.potential(name='rate_link_%d'%rf.id)
        def rate_link(alpha=alpha, beta=beta, rate=rate_stoch):
            return mc.beta_like(rate, alpha, beta)
        rf.vars['rate_link'] = rate_link
        rf.vars['linked-in rate'] = rf.rate_stoch
    else:
        @mc.deterministic(name='realized_rate_%d'%rf.id)
        def realized_rate(alpha=alpha, beta=beta):
            return mc.rbeta(alpha + NEARLY_ZERO, beta + NEARLY_ZERO)
        rf.vars['realized rate'] = realized_rate
        rf.mcmc_fit_stoch = realized_rate


    ########################################################################
    # set up stochs for the priors and observed data
    #
    
    add_priors_to_rf_vars(rf)

    rf.vars['beta_binom_stochs'] = []
    rf.vars['observed_rates'] = []
    for r in rf.rates.all():
        r.numerator = min(r.numerator, r.denominator)

        logit_p = mc.Normal('logit_p_%d' % r.id, 0., 1/(10.)**6,
                            value=mc.logit((1. + r.numerator)/(2. + r.denominator)),
                            verbose=0)

        @mc.deterministic(name='p_%d' % r.id)
        def p(logit_p=logit_p):
            return mc.invlogit(logit_p)

        @mc.potential(name='beta_potential_%d' % r.id)
        def potential_p(p=p,
                        alpha=alpha, beta=beta,
                        a0=r.age_start, a1=r.age_end,
                        pop_vals=r.population()):
            a = rate_for_range(alpha, a0, a1, pop_vals)
            b = rate_for_range(beta, a0, a1, pop_vals)
            return mc.beta_like(p, a, b)

        obs = mc.Binomial("rate_%d" % r.id, value=r.numerator, n=r.denominator, p=p, observed=True)
        rf.vars['beta_binom_stochs'] += [logit_p]
        rf.vars['observed_rates'] += [p, potential_p, obs]

