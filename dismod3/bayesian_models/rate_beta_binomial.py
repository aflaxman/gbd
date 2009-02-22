# model observed rates as binomial draws from a common rate with a beta distribution

from probabilistic_utils import *

MIN_CONFIDENCE = 1.e2
MAX_CONFIDENCE = 1.e10

def setup_rate_model(rf, rate_stoch=None):
    rf.vars = {}

    initial_value=trim(rf.fit['normal_approx'], NEARLY_ZERO, 1. - NEARLY_ZERO)
    add_stoch_to_rf_vars(rf, 'Erf_%d'%rf.id,
                         initial_value,
                         transform='logit')
    Erf = rf.vars['Erf_%d'%rf.id]
    add_stoch_to_rf_vars(rf, 'confidence_%d'%rf.id,
                         initial_value=MIN_CONFIDENCE*np.ones(len(rf.fit['normal_approx'])),
                         transform='log')
    confidence = rf.vars['confidence_%d'%rf.id]

    @mc.deterministic
    def alpha(rate=Erf, confidence=confidence):
        return rate * trim(confidence + MIN_CONFIDENCE, 0, MAX_CONFIDENCE)
    @mc.deterministic
    def beta(rate=Erf, confidence=confidence):
        return (1. - rate) * trim(confidence + MIN_CONFIDENCE, 0, MAX_CONFIDENCE)
    rf.vars['alpha'], rf.vars['beta'] = alpha, beta

    if rate_stoch:
        rf.rate_stoch = rate_stoch
        @mc.potential(name='rate_link_%d'%rf.id)
        def rate_link(alpha=alpha, beta=beta, rate=rate_stoch):
            return mc.beta_like(rate, alpha, beta)
        rf.vars['rate_link'] = rate_link
    else:
        @mc.deterministic(name='rate_%d')
        def rate_stoch(alpha=alpha, beta=beta):
            return mc.rbeta(alpha, beta)
        rf.rate_stoch = rate_stoch
        rf.vars['rate'] = rf.rate_stoch

    rf.map_fit_stoch = Erf
    rf.mcmc_fit_stoch = rf.rate_stoch
    
    tau_smooth_rate = mc.InverseGamma('smooth_rate_tau_%d'%rf.id, 1., 4., value=4.)
    tau_smooth_confidence = mc.InverseGamma('smooth_conf_tau_%d'%rf.id, 1., 4., value=4.)
    tau_smooth_rate = 1./.5**2
    tau_smooth_confidence = 1./.5**2
    rf.vars['hyper-params'] = [tau_smooth_rate, tau_smooth_confidence]

    @mc.potential
    def smooth_rate(f=rf.vars['logit(Erf_%d)'%rf.id], tau=tau_smooth_rate):
        return mc.normal_like(np.diff(f), 0.0, tau)
    rf.vars['smooth prior'] = [smooth_rate]
    @mc.potential
    def smooth_confidence(f=rf.vars['log(confidence_%d)'%rf.id], tau=tau_smooth_confidence):
        return mc.normal_like(np.diff(f), 0.0, tau)
    rf.vars['smooth prior'] += [smooth_confidence]

#    @mc.potential
#    def initially_zero(f=Erf, age_start=0, age_end=5, tau=1./(1e-4)**2):
#        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
#    rf.vars['initially zero prior'] = initially_zero

#    @mc.potential
#    def finally_zero(f=rf.rate_stoch, age_start=90, age_end=100, tau=1./(1e-4)**2):
#        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
#    rf.vars['finally zero prior'] = finally_zero


    rf.vars['observed_rates'] = []
    for r in rf.rates.all():
        r.numerator = min(r.numerator, r.denominator)

        logit_p = mc.Normal('logit_p_%d' % r.id, 0., 1/(10.)**6,
                            value=mc.logit((1. + r.numerator)/(2. + r.denominator)))

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

        rf.vars['observed_rates'] += [logit_p, p, potential_p, obs]

