# model observed rates as binomial draws from a common rate with a beta distribution

from probabilistic_utils import *

MIN_BETA_CONFIDENCE = .1

def setup_rate_model(rf, rate_stoch=None):
    rf.vars = {}

    if rate_stoch:
        rf.rate_stoch = rate_stoch
        rf.vars['rf_%d'%rf.id] = rate_stoch
    else:
        add_stoch_to_rf_vars(rf, 'rf_%d'%rf.id,
                             initial_value=trim(rf.fit['normal_approx'], NEARLY_ZERO, 1. - NEARLY_ZERO),
                             transform='logit')
        rf.rate_stoch = rf.vars['rf_%d'%rf.id]

    add_stoch_to_rf_vars(rf, 'confidence_%d'%rf.id,
                         initial_value=10*np.ones(len(rf.fit['normal_approx'])),
                         transform='log')
    confidence = rf.vars['confidence_%d'%rf.id]

    @mc.deterministic
    def alpha(rate=rf.rate_stoch, confidence=confidence):
        return rate * confidence
    @mc.deterministic
    def beta(rate=rf.rate_stoch, confidence=confidence):
        return (1. - rate) * confidence
    rf.vars['alpha'], rf.vars['beta'] = alpha, beta

    tau_smooth_rate = mc.InverseGamma('smoothing_rate_tau_%d'%rf.id, 0.1, 1.0, value=1.)
    tau_smooth_confidence = mc.InverseGamma('smoothing_confidence_tau_%d'%rf.id, 0.1, 1.0, value=1.)
    rf.vars['hyper-params'] = [tau_smooth_rate, tau_smooth_confidence]

    @mc.potential
    def smooth_rate(f=rf.rate_stoch, tau=tau_smooth_rate):
        return mc.normal_like(np.diff(f), 0.0, tau)

    @mc.potential
    def smooth_confidence(f=confidence, tau=tau_smooth_confidence):
        return mc.normal_like(np.diff(f), 0.0, tau)
    rf.vars['smooth prior'] = smooth_rate, smooth_confidence

    @mc.potential
    def initially_zero(f=rf.rate_stoch, age_start=0, age_end=5, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
    rf.vars['initially zero prior'] = initially_zero

#    @mc.potential
#    def finally_zero(f=rf.rate_stoch, age_start=90, age_end=100, tau=1./(1e-4)**2):
#        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
#    rf.vars['finally zero prior'] = finally_zero


    rf.vars['logit_p'] = []
    rf.vars['p'] = []
    rf.vars['beta_potential'] = []
    rf.vars['observed_rates'] = []
    for r in rf.rates.all():
        r.numerator = min(r.numerator, r.denominator)
        logit_p = mc.Normal('logit_p_%d' % r.id, 0., 1/(10.)**2, value=mc.logit((1. + r.numerator)/(2. + r.denominator)))
        @mc.deterministic(name='p_%d' % r.id)
        def p(logit_p=logit_p):
            return mc.invlogit(logit_p)

        @mc.potential(name='beta_potential_%d' % r.id)
        def potential_p(p=p,
            alpha=alpha, beta=beta, pop_vals=r.population(), a0=r.age_start, a1=r.age_end):
            a,b = rate_for_range(alpha, a0, a1, pop_vals), rate_for_range(beta, a0, a1, pop_vals)
            logp = mc.beta_like(p, a, b)
            if np.isnan(logp) or np.isinf(logp):
                import pdb; pdb.set_trace()
            return logp

        obs = mc.Binomial("rate_%d" % r.id, value=r.numerator, n=r.denominator, p=p, observed=True)
        rf.vars['logit_p'] += [logit_p]
        rf.vars['p'] += [p]
        rf.vars['beta_potential'] += [potential_p]
        rf.vars['observed_rates'] += [obs]

    return vars
