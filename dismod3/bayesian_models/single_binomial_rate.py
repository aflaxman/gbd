# model all observed rates as binomial draws from a single rate function

from probabilistic_utils import *

def setup_rate_model(rf, rate_stoch=None):
    """
    create pymc stochastic variables to model observed rates as
    realizations of binomial random variables with a single underlying
    age-specific rate function, p(a).

    the vars are stored in the rf.vars dict, and the rate stoch is
    also stored in rf.rate_stoch for easy access.

    if rate_stoch is a pymc stochastic with length corresponding to
    rf.fit['out_age_mesh'], it will be used as the rate stoch, which
    allows disease models to use rate models as building blocks
    """
    rf.vars = {}

    if rate_stoch:
        rf.rate_stoch = rate_stoch
        rf.vars['rf_%d'%rf.id] = rate_stoch
    else:
        add_stochs(rf, 'rf_%d'%rf.id,
                   initial_value=trim(rf.fit['normal_approx'], NEARLY_ZERO, 1. - NEARLY_ZERO),
                   transform='logit')
        rf.rate_stoch = rf.vars['rf_%d'%rf.id]

    tau_smooth = mc.InverseGamma('smoothing_tau_%d'%rf.id, 0.1, 1.0, value=1.)
    rf.vars['hyperparams'] = [tau_smooth]

    @mc.potential
    def smooth(f=rf.rate_stoch, tau=tau_smooth):
        return mc.normal_like(np.diff(f), 0.0, tau)

    @mc.potential
    def initially_zero(f=rf.rate_stoch, age_start=0, age_end=5, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
    rf.vars['priors'] = [smooth, initially_zero]

#    @mc.potential
#    def finally_zero(f=rf.rate_stoch, age_start=90, age_end=101, tau=1./(1e-4)**2):
#        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
#    rf.vars['priors'] += [finally_zero]

    rf.vars['observed_rates'] = []
    for r in rf.rates.all():
        r.numerator = min(r.numerator,r.denominator)
        @mc.observed
        @mc.stochastic(name="rate_%d^%d" % (r.id, rf.id))
        def d_stoc(value=(r.numerator,r.denominator,r.age_start,r.age_end),
                   rate=rf.rate_stoch,
                   pop_vals=r.population()):
            n,d,a0,a1 = value
            logp = mc.binomial_like(x=n, n=d,
                                    p=rate_for_range(rate, a0, a1, pop_vals))
            if np.isnan(logp):
                import pdb; pdb.set_trace()
            return logp
        rf.vars['observed_rates'].append(d_stoc)

