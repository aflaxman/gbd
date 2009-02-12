# model all observed rates as binomial draws from a single rate function

from probabilistic_utils import *

def model_vars(asrf):
    vars = {}

    add_rate_stochs(vars, 'asrf_%d'%asrf.id, asrf.fit['age_mesh'], asrf.fit['out_age_mesh'])
    vars['logit(asrf_%d)'%asrf.id].value = NEARLY_ZERO + (1. - NEARLY_ZERO) * mc.logit(np.array(asrf.fit['normal_approx'])[asrf.fit['age_mesh']])

    tau_smooth = mc.InverseGamma('smoothing_tau_%d'%asrf.id, 0.1, 1.0, value=1.)
    vars['hyper params'] = [tau_smooth]

    #import pdb; pdb.set_trace()
    @mc.potential
    def smooth(f=vars['logit(asrf_%d)'%asrf.id], tau=tau_smooth):
        return mc.normal_like(np.diff(f), 0.0, tau)

    @mc.potential
    def initially_zero(f=vars['asrf_%d'%asrf.id], age_start=0, age_end=5, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)

    @mc.potential
    def finally_zero(f=vars['asrf_%d'%asrf.id], age_start=90, age_end=100, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)

    vars['priors'] = [smooth, initially_zero, finally_zero]

    vars['observed_rates'] = []
    for r in asrf.rates.all():
        @mc.observed
        @mc.stochastic(name="rate_%d" % r.id)
        def d_stoc(value=(r.numerator,r.denominator,r.age_start,r.age_end),
                   rate=vars['asrf'],
                   pop_vals=r.population()):
            n,d,a0,a1 = value
            n = max(n,d)
            return mc.binomial_like(x=n, n=d,
                                    p=rate_for_range(rate, a0, a1, pop_vals))
        vars['observed_rates'].append(d_stoc)

    return vars


def save_map(vars, asrf):
    asrf.fit['map'] = list(vars['asrf_%d'%asrf.id].value)
    asrf.save()

def save_mcmc(vars, asrf):
    rate = vars['asrf_%d'%asrf.id].trace()
    trace_len = len(rate)
    
    sr = []
    for ii in asrf.fit['out_age_mesh']:
        sr.append(sorted(rate[:,ii]))
    asrf.fit['mcmc_lower_cl'] = [sr[ii][int(.025*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_median'] = [sr[ii][int(.5*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_upper_cl'] = [sr[ii][int(.975*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_mean'] = list(np.mean(rate, 0))

    asrf.save()
