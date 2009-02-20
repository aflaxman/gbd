# model observed rates as binomial draws from a common rate with a beta distribution

from probabilistic_utils import *

MIN_BETA_CONFIDENCE = .1

def model_vars(asrf):
    vars = {}
    
    age_mesh = asrf.fit['age_mesh']
    out_age_mesh = asrf.fit['out_age_mesh']
    initial_rate = np.array(asrf.fit['normal_approx'])
    ones_mesh = np.ones(len(age_mesh))

    logit_rate = mc.Normal('logit_rate', mu=0.*ones_mesh, tau=1./(1.e3)**2)
    logit_rate.value = mc.logit(.5*.0001 + .5*(1. - .0001) * np.array(initial_rate)[age_mesh])
    
    log_confidence = mc.Normal('log_confidence', mu=12.*ones_mesh, tau=1./(1.e0)**2)
    log_confidence.value = 10*np.ones(len(age_mesh))
    

    @mc.deterministic
    def rate(logit_rate=logit_rate):
        return np.maximum(NEARLY_ZERO, mc.invlogit(gp_interpolate(age_mesh, logit_rate, out_age_mesh)))

    @mc.deterministic
    def confidence(log_confidence=log_confidence):
        return MIN_BETA_CONFIDENCE + np.exp(gp_interpolate(age_mesh, log_confidence, out_age_mesh))

    @mc.deterministic
    def alpha(rate=rate, confidence=confidence):
        return rate * confidence
    @mc.deterministic
    def beta(rate=rate, confidence=confidence):
        return (1. - rate) * confidence

    vars['rate'] = rate
    vars['alpha'], vars['beta'] = alpha, beta
    vars['rate_related'] = [logit_rate, log_confidence, confidence]


    tau_smooth_rate = mc.InverseGamma('smoothing_tau_%d'%asrf.id, 0.1, 1.0, value=1.)
    vars['hyper params'] = [tau_smooth_rate]

    @mc.potential
    def smooth_logit_rate(f=logit_rate, tau=tau_smooth_rate):
        return mc.normal_like(np.diff(f), 0.0, tau)

    @mc.potential
    def smooth_confidence(f=log_confidence, tau=1./(.5)**2):
        return mc.normal_like(np.diff(f), 0.0, tau)
    vars['smooth prior'] = smooth_logit_rate, smooth_confidence

    @mc.potential
    def initially_zero(f=rate, age_start=0, age_end=5, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
    vars['initially zero prior'] = initially_zero

#    @mc.potential
#    def finally_zero(f=rate, age_start=90, age_end=100, tau=1./(1e-4)**2):
#        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
#    vars['finally zero prior'] = finally_zero


    vars['logit_p'] = []
    vars['p'] = []
    vars['beta_potential'] = []
    vars['observed_rates'] = []
    for r in asrf.rates.all():
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
        vars['logit_p'] += [logit_p]
        vars['p'] += [p]
        vars['beta_potential'] += [potential_p]
        vars['observed_rates'] += [obs]

    return vars

def save_map(vars, asrf):
    #import pdb; pdb.set_trace()
    asrf.fit['map'] = list(vars['rate'].value)
    asrf.save()

def save_mcmc(vars, asrf):
    trace_len = len(vars['rate'].trace())
    #rate = np.array(vars['rate'].trace())
    rate = np.array([mc.rbeta(alpha, beta) for alpha, beta in zip(vars['alpha'].trace(), vars['beta'].trace())])
    sr = []
    for ii in asrf.fit['out_age_mesh']:
        sr.append(sorted(rate[:,ii]))
    asrf.fit['mcmc_lower_cl'] = [sr[ii][int(.025*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_median'] = [sr[ii][int(.5*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_upper_cl'] = [sr[ii][int(.975*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_mean'] = list(np.mean(rate, 0))

    asrf.save()
