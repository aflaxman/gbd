# model observed rates as binomial draws from a common rate with a beta distribution

from probabilistic_utils import *

MIN_PARAM_VAL = .001
N_SAMPLES = 25

def model_vars(asrf):
    vars = {}
    
    age_mesh = asrf.fit['age_mesh']
    out_age_mesh = asrf.fit['out_age_mesh']
    initial_rate = np.array(asrf.fit['normal_approx'])
    ones_mesh = np.ones(len(age_mesh))

    log_alpha = mc.Normal('log_alpha', mu=0.*ones_mesh, tau=1./(1.e2)**2)
    log_alpha.value = np.log(1.+100.*initial_rate[age_mesh])

    @mc.deterministic
    def alpha(log_rate=log_alpha):
        return MIN_PARAM_VAL + np.exp(gp_interpolate(age_mesh, log_rate, out_age_mesh))
    vars['log(alpha)'], vars['alpha'] = log_alpha, alpha
    
    log_beta = mc.Normal('log_beta', mu=0.*ones_mesh, tau=1./(1.e2)**2)
    log_beta.value = np.log(1.+100.*(1.-initial_rate[age_mesh]))

    @mc.deterministic
    def beta(log_rate=log_beta):
        return MIN_PARAM_VAL + np.exp(gp_interpolate(age_mesh, log_rate, out_age_mesh))
    vars['log(beta)'], vars['beta'] = log_beta, beta

    @mc.deterministic
    def Ep(alpha=alpha, beta=beta):
        return alpha / (alpha + beta)
    vars['E[p]'] = Ep


    @mc.potential
    def smooth(f=log_alpha, g=log_beta, tau=1./(.5)**2):
        return mc.normal_like(np.diff(f), 0.0, tau) + mc.normal_like(np.diff(g), 0.0, tau)
    vars['smooth prior'] = smooth

    @mc.potential
    def initially_zero(f=Ep, age_start=0, age_end=12, tau=1./(1e-5)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
    vars['initially zero prior'] = initially_zero

    @mc.potential
    def finally_zero(f=Ep, age_start=90, age_end=100, tau=1./(1e-4)**2):
        return mc.normal_like(f[range(age_start, age_end)], 0.0, tau)
    vars['finally zero prior'] = finally_zero


    vars['observed_rates'] = []
    for r in asrf.rates.all():
        @mc.observed
        @mc.stochastic(name="rate_%d" % r.id)
        def obs(value=(r.numerator,r.denominator,r.age_start,r.age_end),
                alpha=alpha, beta=beta,
                pop_vals=r.population()):
            numerator, denominator, a0, a1 = value
            p_samp = mc.rbeta(rate_for_range(alpha, a0, a1, pop_vals),
                              rate_for_range(beta, a0, a1, pop_vals),
                              N_SAMPLES)
            return mc.binomial_like([numerator]*N_SAMPLES, [denominator]*N_SAMPLES, p_samp) / N_SAMPLES
        vars['observed_rates'] += [obs]

    return vars

def save_map(vars, asrf):
    asrf.fit['map'] = list(vars['E[p]'].value)
    asrf.save()

def save_mcmc(vars, asrf):
    trace_len = len(vars['alpha'].trace())
    rate = np.array([mc.rbeta(alpha, beta) for alpha, beta in zip(vars['alpha'].trace(), vars['beta'].trace())])
    sr = []
    for ii in asrf.fit['out_age_mesh']:
        sr.append(sorted(rate[:,ii]))
    asrf.fit['mcmc_lower_cl'] = [sr[ii][int(.025*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_median'] = [sr[ii][int(.5*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_upper_cl'] = [sr[ii][int(.975*trace_len)] for ii in asrf.fit['out_age_mesh']]
    asrf.fit['mcmc_mean'] = list(np.mean(rate, 0))

    asrf.save()
