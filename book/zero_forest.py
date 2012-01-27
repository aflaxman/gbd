""" Fit data with several rate models and generate forest plot"""


import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}

### @export 'zero-inflated-sim-study'
pi_true = .0001
delta_true = 5.
n_pred = 1.e9

iter = 20000
burn = 10000
thin = 10


### @export 'data'
n = pl.array(pl.exp(mc.rnormal(11, 1**-2, size=32)), dtype=int)
k = pl.array(mc.rnegative_binomial(n*pi_true, delta_true), dtype=float)
k[:4] = 0. # zero-inflated model
r = k/n
s = pl.sqrt(r * (1-r) / n)

n_min = min(n)
n_max = max(n)


### @export 'zibb-model'
alpha = mc.Uninformative('alpha', value=4.)
beta = mc.Uninformative('beta', value=1000.)
pi_mean = mc.Lambda('pi_mean', lambda alpha=alpha, beta=beta: alpha/(alpha+beta))
pi = mc.Beta('pi', alpha, beta, value=pl.maximum(1.e-12, pl.minimum(1-1.e-12, r)))
phi = mc.Uniform('phi', lower=0., upper=1., value=.01)

nonzeros = r != 0.
num_nonzeros = nonzeros.sum()
@mc.potential
def obs(pi=pi, phi=phi):
    logp = pl.log(1-phi)*num_nonzeros + mc.binomial_like(r[nonzeros]*n[nonzeros], n[nonzeros], pi[nonzeros])
    for n_i in n[~nonzeros]:
        logp += pl.log(phi + (1-phi) * pl.exp(pl.log(1-pi[~nonzeros]) * n[~nonzeros])).sum()
    return logp

@mc.deterministic
def pred(alpha=alpha, beta=beta, phi=phi):
    if pl.rand() < phi:
        return 0
    else:
        return mc.rbinomial(n_pred, mc.rbeta(alpha, beta)) / float(n_pred)

mcmc = mc.MCMC([alpha, beta, pi_mean, pi, phi, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta, phi])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)
mcmc.sample(iter*10, burn*10, thin*10)

results['Zero-inflated beta binomial'] = dict(pi=pi_mean.stats(), pred=pred.stats(), phi=phi.stats())


### @export 'beta-binomial-model'
alpha = mc.Uninformative('alpha', value=4.)
beta = mc.Uninformative('beta', value=1000.)
pi_mean = mc.Lambda('pi_mean', lambda alpha=alpha, beta=beta: alpha/(alpha+beta))
pi = mc.Beta('pi', alpha, beta, value=pl.maximum(1.e-12, pl.minimum(1-1.e-12, r)))

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(r*n, n, pi)

@mc.deterministic
def pred(alpha=alpha, beta=beta):
    return mc.rbinomial(n_pred, mc.rbeta(alpha, beta)) / float(n_pred)

### @export 'beta-binomial-fit'
mcmc = mc.MCMC([alpha, beta, pi_mean, pi, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)
mcmc.sample(iter*10, burn*10, thin*10)

### @export 'beta-binomial-store'
results['Beta binomial'] = dict(pi=pi_mean.stats(), pred=pred.stats())


### @export 'negative-binomial-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
#delta = mc.Uninformative('delta', value=100.)
mu_log_10_delta = 1.
log_10_delta = mc.Normal('log_10_delta', mu=mu_log_10_delta, tau=.25**-2)
delta = mc.Lambda('delta', lambda x=log_10_delta: .5 + 10**x)

@mc.potential
def obs(pi=pi, delta=delta):
    return mc.negative_binomial_like(r*n, pi*n, delta)

@mc.deterministic
def pred(pi=pi, delta=delta):
    return mc.rnegative_binomial(pi*n_pred, delta) / float(n_pred)

### @export 'negative-binomial-fit-and-store'
mc.MCMC([pi, log_10_delta, delta, obs, pred]).sample(iter, burn, thin)

results['Negative binomial'] = dict(pi=pi.stats(), pred=pred.stats(), delta=delta.stats())


### @export 'log-normal-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)
nonzero_r = pl.maximum(r, .5/n)

@mc.potential
def obs(pi=pi, sigma=sigma, r=nonzero_r):
    return mc.normal_like(pl.log(r), pl.log(pi), 1./((s/r)**2 + sigma**2))

pred_s = pl.sqrt(r * (1-r) / n_pred)
@mc.deterministic
def pred(pi=pi, sigma=sigma):
    s_pred = pl.sqrt(pi*(1-pi)/n_pred)
    return pl.exp(mc.rnormal(pl.log(pi), 1./((s_pred/pi)**2 + sigma**2)))

### @export 'log-normal-fit-and-store'
mc.MCMC([pi, sigma, obs, pred]).sample(iter, burn, thin)

results['Lognormal'] = dict(pi=pi.stats(), pred=pred.stats())


### @export 'offset-log-normal-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
zeta = mc.Uniform('zeta', lower=0, upper=.00005, value=.00001)
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

@mc.potential
def obs(pi=pi, zeta=zeta, sigma=sigma):
    return mc.normal_like(pl.log(r+zeta), pl.log(pi+zeta), 1./((s/(r+zeta))**2 + sigma**2))

@mc.deterministic
def pred(pi=pi, zeta=zeta, sigma=sigma):
    s_pred = pl.sqrt(pi*(1-pi)/n_pred)
    return pl.exp(mc.rnormal(pl.log(pi+zeta),
                    1./((s_pred/(pi+zeta))**2 + sigma**2))) \
                - zeta

### @export 'offset-log-normal-fit-and-store'
mc.MCMC([pi, zeta, sigma, obs, pred]).sample(iter, burn, thin)

results['Offset lognormal'] = dict(pi=pi.stats(), pred=pred.stats())


### @export 'save'
pi_median = []
pi_spread = []
for i, k in enumerate(results):
    pi_median.append(results[k]['pi']['quantiles'][50])
    pi_spread.append(results[k]['pi']['95% HPD interval'][1]-results[k]['pi']['95% HPD interval'][0])
min_est = min(pi_median).round(4)
max_est = max(pi_median).round(4)
min_spread = min(pi_spread).round(4)
max_spread = max(pi_spread).round(4)


book_graphics.save_json('zero_forest.json', vars())

### master graphic of data and models, for zeros subsection of rate models section of stats chapter
book_graphics.forest_plot(r, n,
                          xmax=.0005,
                          model_keys=['Zero-inflated beta binomial', 'Beta binomial', 'Negative binomial', 'Lognormal',  'Offset lognormal'],
                          results=results,
                          pi_true=pi_true,
                          fname='zero_forest.pdf')
