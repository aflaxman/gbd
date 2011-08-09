""" Fit data with several rate models and generate forest plot"""

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
n_pred = 10000

### @export 'data'
# TODO: replace with real data
pi_true = .02
delta_true = 2.
log_n = 10

n = pl.array(pl.exp(mc.rnormal(log_n, 1**-2, size=16)), dtype=int)
r = pl.array(pl.minimum(n, mc.rnegative_binomial(pi_true*n, delta_true)), dtype=float)/n


### @export 'binomial-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(r*n, n, pi)

@mc.deterministic
def pred(pi=pi):
    return mc.rbinomial(n_pred, pi) / float(n_pred)

mc.MCMC([pi, obs, pred]).sample(20000,10000,10)

results['binomial'] = pred.stats()
results['binomial']['trace'] = list(pred.trace())


### @export 'beta-binomial-model'
alpha = mc.Uninformative('alpha', value=1.)
beta = mc.Uninformative('beta', value=99.)
pi = mc.Beta('pi', alpha, beta, value=.01*pl.ones(16))

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(r*n, n, pi)

@mc.deterministic
def pred(alpha=alpha, beta=beta):
    return mc.rbinomial(n_pred, mc.rbeta(alpha, beta)) / float(n_pred)

mcmc = mc.MCMC([alpha, beta, pi, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)
mcmc.sample(200000,100000,100)

results['beta-binomial'] = pred.stats()
results['beta-binomial']['trace'] = list(pred.trace())


### @export 'poisson-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)

@mc.potential
def obs(pi=pi):
    return mc.poisson_like(r*n, pi*n)

@mc.deterministic
def pred(pi=pi):
    return mc.rpoisson(pi*n_pred) / float(n_pred)

mc.MCMC([pi, obs, pred]).sample(20000,10000,10)

results['poisson'] = pred.stats()
results['poisson']['trace'] = list(pred.trace())


### @export 'negative-binomial-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
delta = mc.Uninformative('delta', value=100.)

@mc.potential
def obs(pi=pi, delta=delta):
    return mc.negative_binomial_like(r*n, pi*n, delta)

@mc.deterministic
def pred(pi=pi, delta=delta):
    return mc.rnegative_binomial(pi*n_pred, delta) / float(n_pred)

mc.MCMC([pi, delta, obs, pred]).sample(20000,10000,10)

results['negative-binomial'] = pred.stats()
results['negative-binomial']['trace'] = list(pred.trace())


### @export 'normal-model'
sampling_variance = r * (1-r) / n

pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
sigma_squared = mc.Uninformative('sigma_squared', value=.01)

@mc.potential
def obs(pi=pi, sigma_squared=sigma_squared):
    return mc.normal_like(r, pi, 1./(sampling_variance + sigma_squared))

@mc.deterministic
def pred(pi=pi, sigma_squared=sigma_squared):
    return mc.rnormal(pi, 1./(pi*(1-pi)/n_pred + sigma_squared))

mc.MCMC([pi, sigma_squared, obs, pred]).sample(20000,10000,10)

results['normal'] = pred.stats()
results['normal']['trace'] = list(pred.trace())


### @export 'log-normal-model'
sampling_variance = 0 # TODO: approximate this, something like 1 / (r * (1-r) / n)

pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

@mc.potential
def obs(pi=pi, sigma=sigma):
    return mc.normal_like(pl.log(r), pl.log(pi), 1./(sampling_variance + sigma**2))

@mc.deterministic
def pred(pi=pi, sigma=sigma):
    pred_sampling_variance = 0 # FIXME
    return pl.exp(mc.rnormal(pl.log(pi), 1./(pred_sampling_variance + sigma**2)))

mc.MCMC([pi, sigma, obs, pred]).sample(20000,10000,10)

results['log-normal'] = pred.stats()
results['log-normal']['trace'] = list(pred.trace())


### @export 'offset-log-normal-model'
sampling_variance = 0 # FIXME

pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
zeta = mc.Uniform('zeta', lower=0, upper=.005, value=.001)
#zeta = .001
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

@mc.potential
def obs(pi=pi, zeta=zeta, sigma=sigma):
    return mc.normal_like(pl.log(r+zeta), pl.log(pi+zeta), 1./(sampling_variance + sigma**2))

@mc.deterministic
def pred(pi=pi, zeta=zeta, sigma=sigma):
    pred_sampling_variance = 0 # FIXME
    return pl.exp(mc.rnormal(pl.log(pi+zeta),
                    1./(pred_sampling_variance + sigma**2))) \
                - zeta

mc.MCMC([pi, zeta, sigma, obs, pred]).sample(20000,10000,10)

results['offset-log-normal'] = pred.stats()
results['offset-log-normal']['trace'] = list(pred.trace())


### @export 'save'

def array_to_list(x):
    try:
        return list(x)
    except:
        return None

pl.show()
import simplejson as json
json.dump(dict(vars()), open('schiz_forest.json', 'w'),
          indent=2, skipkeys=True, default=array_to_list)
