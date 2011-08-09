""" Fit data with several rate models and generate forest plot"""


import sys
sys.path += ['..']


import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

results = {}
n_pred = 10000
iter = 20000
burn = 10000
thin = 10

### @export 'data'
dm = dismod3.load_disease_model(15630)
dm.calc_effective_sample_size(dm.data)
some_data = ([d for d in dm.data
              if d['data_type'] == 'prevalence data'
              and d['sex'] == 'male'
              and 15 <= d['age_start'] < 20
              and d['age_end'] == 99
              and d['effective_sample_size'] > 1])
countries = pl.unique([s['region'] for s in some_data])
min_year = min([s['year_start'] for s in some_data])
max_year = max([s['year_end'] for s in some_data])
cy = ['%s-%d'%(s['region'], s['year_start']) for s in some_data]

n = pl.array([s['effective_sample_size'] for s in some_data])
r = pl.array([dm.value_per_1(s) for s in some_data])

### @export 'binomial-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(r*n, n, pi)

@mc.deterministic
def pred(pi=pi):
    return mc.rbinomial(n_pred, pi) / float(n_pred)

### @export 'binomial-fit'
mc.MCMC([pi, obs, pred]).sample(iter, burn, thin)

### @export 'binomial-store'
mc.Matplot.plot(pi)
pl.savefig('ci-prev_meta_analysis-binomial_diagnostic.png')
results['binomial'] = pred.stats()
results['binomial']['trace'] = list(pred.trace())


### @export 'beta-binomial-model'
alpha = mc.Uninformative('alpha', value=4.)
beta = mc.Uninformative('beta', value=1000.)
pi = mc.Beta('pi', alpha, beta, value=r)

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(r*n, n, pi)

@mc.deterministic
def pred(alpha=alpha, beta=beta):
    return mc.rbinomial(n_pred, mc.rbeta(alpha, beta)) / float(n_pred)

### @export 'beta-binomial-fit'
mcmc = mc.MCMC([alpha, beta, pi, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)
mcmc.sample(iter, burn, thin)

### @export 'beta-binomial-store'
#mc.Matplot.plot(alpha)
#mc.Matplot.plot(beta)
mc.Matplot.plot(pi)
pl.savefig('ci-prev_meta_analysis-beta_binomial_diagnostic.png')
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

### @export 'poisson-fit-and-store'
mc.MCMC([pi, obs, pred]).sample(iter, burn, thin)

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

### @export 'negative-binomial-fit-and-store'
mc.MCMC([pi, delta, obs, pred]).sample(iter, burn, thin)

results['negative-binomial'] = pred.stats()
results['negative-binomial']['trace'] = list(pred.trace())


### @export 'negative-binomial_dispersion-prior-exploration'
results['negative-binomial_w-prior'] = {}
for mu_log_10_delta in [1,2,3]:
    pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
    ### @export 'negative-binomial_alt-prior'
    log_10_delta = mc.Normal('log_10_delta', mu=mu_log_10_delta, tau=.25**-2)

    @mc.potential
    def obs(pi=pi, log_10_delta=log_10_delta):
        return mc.negative_binomial_like(r*n, pi*n, 10**log_10_delta)

    ### @export 'negative-binomial_exploration-continues'
    @mc.deterministic
    def pred(pi=pi, log_10_delta=log_10_delta):
        return mc.rnegative_binomial(pi*n_pred, 10**log_10_delta) / float(n_pred)

    mc.MCMC([pi, log_10_delta, obs, pred]).sample(iter, burn, thin)

    results['negative-binomial_w-prior'][log_10_delta] = pred.stats()
    results['negative-binomial_w-prior'][log_10_delta]['trace'] = list(pred.trace())

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

### @export 'normal-fit-and-store'
mc.MCMC([pi, sigma_squared, obs, pred]).sample(iter, burn, thin)

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

### @export 'log-normal-fit-and-store'
mc.MCMC([pi, sigma, obs, pred]).sample(iter, burn, thin)

results['log-normal'] = pred.stats()
results['log-normal']['trace'] = list(pred.trace())


### @export 'offset-log-normal-model'
sampling_variance = 0 # FIXME

pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
zeta = mc.Uniform('zeta', lower=0, upper=.005, value=.001)
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

### @export 'offset-log-normal-fit-and-store'
mc.MCMC([pi, zeta, sigma, obs, pred]).sample(iter, burn, thin)

results['offset-log-normal'] = pred.stats()
results['offset-log-normal']['trace'] = list(pred.trace())


### @export 'save'
pi_median = []
pi_spread = []
for i, k in enumerate('binomial beta-binomial poisson negative-binomial normal log-normal offset-log-normal'.split()):
    pi_median.append(results[k]['quantiles'][50])
    pi_spread.append(results[k]['95% HPD interval'][1]-results[k]['95% HPD interval'][0])
min_est = min(pi_median).round(4)
max_est = max(pi_median).round(4)
min_spread = min(pi_spread).round(4)
max_spread = max(pi_spread).round(4)


def array_to_list(x):
    try:
        return list(x)
    except:
        return None

pl.show()
import simplejson as json
json.dump(dict(vars()), open('schiz_forest.json', 'w'),
          indent=2, skipkeys=True, default=array_to_list)
