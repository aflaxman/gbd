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

s = pl.sqrt(r * (1-r) / n)


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
results['Binomial'] = dict(pi=pi.stats(), pred=pred.stats())


### @export 'beta-binomial-model'
alpha = mc.Uninformative('alpha', value=4.)
beta = mc.Uninformative('beta', value=1000.)
pi_mean = mc.Lambda('pi_mean', lambda alpha=alpha, beta=beta: alpha/(alpha+beta))
pi = mc.Beta('pi', alpha, beta, value=r)

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
#mc.Matplot.plot(alpha)
#mc.Matplot.plot(beta)
mc.Matplot.plot(pi)
pl.savefig('ci-prev_meta_analysis-beta_binomial_diagnostic.png')
results['Beta binomial'] = dict(pi=pi_mean.stats(), pred=pred.stats())


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

results['Poisson'] = dict(pi=pi.stats(), pred=pred.stats())


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

results['Negative binomial'] = dict(pi=pi.stats(), pred=pred.stats())


### @export 'normal-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

@mc.potential
def obs(pi=pi, sigma=sigma):
    return mc.normal_like(r, pi, 1./(s**2 + sigma**2))

@mc.deterministic
def pred(pi=pi, sigma=sigma):
    s_pred = pl.sqrt(pi*(1-pi)/n_pred)
    return mc.rnormal(pi, 1./(s_pred + sigma))

### @export 'normal-fit-and-store'
mc.MCMC([pi, sigma, obs, pred]).sample(iter, burn, thin)

results['Normal'] = dict(pi=pi.stats(), pred=pred.stats())


### @export 'log-normal-model'
pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

@mc.potential
def obs(pi=pi, sigma=sigma):
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
zeta = mc.Uniform('zeta', lower=0, upper=.005, value=.001)
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


book_graphics.save_json('schiz_forest.json', vars())

### data only plot, for computational infrastructure appendix
book_graphics.forest_plot(r, n, data_labels=cy,
                          xmax=.0115,
                          subplot_params=dict(bottom=.1, right=.99, top=.95, left=.15),
                          figparams=book_graphics.quarter_page_params,
                          fname='ci-prev_meta_analysis-schiz_data.png')


### master graphic of data and models, for rate model section of stats chapter
book_graphics.forest_plot(r, n, data_labels=cy,
                          xmax=.0115,
                          model_keys=['Binomial', 'Poisson', 'Beta binomial', 'Negative binomial', 'Normal', 'Lognormal',  'Offset lognormal'],
                          results=results,
                          #subplot_params=dict(bottom=.1, right=.99, top=.95, left=.15),
                          #figparams=book_graphics.quarter_page_params,
                          fname='theory-rate_model-schiz_forest_plot.png')
