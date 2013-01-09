""" Generate plots to demonstrate properties of various rate models"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

n_small = 500
pi_true = .025
delta_true = 5.

n = pl.array(pl.exp(mc.rnormal(10, 1**-2, size=16)), dtype=int)
k = pl.array(mc.rnegative_binomial(n*pi_true, delta_true), dtype=float)
r = k/n

iter = 20000
burn = 10000
thin = 10
results = {}
xmax = .07

### @export 'distribution-comparison'
pl.figure(**book_graphics.quarter_page_params)

ax = pl.axes([.1, .3, .85, .65])
x = pl.arange(0, n_small*pi_true*4, .1)

# plot binomial distribution
y1 = [pl.exp(mc.binomial_like(x_i, n_small, pi_true)) for x_i in x]
pl.step(x, y1, 'k',
        linewidth=1, linestyle='step:', alpha=.8,
        label='Binomial')

# plot poisson distribution
y2 = [pl.exp(mc.poisson_like(x_i, n_small*pi_true)) for x_i in x]
pl.plot(x, y2, 'k',
        linewidth=1, linestyle='steps--', alpha=.8,
        label='Poisson')

pl.legend(loc='upper right', fancybox=True, shadow=True)
pl.yticks([0, .05])
pl.xticks([25, 50, 75], ['','',''])
pl.axis([-.1, n_small*pi_true*4, -.02, 1.1*max(y1)])
pl.xlabel('Count')
#pl.figtext(.11, .94, 'Distribution', ha='left', va='top')

ax = pl.axes([.1, .15, .85, .20])
pl.step(x, pl.array(y1)-y2, color='k')
pl.xticks([0, 25, 50, 75])
pl.yticks([-.001, .001])
pl.axis([-.1, n_small*pi_true*4, -.0015, .0015])
pl.xlabel('Count')
pl.figtext(.11, .34, 'Binomial $-$ Poisson', ha='left', va='top')

pl.savefig('poisson_approx_to_binom.pdf')

### @export 'poisson-model'
n_pred = 1.e9

pi = mc.Uniform('pi', lower=0, upper=1, value=.5)

@mc.potential
def obs(pi=pi):
    return mc.poisson_like(r*n, pi*n)

@mc.deterministic
def pred(pi=pi):
    return mc.rpoisson(pi*n_pred) / float(n_pred)

### @export 'poisson-fit-and-store'
mc.MCMC([pi, obs, pred]).sample(iter, burn, thin, verbose=False, progress_bar=False)

results['Poisson'] = dict(pred=pred.stats(), pi=pi.stats())



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
mc.MCMC([pi, delta, obs, pred]).sample(iter, burn, thin, verbose=False, progress_bar=False)

key = 'Negative Binomial'
results[key] = dict(pred=pred.stats(), pi=pi.stats())

#mc.Matplot.plot(pi)
#mc.Matplot.plot(delta)


model_keys = ['Poisson', 'Negative Binomial']
### @export 'negative-binomial_dispersion-prior-exploration'
for mu_log_10_delta in [1,2,3]:
    ### @export 'negative-binomial_dispersion-alt_prior'
    pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
    log_10_delta = mc.Normal('log_10_delta', mu=mu_log_10_delta, tau=.25**-2)
    delta = mc.Lambda('delta', lambda x=log_10_delta: 5 + 10**x)

    @mc.potential
    def obs(pi=pi, delta=delta):
        return mc.negative_binomial_like(r*n, pi*n, delta)

    @mc.deterministic
    def pred(pi=pi, delta=delta):
        return mc.rnegative_binomial(pi*n_pred, delta) / float(n_pred)

    ### @export 'negative-binomial_dispersion-fit_alt_prior'
    mc.MCMC([pi, log_10_delta, delta, obs, pred]).sample(iter, burn, thin, verbose=False, progress_bar=False)

    key = 'Neg. Binom., $\mu_{\log\delta}=%d$'%mu_log_10_delta
    results[key] = dict(pred=pred.stats(), pi=pi.stats())
    model_keys.append(key)

    #mc.Matplot.plot(pi)
    #mc.Matplot.plot(delta)


book_graphics.save_json('poisson_model.json', vars())
book_graphics.forest_plot(fname='neg_binom_priors.pdf', **vars())

pl.show()