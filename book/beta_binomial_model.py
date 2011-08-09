""" Generate plots to demonstrate properties of various rate models"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

"""
### @export 'binomial-model-funnel'
pi_binomial_funnel = .1

n = pl.exp(mc.rnormal(10, 1**-2, size=10000))
n = pl.array(n.round()+1, dtype=int)
k = mc.rbinomial(n, pi_binomial_funnel)
k = pl.array(k, dtype=float)
r = k/n

pl.figure()
pl.vlines([pi_binomial_funnel], 0, 10*n.max(),
          linewidth=2, linestyle='--', color='black', zorder=10)
pl.plot(r, n, 'o',
        mew=0, alpha=.25)

# TODO: overlay real data with same mean in green squares need to
# select data that has same or similar age intervals, and describe it
# in the text
pl.xlabel('Rate (Per PY)')
pl.ylabel('Study Size (PY)')
pl.axis([r.min()-.01, r.max()+.01, -1000, 101000])
pl.savefig('binomial-model-funnel.png')
pl.show()

### @export 'binomial-model-problem'
pop_A_prev = .1
pop_A_N = 1000
pop_B_prev = .3
pop_B_N = 1000

pi = mc.Uninformative('pi', value=pop_A_prev)
@mc.potential
def obs(pi=pi):
    return pop_A_prev*pop_A_N*pl.log(pi) + (1-pop_A_prev)*pop_A_N*pl.log(1-pi) \
        + pop_B_prev*pop_B_N*pl.log(pi) + (1-pop_B_prev)*pop_B_N*pl.log(1-pi)
pop_C_N = 1000
pop_C_k = mc.Binomial('pop_C_k', pop_C_N, pi)
mc.MCMC([pi, obs, pop_C_k]).sample(20000,10000,2)

pop_C_prev = pop_C_k.stats()['mean'] / float(pop_C_N)
pop_C_prev = pl.round_(pop_C_prev, 2)

pop_C_ui = pop_C_k.stats()['95% HPD interval'] / float(pop_C_N)
pop_C_ui = '[%.0f, %.0f]' % tuple(pop_C_ui*100)


### @export 'binomial-model-ppc'
pi_binomial_ppc = .1
n = pl.exp(mc.rnormal(6, 1**-2, size=16))
n = pl.array(n.round()+1, dtype=int)
k = pl.minimum(n, mc.rnegative_binomial(pi_binomial_ppc*n, 2.))
k = pl.array(k, dtype=float)
r = k/n

pi = mc.Uninformative('pi', value=.5)
mc.binomial_like(k, n, pi)
@mc.potential
def obs(pi=pi):
    return mc.binomial_like(k, n, pi)
@mc.deterministic
def pred(pi=pi):
    return mc.rbinomial(n, pi)

mc.MCMC([pi, obs, pred]).sample(20000,10000,10)

pl.figure(figsize=(8.5,3), dpi=120)
sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(pred.trace()))
for i in sorted_indices:
    pl.plot(i+jitter, pred.trace()[:,i]/float(n[i]), 'bo', mew=0, alpha=.25)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='gs', mew=0, ms=10)

pl.xticks([])
pl.ylabel('Rate (per PY)')
pl.savefig('binomial-model-ppc.png')
"""


### @export 'beta-distribution'
quarter_page_params = dict(figsize=(8.5,3), dpi=120)
pl.figure(**quarter_page_params)

d_pi = .01
pi = pl.arange(d_pi, 1, d_pi)
def plot_beta(alpha, beta, **params):
    pl.plot(pi, [pl.exp(mc.beta_like(pi_i, alpha, beta)) for pi_i in pi],
            label=r'$\alpha = %.1f$, $\beta = %.1f$'%(alpha, beta),
            linewidth=2,
            alpha=.9,
            **params)

def decorate(mean):
    pl.legend(loc='upper center', bbox_to_anchor=(.5,-.1))
    xmin, xmax, ymin, ymax = pl.axis()
    pl.vlines([mean], -ymax, ymax*10, linestyle='dashed', zorder=20)
    pl.xticks([0, .5, 1])
    pl.yticks([])
    pl.axis([-.01, 1.01, -ymax*.05, ymax*1.01])

pl.subplot(1,2,1)
plot_beta(.5,.5)
plot_beta(1,1)
plot_beta(2,2)
plot_beta(10,10)
decorate(mean=.5)

pl.subplot(1,2,2)
plot_beta(.5,1.5)
plot_beta(1,3)
plot_beta(2,6)
plot_beta(10,30)
decorate(mean=.25)

pl.subplots_adjust(top=.95, bottom=.6)
pl.savefig('beta-distribution.png')

### @export 'save-vars'
pl.show()
import simplejson as json
json.dump(dict(vars()), open('rate_model.json', 'w'),
          indent=2, skipkeys=True, default=lambda x: None)





