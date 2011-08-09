""" Generate plots to demonstrate properties of various rate models"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

quarter_page_params = dict(figsize=(8.5,3), dpi=120)
half_page_params = dict(figsize=(8.5, 5.5), dpi=120)
"""

### @export 'binomial-model-funnel'
def decorate_funnel():
    # TODO: overlay real data with same mean in green squares need to
    # select data that has same or similar age intervals, and describe it
    # in the text
    pl.xlabel('Rate (Per PY)')
    pl.ylabel('Study Size (PY)')
    pl.axis([-.0005, .0405, -5000, 105000])

pi_binomial_funnel = .02

n = pl.exp(mc.rnormal(10, 1**-2, size=10000))
n = pl.array(n.round()+1, dtype=int)
k = mc.rbinomial(n, pi_binomial_funnel)
k = pl.array(k, dtype=float)
r = k/n

pl.figure(**half_page_params)
pl.vlines([pi_binomial_funnel], 0, 10*n.max(),
          linewidth=2, linestyle='--', color='black', zorder=10)
pl.plot(r, n, 'o',
        mew=0, alpha=.25)

decorate_funnel()


### @export 'binomial-model-problem'
pop_A_prev = .01
pop_A_N = 1000
pop_B_prev = .03
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
pi_binomial_ppc = .02
ppc_example_dispersion = 2.
log_ppc_n = 10
n = pl.exp(mc.rnormal(log_ppc_n, 1**-2, size=16))
n = pl.array(n.round()+1, dtype=int)
k = pl.minimum(n, mc.rnegative_binomial(pi_binomial_ppc*n, ppc_example_dispersion))
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

pl.figure(**quarter_page_params)
sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(pred.trace()))
for i in sorted_indices:
    pl.plot(i+jitter, pred.trace()[:,i]/float(n[i]), 'bo', mew=0, alpha=.25, zorder=-100)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='gs', mew=0, ms=10)

pl.xticks([])
pl.yticks([0, .02, .04, .06, .08])
pl.ylabel('Rate (per PY)')
pl.axis([-.5, 15.5,-.001,.081])
pl.savefig('binomial-model-ppc.png')


### @export 'beta-distribution'
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


### @export 'beta-binomial-model-funnel'
def plot_beta_binomial_funnel(alpha, beta):
    pi_true = alpha/(alpha+beta)
    pi = mc.rbeta(alpha, beta, size=10000)

    n = pl.exp(mc.rnormal(10, 1**-2, size=10000))
    n = pl.array(n.round()+1, dtype=int)

    k = mc.rbinomial(n, pi)
    k = pl.array(k, dtype=float)
    r = k/n
    pl.vlines([pi_true], 0, 10*n.max(),
              linewidth=2, linestyle='--', color='black', zorder=10)
    pl.plot(r, n, 'o',
            mew=0, alpha=.25)

    decorate_funnel()
    pl.title(r'$\alpha=%d$, $\beta=%d$' % (alpha, beta))

pl.figure(**half_page_params)
pl.subplots_adjust(wspace=.4)
pl.subplot(2,2,1)
plot_beta_binomial_funnel(20., 980.)

pl.subplot(2,2,2)
plot_beta_binomial_funnel(2000., 98000.)

pl.subplot(2,1,2)
n = pl.exp(mc.rnormal(log_ppc_n, 1**-2, size=16))
n = pl.array(n.round()+1, dtype=int)
k = pl.minimum(n, mc.rnegative_binomial(.02*n, ppc_example_dispersion))
k = pl.array(k, dtype=float)
r = k/n

alpha = mc.Uninformative('alpha', value=1.)
beta = mc.Uninformative('beta', value=99.)
pi = mc.Beta('pi', alpha, beta, value=.01*pl.ones(16))
@mc.potential
def obs(pi=pi):
    return mc.binomial_like(k, n, pi)
@mc.deterministic
def pred(pi=pi, alpha=alpha, beta=beta):
    return [mc.rbinomial(n, pi), mc.rbinomial(n, mc.rbeta(alpha, beta))]

mcmc = mc.MCMC([alpha, beta, pi, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)  # TODO: consider making this a Gibbs step
mcmc.sample(200000,100000,100)

sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(pred.trace()))
for i,s_i in enumerate(sorted_indices):
    pl.plot(s_i+jitter, pred.trace()[:, 1, i]/float(n[i]), 'ro', mew=0, alpha=.25, zorder=-100)
    pl.plot(s_i+jitter, pred.trace()[:, 0, i]/float(n[i]), 'bo', mew=0, alpha=.25, zorder=-99)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='gs', mew=0, ms=10)

pl.xticks([])
pl.ylabel('Rate (per PY)')

pl.savefig('beta-binomial-funnel.png')

mc.Matplot.plot(alpha)
mc.Matplot.plot(beta)
mc.Matplot.plot(pi)

"""


### @export 'beta-binomial-model-fixes-problem'
pop_A_prev = .01
pop_A_N = 1000
pop_B_prev = .03
pop_B_N = 1000

alpha = mc.Uninformative('alpha', value=10*pop_A_prev)
beta = mc.Uninformative('beta', value=10*(1-pop_A_prev))
pi = mc.Beta('pi', alpha, beta, value=[pop_A_prev, pop_B_prev, .02])
@mc.potential
def obs(pi=pi):
    return pop_A_prev*pop_A_N*pl.log(pi[0]) + (1-pop_A_prev)*pop_A_N*pl.log(1-pi[0]) \
        + pop_B_prev*pop_B_N*pl.log(pi[1]) + (1-pop_B_prev)*pop_B_N*pl.log(1-pi[1])
pop_C_N = 1000
pop_C_k = mc.Binomial('pop_C_k', pop_C_N, pi[2])
mc.MCMC([alpha, beta, pi, obs, pop_C_k]).sample(20000,10000,20)

pop_C_prev = pop_C_k.stats()['mean'] / float(pop_C_N)
pop_C_prev = pl.round_(pop_C_prev, 2)
print pop_C_prev

pop_C_ui = pop_C_k.stats()['95% HPD interval'] / float(pop_C_N)
pop_C_ui = '[%.0f, %.0f]' % tuple(pop_C_ui*100)
print pop_C_ui


### @export 'save-vars'
pl.show()
import simplejson as json
json.dump(dict(vars()), open('rate_model.json', 'w'),
          indent=2, skipkeys=True, default=lambda x: None)





