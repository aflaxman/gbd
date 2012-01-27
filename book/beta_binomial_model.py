""" Generate plots to demonstrate properties of beta-binomial model"""

### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

### @export 'beta-distribution'
pl.figure(**book_graphics.quarter_page_params)

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
import simplejson as json
schiz = json.load(open('schiz_forest.json'))

def plot_beta_binomial_funnel(alpha, beta):
    pi_true = alpha/(alpha+beta)
    pi = mc.rbeta(alpha, beta, size=10000)

    n = pl.exp(mc.rnormal(10, 2**-2, size=10000))
    k = mc.rbinomial(pl.array(n, dtype=int), pi)
    r = k/n
    pl.vlines([pi_true], 0, 10*n.max(),
              linewidth=2, linestyle='--', color='black', zorder=10)
    pl.plot(r, n, 'o',
            mew=0, alpha=.25)

    pl.semilogy(schiz['r'], schiz['n'], 'gs', mew=1, mec='white', ms=4,
                label='Data')

    pl.xlabel('Rate (Per PY)')
    pl.ylabel('Study Size (PY)')
    pl.axis([-.0001, .0101, 50., 1500000])
    pl.title(r'$\alpha=%d$, $\beta=%d$' % (alpha, beta))

pl.figure(**book_graphics.half_page_params)
pl.subplots_adjust(wspace=.4)
pl.subplot(2,2,1)
plot_beta_binomial_funnel(4, 996)

pl.subplot(2,2,2)
plot_beta_binomial_funnel(40., 9960.)

pl.subplot(2,1,2)
r = pl.array(schiz['r'])
n = schiz['n']
k = r*n

alpha = mc.Uninformative('alpha', value=1.)
beta = mc.Uninformative('beta', value=999.)
pi = mc.Beta('pi', alpha, beta, value=.001*pl.ones(16))
pi_mean = mc.Lambda('pi_mean', lambda alpha=alpha, beta=beta: alpha/(alpha+beta))

@mc.potential
def obs(pi=pi):
    return mc.binomial_like(k, n, pi)
@mc.deterministic
def pred(pi=pi, alpha=alpha, beta=beta):
    return mc.rbinomial(n, pi)

mcmc = mc.MCMC([alpha, beta, pi, pi_mean, obs, pred])
mcmc.use_step_method(mc.AdaptiveMetropolis, [alpha, beta])
mcmc.use_step_method(mc.AdaptiveMetropolis, pi)  # TODO: consider making this a Gibbs step
mcmc.sample(200000,100000,100)

sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(pred.trace()))
for i,s_i in enumerate(sorted_indices):
    pl.plot(s_i+jitter, pred.trace()[:, i]/float(n[i]), 'bo', mew=0, alpha=.25, zorder=-99)
    #pl.plot(s_i+jitter, pi_mean.trace(), 'ro', mew=0, alpha=.25, zorder=-100)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='gs', mew=1, mec='white', ms=5)

pl.xticks([])
pl.ylabel('Rate (per PY)')
pl.axis([-.5, 15.5,-.0001,.0121])
pl.savefig('beta-binomial-funnel.png')

mc.Matplot.plot(alpha)
mc.Matplot.plot(beta)
mc.Matplot.plot(pi)



### @export 'beta-binomial-model-fixes-problem'
pop_A_prev = .002
pop_A_N = 50000
pop_B_prev = .006
pop_B_N = 50000

alpha = mc.Uninformative('alpha', value=10*pop_A_prev)
beta = mc.Uninformative('beta', value=10*(1-pop_A_prev))
pi = mc.Beta('pi', alpha, beta, value=[pop_A_prev, pop_B_prev, .02])
@mc.potential
def obs(pi=pi):
    return pop_A_prev*pop_A_N*pl.log(pi[0]) + (1-pop_A_prev)*pop_A_N*pl.log(1-pi[0]) \
        + pop_B_prev*pop_B_N*pl.log(pi[1]) + (1-pop_B_prev)*pop_B_N*pl.log(1-pi[1])
pop_C_N = 50000
pop_C_k = mc.Binomial('pop_C_k', pop_C_N, pi[2])
mc.MCMC([alpha, beta, pi, obs, pop_C_k]).sample(200000,100000,20)

pop_C_prev = pop_C_k.stats()['quantiles'][50] / float(pop_C_N)
pop_C_prev_per_1000 = '%.0f' % (pop_C_prev*1000)
print pop_C_prev_per_1000

pop_C_ui = pop_C_k.stats()['95% HPD interval'] / float(pop_C_N)
pop_C_ui_per_1000 = '[%.0f, %.0f]' % tuple(pop_C_ui*1000)
print pop_C_ui_per_1000


### @export 'save-vars'
book_graphics.save_json('beta_binomial_model.json', vars())




