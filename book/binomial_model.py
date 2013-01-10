""" Generate plots to demonstrate properties of various rate models"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)

### @export 'binomial-model-funnel'
pi_binomial_funnel = .004

n = pl.exp(mc.rnormal(10, 2**-2, size=10000))
k = mc.rbinomial(pl.array(n.round(), dtype=int), pi_binomial_funnel)
r = k/n

pl.figure(**book_graphics.half_page_params)
pl.vlines([pi_binomial_funnel], .1*n.min(), 10*n.max(),
          linewidth=2, linestyle='-', color='w', zorder=9, label='_nolegend_')
pl.vlines([pi_binomial_funnel], .1*n.min(), 10*n.max(),
          linewidth=1, linestyle='--', color='k', zorder=10, label='$\pi$')
pl.plot(r, n, 'ko',
        mew=0, alpha=.25,
        label='Predicted distribution')


import simplejson as json
schiz = json.load(open('schiz_forest.json'))
pl.semilogy(schiz['r'], schiz['n'], 'ks', mew=1, mec='white', ms=8,
            label='Observed value')


pl.xlabel('Rate ($r$)', fontsize='x-large')
pl.ylabel('Study Size ($n$)', fontsize='x-large')
pl.axis([-.0001, .0101, 50., 1500000])
pl.legend(numpoints=1, fancybox=True, shadow=True, prop={'size':'x-large'})
pl.subplots_adjust(bottom=.13, top=.93)
pl.savefig('graphics/binomial-model-funnel.pdf')
pl.savefig('graphics/binomial-model-funnel.png')


### @export 'binomial-model-problem'
n = 50000
pop_A_prev = .002
pop_A_N = n
pop_B_prev = .006
pop_B_N = n

pi = mc.Uninformative('pi', value=pop_A_prev)
@mc.potential
def obs(pi=pi):
    return pop_A_prev*pop_A_N*pl.log(pi) + (1-pop_A_prev)*pop_A_N*pl.log(1-pi) \
        + pop_B_prev*pop_B_N*pl.log(pi) + (1-pop_B_prev)*pop_B_N*pl.log(1-pi)
pop_C_N = n
pop_C_k = mc.Binomial('pop_C_k', pop_C_N, pi)
mc.MCMC([pi, obs, pop_C_k]).sample(20000,10000,2, verbose=False, progress_bar=False)

pop_C_prev = pop_C_k.stats()['quantiles'][50] / float(pop_C_N)
pop_C_prev_per_1000 = '%.0f' % (pop_C_prev*1000)
print pop_C_prev_per_1000

pop_C_ui = pop_C_k.stats()['95% HPD interval'] / float(pop_C_N)
pop_C_ui_per_1000 = '[%.0f, %.0f]' % tuple(pop_C_ui*1000)
print pop_C_ui_per_1000




### @export 'binomial-model-ppc'
r = pl.array(schiz['r'])
n = pl.array(schiz['n'], dtype=int)
k = r*n

pi = mc.Uninformative('pi', value=.5)
mc.binomial_like(k, n, pi)
@mc.potential
def obs(pi=pi):
    return mc.binomial_like(k, n, pi)
@mc.deterministic
def pred(pi=pi):
    return mc.rbinomial(n, pi)

mc.MCMC([pi, obs, pred]).sample(20000,10000,10, verbose=False, progress_bar=False)

pl.figure(**book_graphics.quarter_page_params)
sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(pred.trace()))
for i, s_i in enumerate(sorted_indices):
    if i == 0:
        label = 'Predicted distribution'
    else:
        label = '_nolabel_'
    pl.plot(s_i+jitter, pred.trace()[:,i]/float(n[i]), 'ko', mew=0, alpha=.25, zorder=-100, label=label)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ks', mew=1, ms=5, mec='white', label='Observed value')

pl.xticks([])
pl.yticks([0, .002, .004, .006, .008, .01], fontsize='large')
pl.ylabel('Rate ($r$)', fontsize='x-large')
pl.axis([-.5, 15.5,-.0001,.0121])
pl.legend(loc='upper left', numpoints=1, fancybox=True, shadow=True)
pl.savefig('graphics/binomial-model-ppc.pdf')
pl.savefig('graphics/binomial-model-ppc.png')


### @export 'save-vars'
book_graphics.save_json('binomial_model.json', vars())

pl.show()


