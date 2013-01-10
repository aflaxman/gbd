""" Generate plots to demonstrate properties of negative binomial model"""

### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import rate_model
reload(rate_model)

import book_graphics
reload(book_graphics)


### @export 'negative-binomial-model-funnel'
import simplejson as json
schiz = json.load(open('schiz_forest.json'))

def plot_funnel(pi_true, delta_str):
    n = pl.exp(mc.rnormal(10, 2**-2, size=10000))
    delta = float(delta_str)*pl.ones_like(n)
    p = pi_true*pl.ones_like(n)

    # old way:
    #delta = delta * p * n

    nb = rate_model.neg_binom('funnel', p, delta, p, n)
    r = nb['p_pred'].value

    pl.vlines([pi_true], .1*n.min(), 10*n.max(),
              linewidth=2, linestyle='-', color='w', zorder=9)
    pl.vlines([pi_true], .1*n.min(), 10*n.max(),
              linewidth=1, linestyle='--', color='black', zorder=10)
    pl.plot(r, n, 'ko',
            mew=0, alpha=.25)

    pl.semilogy(schiz['r'], schiz['n'], 'ks', mew=1, mec='white', ms=4,
                label='Observed values')

    pl.xlabel('Rate (per PY)')
    pl.ylabel('Study size (PY)')
    pl.axis([-.0001, .0101, 50., 15000000])
    pl.title(r'$\delta = %s$'%delta_str)

pl.figure(figsize=(11, 8.5), dpi=120)
pl.subplots_adjust(wspace=.4)
pl.subplot(2,2,1)
plot_funnel(.004, '5')
pl.subplot(2,2,2)
plot_funnel(.004, '50')

pl.subplot(2,1,2)
r = pl.array(schiz['r'])
n = pl.array(schiz['n'])

pi = mc.Uniform('pi', 0, 1, value=.001)
delta = mc.Uniform('delta', 0, 100, value=.0001)
nb = rate_model.neg_binom_model('funnel', pi*pl.ones_like(n), delta*pl.ones_like(n), r, n)
# old way:
#nb = rate_model.neg_binom_model('funnel', pi, delta*r*n, r, n)

mcmc = mc.MCMC([pi, delta, nb])
mcmc.sample(20000, 10000, 10, verbose=False, progress_bar=False)

sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(nb['p_pred'].trace()))
for i,s_i in enumerate(sorted_indices):
    pl.plot(s_i+jitter, nb['p_pred'].trace()[:, i], 'ko', mew=0, alpha=.25, zorder=-99)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ws', mew=1, mec='white', ms=5, elinewidth=3, capsize=0)
pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ks', mew=1, mec='white', ms=5)

pl.xticks([])
pl.ylabel('Rate (per PY)')
pl.axis([-.5, 15.5,-.0001,.0121])
pl.savefig('graphics/negative-binomial-funnel.pdf')
pl.savefig('graphics/negative-binomial-funnel.png')


#mc.Matplot.plot(pi)
#mc.Matplot.plot(delta)
print 'pi', pi.stats()['quantiles']
print 'delta', delta.stats()['quantiles']

pl.show()