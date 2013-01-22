""" Generate plots to demonstrate properties of offset log-normal model"""

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


### @export 'offset-log-normal-funnel'
import simplejson as json
schiz = json.load(open('schiz_forest.json'))

def plot_funnel(pi_true, sigma_str):
    n = pl.exp(mc.rnormal(10, 2**-2, size=10000))
    sigma = float(sigma_str)*pl.ones_like(n)
    p = pi_true*pl.ones_like(n)

    oln = rate_model.offset_log_normal('funnel', p, sigma, p, pl.sqrt(p*(1-p)/n))
    r = oln['p_pred'].value

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
    pl.xticks([0, .005, .01])
    pl.axis([-.0001, .0101, 50., 15000000])
    pl.title(r'$\sigma = %s$'%sigma_str)

pl.figure(figsize=(11, 8.5), dpi=120)
pl.subplots_adjust(wspace=.4)
pl.subplot(2,2,1)
plot_funnel(.004, '0.5')
pl.subplot(2,2,2)
plot_funnel(.004, '0.1')

pl.subplot(2,1,2)
r = pl.array(schiz['r'])
n = pl.array(schiz['n'])

pi = mc.Uniform('pi', 0, 1, value=.001)
sigma = mc.Uniform('sigma', 0, 100, value=.0001)
oln = rate_model.offset_log_normal('funnel', pi*pl.ones_like(n), sigma*pl.ones_like(n), r, pl.sqrt(r*(1-r)/n))

mcmc = mc.MCMC([pi, sigma, oln])
mcmc.sample(20000, 10000, 10, verbose=False, progress_bar=False)

sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(oln['p_pred'].trace()))
for i,s_i in enumerate(sorted_indices):
    pl.plot(s_i+jitter, oln['p_pred'].trace()[:, i], 'ko', mew=0, alpha=.25, zorder=-99)

pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ws', mew=1, mec='white', ms=5, elinewidth=3, capsize=0)
pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ks', mew=1, mec='white', ms=5)

pl.xticks([])
pl.ylabel('Rate (per PY)')
pl.axis([-.5, 15.5,-.0001,.0121])
pl.subplots_adjust(hspace=.35)
pl.savefig('book/graphics/offset-log-normal-funnel.pdf')
pl.savefig('book/graphics/offset-log-normal-funnel.png')


#mc.Matplot.plot(pi)
#mc.Matplot.plot(delta)
print 'pi', pi.stats()['quantiles']
print 'sigma', sigma.stats()['quantiles']
