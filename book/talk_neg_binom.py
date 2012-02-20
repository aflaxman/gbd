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


colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f0', '#ffff33']

### @export 'negative-binomial-model-funnel'
import simplejson as json
schiz = json.load(open('schiz_forest.json'))

def plot_joint_density(X, Y, bounds=None):
    if bounds:
        X_min, X_max, Y_min, Y_max = bounds
    else:
        X_min = X.min()
        X_max = X.max()
        Y_min = Y.min()
        Y_max = Y.max()

    pl.plot(X, Y,
         linestyle='none', marker='o', color='green', mec='green',
         alpha=.75, zorder=-99)

    import scipy.stats
    gkde = scipy.stats.gaussian_kde([X, Y])
    x,y = pl.mgrid[X_min:X_max:(X_max-X_min)/25.,
                   Y_min:Y_max:(Y_max-Y_min)/25.]
    z = pl.array(gkde.evaluate([x.flatten(), y.flatten()])).reshape(x.shape)
    pl.contour(x, y, z, linewidths=2)

    pl.axis([X_min, X_max, Y_min, Y_max])


def plot_funnel(pi_true, delta_str):
    delta = float(delta_str)
    n = pl.exp(mc.rnormal(10, 2**-2, size=10000))
    p = pi_true*pl.ones_like(n)

    # old way:
    #delta = delta * p * n

    nb = rate_model.neg_binom_model('funnel', pi_true, delta, p, n)
    r = nb['p_pred'].value

    pl.vlines([pi_true], .1*n.min(), 10*n.max(),
              linewidth=5, linestyle='--', color='black', zorder=10)
    pl.plot(r, n, 'o', color=colors[0], ms=10,
            mew=0, alpha=.25)

    pl.semilogy(schiz['r'], schiz['n'], 's', mew=1, mec='white', ms=15,
                color=colors[1],
                label='Observed Values')

    pl.xlabel('Rate (Per 1000 PY)', size=32)
    pl.ylabel('Study Size (PY)', size=32)
    pl.axis([-.0001, .0101, 50., 15000000])
    pl.title(r'$\delta = %s$'%delta_str, size=48)
    pl.xticks([0, .005, .01], [0, 5, 10], size=30)
    pl.yticks(size=30)


pl.figure(figsize=(17., 11), dpi=72)
pl.subplots_adjust(wspace=.41, hspace=.45)
pl.subplot(2,2,1)
plot_funnel(.004, '5')
pl.subplot(2,2,2)
plot_funnel(.004, '50')

pl.subplot(2,1,2)
r = pl.array(schiz['r'])
n = pl.array(schiz['n'])

pi = mc.Uniform('pi', 0, 1, value=.001)
delta = mc.Uniform('delta', 0, 100, value=.0001)
nb = rate_model.neg_binom_model('funnel', pi, delta, r, n)
# old way:
#nb = rate_model.neg_binom_model('funnel', pi, delta*r*n, r, n)

mcmc = mc.MCMC([pi, delta, nb])
mcmc.sample(20000, 10000, 10)

sorted_indices = r.argsort().argsort()
jitter = mc.rnormal(0, .1**-2, len(nb['p_pred'].trace()))
for i,s_i in enumerate(sorted_indices):
    pl.plot(s_i+jitter, nb['p_pred'].trace()[:, i], 'o', mew=0, alpha=.25, zorder=-99, color=colors[2], ms=10)

#pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='ws', mew=1, mec='white', ms=15, elinewidth=10, capsize=0)
pl.errorbar(sorted_indices, r, yerr=1.96*pl.sqrt(r*(1-r)/n), fmt='s', mew=1, mec='white', ms=15, elinewidth=10, capsize=0,
            color=colors[3])

pl.xticks([])
pl.yticks([0, .005, .01], [0, 5, 10], size=30)
pl.title('Posterior Predictions', size=48)
pl.ylabel('Rate (per 1000 PY)', size=32)
pl.axis([-.5, 15.5,-.0001,.0121])
pl.savefig('/media/windows/t/negative-binomial-funnel.png')
pl.show()

#mc.Matplot.plot(pi)
#mc.Matplot.plot(delta)
print 'pi', pi.stats()['quantiles']
print 'delta', delta.stats()['quantiles']
