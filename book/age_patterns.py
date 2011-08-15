""" Explore features of age pattern priors"""


### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc

import dismod3
import book_graphics
reload(book_graphics)


### @export 'models-of-varying-smoothness'
data = pl.array([[10., .1, .025],
                 [50., .25, .025]])
smoothness=200 # very smooth
smoothness=100 # moderately smooth
smoothness=50 # slightly smooth

## setup lognormal prior on intercept
gamma = mc.Normal('gamma', 0.,
                  2.**-2,
                  value=pl.log(data[:,1].mean()))
mu = mc.Lambda('mu', lambda gamma=gamma: pl.exp(gamma))

## setup Gaussian Process prior on age pattern
mesh = pl.arange(0,101,1)
@mc.deterministic
def M(mu=mu):
    return mc.gp.Mean(lambda x: mu*pl.ones(len(x)))
C = mc.gp.FullRankCovariance(mc.gp.matern.euclidean,
                             amp=data[:,1].max(),
                             scale=smoothness,
                             diff_degree=2)
sm = mc.gp.GPSubmodel('sm', M, C, mesh,
                      init_vals=mu.value*pl.ones_like(mesh))

## condition on rate being positive
@mc.potential
def positive(f=sm.f_eval):
    if pl.any(f < 0.):
        return -pl.inf
    else:
        return 0.

## likelihood of observed data, using normal model for simplicity
@mc.deterministic
def data_expected(mu=mu, sm=sm):
    return sm.f(data[:,0])
@mc.observed
def data_obs(data_expected=data_expected, value=data):
    return mc.normal_like(value[:,1],
                          data_expected,
                          value[:,2]**-2)

## generate good initial values with MAP fit
mc.MAP([gamma, data_obs]).fit(method='fmin_powell', verbose=1)
mc.MAP([sm.f_eval, data_obs]).fit(method='fmin_powell', verbose=1)
mc.MAP([gamma, data_obs]).fit(method='fmin_powell', verbose=1)

## combine all PyMC nodes in to an MCMC object
mcmc = mc.MCMC([gamma, mu, M, sm, positive, data_expected, data_obs])
mcmc.use_step_method(mc.gp.GPParentAdaptiveMetropolis, gamma)
mcmc.sample(10000, 5000, 10)

### @export 'plot-varying-smoothing'
mc.Matplot.plot(data_expected)
mc.Matplot.plot(mu)

pl.figure()
ages = pl.arange(0,101,1)

for i, f_n in enumerate(sm.f_eval.trace()):
    pl.plot(ages, f_n,
         '-', color='grey',
            zorder=-1)
    if i % (len(sm.f_eval.trace())/4) == 0:
        pl.plot(ages, f_n,
                '-', linewidth=2,
                color=pl.cm.spectral(.05+float(i)/len(sm.f_eval.trace())),
                zorder=0)
        
pl.plot(ages, sm.f_eval.stats()['quantiles'][50],
        'k-', linewidth=3, label='Median')
pl.plot(ages, sm.f_eval.stats()['mean'],
        'k--', linewidth=3, label='Mean')
pl.errorbar(data[:,0], data[:,1], yerr=data[:,2]*1.96,
            fmt='gs',
            mec='white', mew=1, ms=5)
pl.legend()

pl.savefig('smoothness_priors.pdf')

pl.show()

### @export 'store-results'

book_graphics.save_json('age_patterns.json', vars())



"""

#eps_p_f = mc.Lambda('eps_p_f', lambda mu=mu: data[:,1]-mu)
#sm.f_eval.children.add(eps_p_f) # WORKAROUND: this shouldn't be necess.

mcmc = mc.MCMC([sm, gamma, mu, f, positive, obs])
mcmc.use_step_method(mc.gp.GPParentAdaptiveMetropolis, gamma)
import step_methods
reload(step_methods)
#mcmc.use_step_method(step_methods.GPEvaluationGibbs, sm,
#                     eps_p_f=eps_p_f, ti=data[:,0], V=data[:,2]**2)
"""
