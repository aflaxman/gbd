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
## setup lognormal prior on intercept
gamma = mc.Normal('gamma', 0., 2.**-2, value=0.)
mu = mc.Lambda('mu', lambda gamma=gamma: pl.exp(gamma))

## setup Gaussian Process prior on age pattern
mesh = pl.arange(0,101,1) #dismod3.settings.gbd_ages
M = mc.gp.Mean(lambda x: pl.zeros(len(x)))
C = mc.gp.FullRankCovariance(mc.gp.matern.euclidean,
                             amp=.1,
                             scale=20., 
                             diff_degree=2)
sm = mc.gp.GPSubmodel('sm', M, C, mesh,
                      init_vals=pl.zeros_like(mesh))

## combine to get age-specific mean of rate
@mc.deterministic
def f(sm=sm, mu=mu):
    return mu + sm.f(pl.arange(101))

## condition on rate being positive
@mc.potential
def positive(f=f):
    if pl.any(f < 0.):
        return -pl.inf
    else:
        return 0.

## likelihood of observed data, using normal model for simplicity
data = pl.array([[ 0., .1, .005],
                 [50., .5, .005]])
@mc.deterministic
def data_expected(f=f):
    return [f[data[0,0]], f[data[1,0]]],  # TODO: make this line more elegant
@mc.observed
def data_obs(data_expected=data_expected, value=data):
    return mc.normal_like(value[:,1],
                          data_expected,
                          value[:,2]**-2)

## generate good initial values with MAP fit
mc.MAP([gamma, data_obs]).fit(method='fmin_powell', verbose=1)

## combine all PyMC nodes in to an MCMC object
mcmc = mc.MCMC([sm, gamma, mu, f, positive, data_expected, data_obs])

#mcmc.use_step_method(mc.gp.GPParentAdaptiveMetropolis, gamma)
#mcmc.use_step_method(mc.AdaptiveMetropolis, [gamma, sm.f_eval])
mcmc.sample(10000, 5000, 10)

### @export 'plot-varying-smoothing'
mc.Matplot.plot(data_expected)
mc.Matplot.plot(mu)

pl.figure()
ages = pl.arange(0,101,1)

for f_n in f.trace():
    pl.plot(ages, f_n,
         '-', color='grey',
            zorder=-1)
pl.plot(ages, f.stats()['quantiles'][50],
        'k-', linewidth=3)
pl.errorbar(data[:,0], data[:,1], yerr=data[:,2]*1.96,
            fmt='gs',
            mec='white', mew=1, ms=10)

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
