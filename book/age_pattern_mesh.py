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
in_mesh = dismod3.settings.gbd_ages
out_mesh = pl.arange(101)
data = pl.array([[10., 1, .25],
                 [50., 2.5, .25]])

scale = dict(Very=dismod3.utils.rho['very'],
             Moderately=dismod3.utils.rho['moderately'],
             Slightly=dismod3.utils.rho['slightly'])

for kind in ['linear', 'zero']:
    pl.figure(**book_graphics.quarter_page_params)
    for col, smoothness in enumerate(['Slightly', 'Moderately', 'Very']):

        ## setup lognormal prior on intercept
        gamma = mc.Normal('gamma', 0.,
                          2.**-2,
                          value=pl.log(data[:,1].mean()))
        mu = mc.Lambda('mu', lambda gamma=gamma: pl.exp(gamma))

        ## setup Gaussian Process prior on age pattern
        @mc.deterministic
        def M(mu=mu):
            return mc.gp.Mean(lambda x: mu*pl.ones(len(x)))
        C = mc.gp.FullRankCovariance(mc.gp.matern.euclidean,
                                     amp=data[:,1].max(),
                                     scale=scale[smoothness],
                                     diff_degree=2)
        sm = mc.gp.GPSubmodel('sm', M, C, in_mesh,
                              init_vals=mu.value*pl.ones_like(in_mesh))

        @mc.deterministic
        def f(f=sm.f_eval):
            return dismod3.utils.interpolate(in_mesh, f, out_mesh, kind=kind)

        ## condition on rate being positive
        @mc.potential
        def positive(f=f):
            if pl.any(f < 0.):
                return -pl.inf
            else:
                return 0.

        ## likelihood of observed data, using normal model for simplicity
        @mc.deterministic
        def data_expected(f=f):
            return [f[data[0,0]], f[data[1,0]]]
        @mc.observed
        def data_obs(data_expected=data_expected, value=data):
            return mc.normal_like(value[:,1],
                                  data_expected,
                                  value[:,2]**-2)

        ## generate good initial values with MAP fit
        mc.MAP([gamma, data_obs]).fit(method='fmin_powell', verbose=1)
        mc.MAP([sm.f_eval, data_obs]).fit(method='fmin_powell', verbose=1)
        mc.MAP([gamma, data_obs]).fit(method='fmin_powell', verbose=1)

        ## sample from posterior distribution with MCMC
        mcmc = mc.MCMC([gamma, mu, M, sm, f, positive, data_expected, data_obs])
        mcmc.use_step_method(mc.gp.GPParentAdaptiveMetropolis, gamma)
        mcmc.sample(5000, 4000, 10)

        ### @export 'plot-varying-smoothing'
        pl.subplot(1, 3, col+1)
        if kind == 'zero':
            my_plot = pl.step
        else:
            my_plot = pl.plot
        for n, f_n in enumerate(f.trace()):
            my_plot(out_mesh, f_n,
                 '-', color='grey', linewidth=1,
                    zorder=-1)
            if n % 250 == 0:
                my_plot(out_mesh, f_n,
                        '-', color='white', linewidth=3,
                        zorder=1)
                my_plot(out_mesh, f_n,
                        '-', color='black', linewidth=1,
                        zorder=2)

        my_plot(out_mesh, f.stats()['95% HPD interval'],
                'k--', linewidth=1, label='Uncertainty interval')
        pl.errorbar(data[:,0], data[:,1], yerr=data[:,2]*1.96,
                    fmt='gs',
                    mec='white', mew=0, ms=10)

        pl.axis([0, 100, 0, 10])
        pl.text(10, 9, '%s ($\\rho = %d$)' % (smoothness, scale[smoothness]), ha='left', va='top')
        if col == 0:
            pl.ylabel('Rate (per 1000 PY)')
        else:
            pl.yticks([])
        pl.xticks([25,50,75])
        pl.xlabel('Age (Years)')

    pl.subplots_adjust(left=.1, bottom=.2, top=.95, right=.95, wspace=0)
    pl.savefig('smoothness_%s_priors.png'%kind)
    pl.savefig('smoothness_%s_priors.pdf'%kind)

