### @export 'setup'

import sys
sys.path += ['..', '../..']

import pylab as pl
import pymc as mc
import pandas

import age_pattern
import expert_prior_model
reload(age_pattern)
import data_simulation

import graphics
import dismod3
reload(dismod3)

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f0', '#ffff33']

### @export 'initialize'
df = pandas.read_csv('ssas_mx.csv', index_col=None)

ages = pl.arange(101)
knots = [0, 15, 60, 100]
import scipy.interpolate
Y_true = pl.exp(scipy.interpolate.interp1d(knots, pl.log([1.2, .3, .6, 1.5]), kind='linear')(ages))

N = 50
tau = .1**-2
X = pl.array(mc.runiform(pl.arange(0., 100., 100./N), 100./N + pl.arange(0., 100., 100./N), size=N), dtype=int)
Y = mc.rnormal(Y_true[X], tau)

### @export 'initial-rates'
pl.figure(figsize=(17., 11), dpi=72)

dismod3.graphics.plot_data_bars(df, 'talk')
pl.semilogy([0], [.1], '-')

pl.title('All-cause mortality rate\nin 1990 for females\nin sub-Saharan Africa, Southern.', size=55)
pl.ylabel('Rate (Per PY)', size=48)
pl.xlabel('Age (Years)', size=48)

pl.subplots_adjust(.1, .175, .98, .7)
pl.axis([-5, 105, 2.e-4, .8])
pl.xticks(size=30)
pl.yticks(size=30)
pl.grid()
pl.show()
pl.savefig('/media/windows/t/ssas_mx.png')


### @export 'spline_fig'
pl.figure(figsize=(17., 11), dpi=72)

pl.ylabel('$\mu(a)$', size=48)
pl.xlabel('$a$', size=48)
pl.xticks(size=30)
pl.yticks(size=30)
pl.grid()

pl.subplots_adjust(.1, .175, .98, .875, .275)

pl.plot(ages, Y_true, 'k-', linewidth=10, label='Truth')
pl.axis([-5, 105, 0., 1.7])
pl.savefig('/media/windows/t/splines-fig-1.png')

pl.plot(X, Y, 'kx', ms=15, mew=8, label='Simulated Data')
pl.axis([-5, 105, 0., 1.7])
pl.savefig('/media/windows/t/splines-fig-2.png')

for params in [dict(label='Piecewise Constant', interpolation_method='zero', linestyle='steps-mid-', color=colors[0]),
               dict(label='Piecewise Linear', interpolation_method='linear', linestyle='-', color=colors[1]),
               dict(label='Cubic', interpolation_method='cubic', linestyle='-', color=colors[2])]:

    interpolation_method=params.pop('interpolation_method')
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf, interpolation_method=interpolation_method)
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages, vars['mu_age'].value, 'k', linewidth=10, **params)

    pl.axis([-5, 105, 0., 1.7])
    pl.savefig('/media/windows/t/splines-fig-%s.png'%interpolation_method)

### @export 'monotone-spline_fig'

pl.figure(figsize=(17., 11), dpi=72)
pl.subplots_adjust(.1, .175, .98, .875, .275)

pl.ylabel('$\mu(a)$', size=48)
pl.xlabel('$a$', size=48)
pl.xticks(size=30)
pl.yticks(size=30)
pl.grid()

pl.plot(X, Y, 'kx', ms=15, mew=8)

pl.axis([-5, 105, 0., 1.7])
pl.savefig('/media/windows/t/prior-fig-1.png')


for params in [dict(label='$\mu(a)$ unconstrained', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=0)), color=colors[0]),
               dict(label='$\mu(a)$ decreasing for $a \leq 50$', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=50)), color=colors[1]),
               dict(label='$\mu(a)$ increasing for $a \leq 50$', value=dict(increasing=dict(age_start=0, age_end=50), decreasing=dict(age_start=0, age_end=0)), color=colors[2]),
               ]:
    #
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.derivative_constraints('t',
                                                          params.pop('value'),
                                                          vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    #
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages, vars['mu_age'].value, 'w-', linewidth=10, **params)

    pl.axis([-5, 105, 0., 1.7])
    pl.savefig('/media/windows/t/prior-fig-%s.png'%params['color'])

pl.show()
