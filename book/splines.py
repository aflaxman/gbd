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
import book_graphics
reload(book_graphics)


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
pl.figure(**book_graphics.quarter_page_params)

graphics.plot_data_bars(df)
pl.semilogy([0], [.1], '-')

pl.ylabel('Rate (Per PY)')
pl.xlabel('Age (Years)')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 2.e-4, .8])

min_mx = df['value'].min()
age_min_mx = df['age_start'][df['value'].argmin()]
max_mx = df['value'].max()
age_max_mx = df['age_start'][df['value'].argmax()]

pl.savefig('ssas-mx_female_1990.pdf')


### @export 'spline_fig'

pl.figure(**book_graphics.half_page_params)

pl.plot(ages, Y_true, 'k-', label='Truth')
pl.plot(X, Y, 'kx', ms=4, mew=2, label='Simulated Data')

for params in [dict(label='Piecewise Constant', interpolation_method='zero', linestyle='steps-mid--'),
               dict(label='Piecewise Linear', interpolation_method='linear', linestyle='-.'),
               dict(label='Cubic', interpolation_method='cubic', linestyle=':')]:

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf, interpolation_method=params.pop('interpolation_method'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages, vars['mu_age'].value, 'w', linewidth=3, **params)
    pl.plot(ages, vars['mu_age'].value, 'k', **params)


pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'small'})
pl.ylabel('$\mu(a)$')
pl.xlabel('$a$')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 0., 1.7])
pl.savefig('splines-fig.pdf')

### @export 'spline_fig'

pl.figure(**book_graphics.half_page_params)

pl.plot(X, Y, 'kx', ms=4, mew=2, label='Simulated Data')

for params in [dict(label='$\sigma = 10^{-1}$', smoothing=.1, marker='s'),
               dict(label='$\sigma = 10^{-2}$', smoothing=.01, marker='o'),
               dict(label='$\sigma = 10^{-3}$', smoothing=.001, marker='^')]:
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=params.pop('smoothing'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'medium'})
pl.ylabel('$\mu(a)$')
pl.xlabel('$a$')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 0., 1.7])
pl.savefig('smoothing-splines.pdf')


### @export 'level_value-spline_fig'

pl.figure(**book_graphics.half_page_params)

pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$\mu(a) = .1$', value=.1, marker='s'),
               dict(label='$\mu(a) = .5$', value=.5, marker='o'),
               dict(label='$\mu(a) = 1$', value=1., marker='^')]:

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.level_constraints('t',
                                                     dict(level_value=dict(age_before=15, age_after=101, value=params.pop('value')),
                                                          level_bounds=dict(upper=pl.inf, lower=-pl.inf)),
                                                     vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'medium'}, title='For $a < 15$:')
pl.ylabel('$\mu(a)$')
pl.xlabel('$a$')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 0., 1.7])
pl.savefig('level_value-smoothing-splines.pdf')


### @export 'level_bound-spline_fig'

pl.figure(**book_graphics.half_page_params)

pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$.2 \leq \mu(a) \leq 1.5$', value=1.5, marker='s'),
               dict(label='$.2 \leq \mu(a) \leq 1.0$', value=1.0, marker='o'),
               dict(label='$.2 \leq \mu(a) \leq 0.8$', value=0.8, marker='^')]:
    #
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.level_constraints('t',
                                                     dict(level_value=dict(age_before=0, age_after=101, value=0.),
                                                          level_bounds=dict(upper=params.pop('value'), lower=.2)),
                                                     vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    #
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=1)
    pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'medium'})
pl.ylabel('$\mu(a)$')
pl.xlabel('$a$')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 0., 1.7])
pl.savefig('level_bound-smoothing-splines.pdf')



### @export 'monotone-spline_fig'

pl.figure(**book_graphics.half_page_params)

pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$\mu(a)$ unconstrained', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=0)), marker='s'),
               dict(label='$\mu(a)$ decreasing for $a \leq 50$', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=50)), marker='o'),
               dict(label='$\mu(a)$ increasing for $a \leq 50$', value=dict(increasing=dict(age_start=0, age_end=50), decreasing=dict(age_start=0, age_end=0)), marker='^'),
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
    pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'medium'})
pl.ylabel('$\mu(a)$')
pl.xlabel('$a$')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 0., 1.7])
pl.savefig('monotone-smoothing-splines.pdf')

