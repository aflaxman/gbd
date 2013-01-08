# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

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

# adapted from pymc.MCMC.py
def _calc_dic(self):
    """Calculates deviance information Criterion"""
    # Find mean deviance
    mean_deviance = np.mean(self.db.trace('deviance')(), axis=0)
    # Set values of all parameters to their mean
    for stochastic in self.stochastics:

        # Calculate mean of paramter
        try:
            mean_value = np.mean(self.db.trace(stochastic.__name__)(), axis=0)

            # Set current value to mean
            stochastic.value = mean_value

        except KeyError:
            print_("No trace available for %s. DIC value may not be valid." % stochastic.__name__)

    # Return twice deviance minus deviance at means
    return 2*mean_deviance - self.deviance



# <codecell>

### @export 'initialize'
df = pandas.read_csv('book/ssas_mx.csv', index_col=None)
df['age_end'] += 1

# <codecell>

### @export 'initial-rates'
pl.figure(**book_graphics.quarter_page_params)

graphics.plot_data_bars(df)
pl.semilogy([0], [.1], '-')

pl.ylabel('Rate (per 1)', fontsize='xx-large')
pl.yticks(fontsize='large')
pl.xlabel('Age (years)', fontsize='xx-large')
pl.xticks(fontsize='large')

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 2.e-4, .8])

min_mx = df['value'].min()
age_min_mx = df['age_start'][df['value'].argmin()]
max_mx = df['value'].max()
age_max_mx = df['age_start'][df['value'].argmax()]
pl.grid()
pl.subplots_adjust(bottom=.2)
pl.savefig('ssas-mx_female_1990.pdf')

# <codecell>

print 'They vary %.0f-fold between the minimum in the %d to %d year olds, and the maximum at oldest ages.' % (max_mx/min_mx, age_min_mx, age_min_mx+4)

# <codecell>

print '$m_{\all}$ is as low as %.0f per 10,000 PY at age %d, but rises to %.0f per 10,000 PY at age %d.' % (min_mx*10000, age_min_mx, max_mx*10000, age_max_mx)

# <codecell>

mc.np.random.seed(1234567)

ages = pl.arange(101)
knots = [0, 15, 60, 100]
import scipy.interpolate
Y_true = pl.exp(scipy.interpolate.interp1d(knots, pl.log([1.2, .3, .6, 1.5]), kind='linear')(ages))

N = 50
tau = .1**-2
X = pl.array(mc.runiform(pl.arange(0., 100., 100./N), 100./N + pl.arange(0., 100., 100./N), size=N), dtype=int)
Y = mc.rnormal(Y_true[X], tau)

# <codecell>

### @export 'spline_fig'

pl.figure(**book_graphics.three_quarter_page_params)
pl.subplot(2,1,1)
pl.plot(ages, Y_true, 'k:', label='Truth')
pl.plot(X, Y, 'kx', ms=4, mew=2, label='Simulated data')

results = pandas.DataFrame(pl.zeros((2,3)),
                           index=['Piecewise constant','Piecewise linear'], 
                           columns=['AIC','BIC','DIC'])
                           
for params in [dict(label='Piecewise constant', interpolation_method='zero', linestyle='steps-mid--'),
               dict(label='Piecewise linear', interpolation_method='linear', linestyle='-'),]:
               #dict(label='Cubic', interpolation_method='cubic', linestyle='-.')]:

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf, interpolation_method=params.pop('interpolation_method'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    mc.MAP(vars).fit(method='fmin_powell', verbose=0)
    #pl.plot(ages, vars['mu_age'].value, 'w', linewidth=3, **params)
    pl.plot(ages, vars['mu_age'].value, 'k', **params)
    
    res = mc.MAP(vars)
    res.fit(method='fmin_powell', verbose=0)
    results.ix[params['label'],'AIC'] = res.AIC
    results.ix[params['label'],'BIC'] = res.BIC
    #results.ix[params[label],'DIC'] = res.DIC

def decorate_figure():
    pl.legend(loc='upper center', bbox_to_anchor=(.5,-.45), fancybox=True, shadow=True)
    #pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'x-large'})
    pl.xticks(fontsize='x-large')
    pl.xlabel('$a$', fontsize='xx-large')
    pl.yticks([0., .5, 1., 1.5], fontsize='x-large')
    pl.ylabel('$h(a)$', rotation=0, fontsize='xx-large')
    
    pl.subplots_adjust(.1, .1, .98, .875, .275) # l b r t w h
    pl.axis([-5, 105, 0., 1.7])
    pl.grid()
    
decorate_figure()
pl.savefig('splines-fig.pdf')

# <codecell>

### @export 'spline_fig'
knots = range(0, 101, 5)
pl.figure(**book_graphics.three_quarter_page_params)
pl.subplot(2,1,1)
pl.plot(X, Y, 'kx', ms=4, mew=2, label='Simulated data')

for params in [dict(label='$\sigma = 10^{-1}$', smoothing=.1, marker='s'),
               dict(label='$\sigma = 10^{-2}$', smoothing=.01, marker='o'),
               dict(label='$\sigma = 10^{-3}$', smoothing=.001, marker='^')]:
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=params.pop('smoothing'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=0)
    #pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)



decorate_figure()
pl.savefig('smoothing-splines.pdf')

# <codecell>

### @export 'level_value-spline_fig'

knots = [0, 15, 60, 100]
#knots = range(0,101,10)
pl.figure(**book_graphics.three_quarter_page_params)
pl.subplot(2,1,1)
pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$h(a) =$ $0.1$', value=.1, marker='s'),
               dict(label='$h(a) =$ $0.5$', value=.5, marker='o'),
               dict(label='$h(a) =$ $1$', value=1., marker='^')]:

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.level_constraints('t',
                                                     dict(level_value=dict(age_before=15, age_after=101, value=params.
pop('value')),
                                                          level_bounds=dict(upper=pl.inf, lower=-pl.inf)),
                                                     vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', verbose=0)
    #pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
#pl.legend(loc='lower right', fancybox=True, shadow=True, prop={'size':'xx-large'}, title='For $a < 15$:')
pl.savefig('level_value-smoothing-splines.pdf')

# <codecell>

### @export 'level_bound-spline_fig'

pl.figure(**book_graphics.three_quarter_page_params)
pl.subplot(2,1,1)
pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$0.2 \leq$ $h(a) \leq$ $1.5$', value=1.5, marker='s'),
               dict(label='$0.2 \leq$ $h(a) \leq$ $1.0$', value=1.0, marker='o'),
               dict(label='$0.2 \leq$ $h(a) \leq$ $0.8$', value=0.8, marker='^')]:
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
    mc.MAP(vars).fit(method='fmin_powell', verbose=0)
    #pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)


decorate_figure()
pl.yticks([0, .5, .8, 1., 1.5])
pl.savefig('level_bound-smoothing-splines.pdf')

# <codecell>

### @export 'monotone-spline_fig'

pl.figure(**book_graphics.three_quarter_page_params)
pl.subplot(2,1,1)
pl.plot(X, Y, 'kx', ms=4, mew=2)

for params in [dict(label='$h(a)$ unconstrained', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=0)), marker='s'),
               dict(label='$h(a)$ decreasing for $a \leq$ $50$', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=50)), marker='o'),
               dict(label='$h(a)$ increasing for $a \leq$ $50$', value=dict(increasing=dict(age_start=0, age_end=50), decreasing=dict(age_start=0, age_end=0)), marker='^'),
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
    mc.MAP(vars).fit(method='fmin_powell', verbose=0)
    #pl.plot(ages[knots], vars['mu_age'].value, 'w-', linewidth=3, **params)
    pl.plot(ages[knots], vars['mu_age'].value[knots], 'k-', mec='w', **params)



decorate_figure()
pl.savefig('monotone-smoothing-splines.pdf')

# <codecell>


