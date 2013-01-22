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

# <codecell>

book_graphics.set_font()

### @export 'initialize'
df = pandas.read_csv('/home/j/Project/dismod/gbd/data/ssas_mx.csv', index_col=None)
df['age_end'] += 1

# <codecell>

### @export 'initial-rates'
pl.figure(**book_graphics.quarter_page_params)

graphics.plot_data_bars(df)
pl.semilogy([0], [.1], '-')

pl.ylabel('Rate (Per 1)')#
pl.yticks(fontsize='large')
pl.xlabel('Age (Years)')
pl.xticks()

pl.subplots_adjust(.1, .175, .98, .875, .275)
pl.axis([-5, 105, 2.e-4, .8])

min_mx = df['value'].min()
age_min_mx = df['age_start'][df['value'].argmin()]
max_mx = df['value'].max()
age_max_mx = df['age_start'][df['value'].argmax()]

pl.savefig('book/graphics/ssas-mx_female_1990.pdf')

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

def decorate_figure():
    pl.legend(loc='lower right', fancybox=True, shadow=True)#, prop={'size':'x-large'})
    #pl.xticks(fontsize='x-large')
    pl.yticks([0., .5, 1., 1.5])#, fontsize='x-large')
    pl.ylabel('$h(a)$', rotation=0)#, fontsize='xx-large')
    
    pl.subplots_adjust(.1, .175, .98, .875, .275)
    pl.axis([-5, 105, 0., 1.7])
    

# <codecell>

fig = pl.figure(**book_graphics.three_quarter_page_params)

for i, params in enumerate([dict(label='Piecewise Constant', subt='(a)', interpolation_method='zero', linestyle='steps-mid-'),
                            dict(label='Piecewise Linear', subt='(b)', interpolation_method='linear', linestyle='-'),]):

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf, interpolation_method=params.pop('interpolation_method'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    mc.MAP(vars).fit(method='fmin_powell', tol=.00001, verbose=0)
    #pl.plot(ages, vars['mu_age'].value, 'w', linewidth=3, **params)
    
    if i == 0:
        ax1 = fig.add_subplot(2,1,i+1)
        ax1.plot(X, Y, 'ks', ms=4, mew=2)
        ax1.plot(ages, vars['mu_age'].value, 'k', linewidth=2, linestyle=params['linestyle'])#, **params)
    else:
        ax2 = fig.add_subplot(2,1,i+1, sharex=ax1)
        ax2.plot(X, Y, 'ks', ms=4, mew=2)
        ax2.plot(ages, vars['mu_age'].value, 'k', linewidth=2, linestyle=params['linestyle'])#, **params)
        pl.xlabel('$a$')
    decorate_figure()
    pl.setp(ax1.get_xticklabels(), visible=False)
    book_graphics.subtitle(params['subt'] + ' ' + params['label'])
pl.subplots_adjust(hspace=.1, top=.99, bottom=.11)
pl.savefig('book/graphics/splines-fig.pdf')

# <codecell>

### @export 'spline_fig'
knots = range(0, 101, 5)
fig1_data = {}

for i, params in enumerate([dict(label=r'$\sigma = 0.5$', subt='(a)', smoothing=.5),
                            dict(label=r'$\sigma = 0.05$', subt='(b)', smoothing=.05),
                            dict(label=r'$\sigma = 0.005$', subt='(c)', smoothing=.005)]):
    
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=params.pop('smoothing'))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', tol=.00001, verbose=0)
    
    fig1_data[params['subt']] = vars['mu_age'].value[knots]

fig1 = pl.figure(**book_graphics.full_minus_page_params)

ax11 = fig1.add_subplot(3,1,1)
ax11.plot(X, Y, 'ks', ms=4, mew=2)
ax11.plot(ages[knots], fig1_data['(a)'], 'k-', linewidth=2)
decorate_figure()
book_graphics.subtitle('(a) ' + r'$\sigma = 0.5$' )

ax12 = fig1.add_subplot(3,1,2, sharex=ax11)
ax12.plot(X, Y, 'ks', ms=4, mew=2)
ax12.plot(ages[knots], fig1_data['(b)'], 'k-', linewidth=2)
decorate_figure()
book_graphics.subtitle('(b) ' + r'$\sigma = 0.05$')

ax13 = fig1.add_subplot(3,1,3, sharex=ax11)
ax13.plot(X, Y, 'ks', ms=4, mew=2)
ax13.plot(ages[knots], fig1_data['(c)'], 'k-', linewidth=2)
decorate_figure()
book_graphics.subtitle('(c) ' + r'$\sigma = 0.005$')
pl.xlabel('$a$')

pl.setp(ax11.get_xticklabels(), visible=False)
pl.setp(ax12.get_xticklabels(), visible=False) 

pl.subplots_adjust(hspace=.1, top=.97, bottom=.08)

pl.savefig('book/graphics/smoothing-splines.pdf')

# <codecell>

### @export 'level_value-spline_fig'

knots = [0, 15, 60, 100]
#knots = range(0,101,10)

fig2_data = {}
for i, params in enumerate([dict(label='$h(a) = .1$ for $a<15$', subt='(a)', value=.1),
                            dict(label='$h(a) = .5$ for $a<15$', subt='(b)', value=.5),
                            dict(label='$h(a) = 1$ for $a<15$', subt='(c)', value=1.)]):
    
    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.level_constraints('t',
                                                     dict(level_value=dict(age_before=15, age_after=101, value=params.pop('value')),
                                                          level_bounds=dict(upper=pl.inf, lower=-pl.inf)),
                                                     vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)

    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', tol=.00001, verbose=0)
    
    fig2_data[params['subt']] = vars['mu_age'].value[knots]
    
fig2 = pl.figure(**book_graphics.full_minus_page_params)

ax21 = fig2.add_subplot(3,1,1)
ax21.plot(X, Y, 'ks', ms=4, mew=2)
ax21.plot(ages[knots], fig2_data['(a)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(a) ' + '$h(a) = 0.1$ for $a<15$')

ax22 = fig2.add_subplot(3,1,2,sharex=ax21)
ax22.plot(X, Y, 'ks', ms=4, mew=2)
ax22.plot(ages[knots], fig2_data['(b)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(b) ' + '$h(a) = 0.5$ for $a<15$')

ax23 = fig2.add_subplot(3,1,3,sharex=ax21)
ax23.plot(X, Y, 'ks', ms=4, mew=2)
ax23.plot(ages[knots], fig2_data['(c)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(c) ' + '$h(a) = 1$ for $a<15$')
pl.xlabel('$a$')

pl.setp(ax21.get_xticklabels(), visible=False)
pl.setp(ax22.get_xticklabels(), visible=False) 

pl.subplots_adjust(hspace=.1, top=.97, bottom=.08)

pl.savefig('book/graphics/level_value-smoothing-splines.pdf')

# <codecell>

i

# <codecell>

### @export 'level_bound-spline_fig'

fig3_data = {}
for i, params in enumerate([dict(label='$.2 \leq h(a) \leq 1.5$', subt='(a)', value=1.5),
                            dict(label='$.2 \leq h(a) \leq 1.0$', subt='(b)', value=1.0),
                            dict(label='$.2 \leq h(a) \leq 0.8$', subt='(c)', value=0.8)]):

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
    mc.MAP(vars).fit(method='fmin_powell', tol=.00001, verbose=0)

    fig3_data[params['subt']] = vars['mu_age'].value[knots]
    
fig3 = pl.figure(**book_graphics.full_minus_page_params)

ax31 = fig3.add_subplot(3,1,1)
ax31.plot(X, Y, 'ks', ms=4, mew=2)
ax31.plot(ages[knots], fig3_data['(a)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([0, .5, 1., 1.5])
pl.hlines([1.5], -5, 105, linewidth=1, linestyle='dotted', color='k')
book_graphics.subtitle('(a) ' + '$.2 \leq h(a) \leq 1.5$')

ax32 = fig3.add_subplot(3,1,2,sharex=ax31)
ax32.plot(X, Y, 'ks', ms=4, mew=2)
ax32.plot(ages[knots], fig3_data['(b)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([0, .5, 1., 1.5])
pl.hlines([1.], -5, 105, linewidth=1, linestyle='dotted', color='k')
book_graphics.subtitle('(b) ' + '$.2 \leq h(a) \leq 1.0$')

ax33 = fig3.add_subplot(3,1,3,sharex=ax31)
ax33.plot(X, Y, 'ks', ms=4, mew=2)
ax33.plot(ages[knots], fig3_data['(c)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([0, .5, 1., 1.5])
pl.hlines([.8], -5, 105, linewidth=1, linestyle='dotted', color='k')
book_graphics.subtitle('(c) ' + '$.2 \leq h(a) \leq 0.8$')
pl.xlabel('$a$')

pl.setp(ax31.get_xticklabels(), visible=False)
pl.setp(ax32.get_xticklabels(), visible=False) 

pl.subplots_adjust(hspace=.1, top=.97, bottom=.08)
pl.savefig('book/graphics/level_bound-smoothing-splines.pdf')

# <codecell>

### @export 'monotone-spline_fig'

fig4_data = {}
for i, params in enumerate([dict(label='$h(a)$ unconstrained', subt='(a)', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=0))),
                            dict(label='$h(a)$ decreasing for $a \leq 50$', subt='(b)', value=dict(increasing=dict(age_start=0, age_end=0), decreasing=dict(age_start=0, age_end=50))),
                            dict(label='$h(a)$ increasing for $a \leq 50$', subt='(c)', value=dict(increasing=dict(age_start=0, age_end=50), decreasing=dict(age_start=0, age_end=0)))]):

    vars = age_pattern.age_pattern('t', ages=ages, knots=knots, smoothing=pl.inf)
    vars.update(expert_prior_model.derivative_constraints('t',
                                                          params.pop('value'),
                                                          vars['mu_age'], ages))
    vars['mu_pred'] = mc.Lambda('mu_pred', lambda mu_age=vars['mu_age'], X=X : mu_age[X])
    vars['Y'] = mc.Normal('Y', mu=vars['mu_pred'], tau=tau, value=Y, observed=True)
    #
    for i, k_i in enumerate(knots):
        vars['gamma'][i].value = Y_true[k_i]
    mc.MAP(vars).fit(method='fmin_powell', tol=.00001, verbose=0)
    
    fig4_data[params['subt']] = vars['mu_age'].value[knots]

fig4 = pl.figure(**book_graphics.full_minus_page_params)

ax41 = fig4.add_subplot(3,1,1)
ax41.plot(X, Y, 'ks', ms=4, mew=2)
ax41.plot(ages[knots], fig4_data['(a)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(a) ' + '$h(a)$ unconstrained')

ax42 = fig4.add_subplot(3,1,2,sharex=ax41)
ax42.plot(X, Y, 'ks', ms=4, mew=2)
ax42.plot(ages[knots], fig4_data['(b)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(b) ' + '$h(a)$ decreasing for $a \leq 50$')

ax43 = fig4.add_subplot(3,1,3,sharex=ax41)
ax43.plot(X, Y, 'ks', ms=4, mew=2)
ax43.plot(ages[knots], fig4_data['(c)'], 'k-', linewidth=2)
decorate_figure()
pl.yticks([.1, .5, 1., 1.5])
book_graphics.subtitle('(c) ' + '$h(a)$ increasing for $a \leq 50$')
pl.xlabel('$a$')

pl.setp(ax41.get_xticklabels(), visible=False)
pl.setp(ax42.get_xticklabels(), visible=False) 

pl.subplots_adjust(hspace=.1, top=.97, bottom=.08)

pl.savefig('book/graphics/monotone-smoothing-splines.pdf')

# <codecell>

pl.show()
