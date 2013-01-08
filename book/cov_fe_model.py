""" Fixed Effect Covariate model example
"""

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import scipy.interpolate

import graphics
import book_graphics
reload(book_graphics)
import data_simulation

import dismod3
reload(dismod3)

pi_true = scipy.interpolate.interp1d([0, 20, 40, 60, 100], [.4, .425, .6, .5, .4])
beta_true = .3
delta_true = 50.
N = 30

# start with a simple model with N rows of data
model = data_simulation.simple_model(N)


# set covariate to 0/1 values randomly
model.input_data['x_cov'] = 1. * mc.rcategorical([.5, .5], size=N)

# record the true age-specific rates
model.ages = pl.arange(0, 101, 1)
model.pi_age_true = pi_true(model.ages)


# choose age groups randomly
age_width = pl.zeros(N)
age_mid = mc.runiform(age_width/2, 100-age_width/2, size=N)
age_start = pl.array(age_mid - age_width/2, dtype=int)
age_end = pl.array(age_mid + age_width/2, dtype=int)

model.input_data['age_start'] = age_start
model.input_data['age_end'] = age_end


# choose effective sample size uniformly at random
n = mc.runiform(100, 10000, size=N)
model.input_data['effective_sample_size'] = n


# find true rate, with covariate
p = model.pi_age_true[age_start] * pl.exp(model.input_data['x_cov']*beta_true)


# sample observed rate values from negative binomial distribution
model.input_data['true'] = p
model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true) / n

print model.input_data.drop(['standard_error', 'upper_ci', 'lower_ci'], axis=1)






# Create age-group model
## Spline model to represent age-specific rate
model.vars += dismod3.age_pattern.spline(name='sim', ages=model.ages,
                                         knots=pl.arange(0,101,20),
                                         smoothing=pl.inf,
                                         interpolation_method='linear')

## Midpoint model to represent age-group data
model.vars += dismod3.age_group.midpoint_approx(name='sim', ages=model.ages,
                                                mu_age=model.vars['mu_age'], 
                                                age_start=model.input_data['age_start'], age_end=model.input_data['age_end'])

## FE model for sim covariate
model.vars += dismod3.covariates.mean_covariate_model(name='sim',
                                                     mu=model.vars['mu_interval'],
                                                     input_data=model.input_data,
                                                     parameters=model.parameters,
                                                     model=model,
                                                     root_area='all',
                                                     root_sex='total',
                                                     root_year='all')

## Uniform prior on negative binomial rate model over-dispersion term
model.vars += {'delta': mc.Uniform('delta_sim', 0., 1000., value=10.)}
model.vars += {'mu_age_1': mc.Lambda('mu_age_1_sim', lambda beta_0=model.vars['beta'][0], mu_age=model.vars['mu_age']: pl.exp(beta_0)*mu_age)}

## Negative binomial rate model
model.vars += dismod3.rate_model.neg_binom(name='sim',
                                           pi=model.vars['pi'],
                                           delta=model.vars['delta'],
                                           p=model.input_data['value'],
                                           n=model.input_data['effective_sample_size'])



#print "Fitting vars:"
#print model.vars.describe()
model.map = mc.MAP(model.vars)
model.map.fit(method='fmin_powell', verbose=False)

model.mcmc = mc.MCMC(model.vars)
model.mcmc.use_step_method(mc.AdaptiveMetropolis, model.vars['gamma'])
model.mcmc.sample(20000, 10000, 100, verbose=False, progress_bar=False)

# Always check model convergence
#mc.Matplot.plot(model.mcmc)
dismod3.graphics.plot_acorr(model.vars)
dismod3.graphics.plot_trace(model.vars)





pl.figure(**book_graphics.half_page_params)

for j in range(2):
    df = model.input_data[model.input_data['x_cov'] == j]
    pl.plot(df['age_start'].__array__(), df['value'].__array__(),
            color='k', linestyle='none', marker='ox'[j], mew=1+j, mec='wk'[j], ms=5,
            label='Observed ($x_i =$ $%d$)'%j)

pl.plot(model.ages[::10], model.pi_age_true[::10], 'w-', linewidth=3)
pl.plot(model.ages[::10], model.pi_age_true[::10], 'sk--', label='Truth ($x=$ $0$)')
pl.plot(model.ages[::10], model.pi_age_true[::10]*pl.exp(beta_true), 'w-', linewidth=3)
pl.plot(model.ages[::10], model.pi_age_true[::10]*pl.exp(beta_true), '^k--', label='Truth ($x=$ $1$)')


pl.plot(model.ages[::10], model.vars['mu_age'].stats()['mean'][::10], 'w-', linewidth=3)
pl.plot(model.ages[::10], model.vars['mu_age'].stats()['mean'][::10], 'sk-', label='Posterior mean ($x=$ $0$)')
pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,0], 'k:')
pl.plot(model.ages, model.vars['mu_age'].stats()['95% HPD interval'][:,1], 'k:')

pl.plot(model.ages[::10], model.vars['mu_age_1'].stats()['mean'][::10], 'w-', linewidth=3)
pl.plot(model.ages[::10], model.vars['mu_age_1'].stats()['mean'][::10], '^k-', label='Posterior mean ($x=$ $1$)')
pl.plot(model.ages, model.vars['mu_age_1'].stats()['95% HPD interval'][:,0], 'k:')
pl.plot(model.ages, model.vars['mu_age_1'].stats()['95% HPD interval'][:,1], 'k:')


pl.legend(fancybox=True, shadow=True, loc='lower center', ncol=3)
pl.xlabel('Age (years)')
pl.ylabel('Rate (per PY)')
pl.axis([-5, 105, 0., 1.])

pl.subplots_adjust(.1, .1, .98, .98, .275, 0)
pl.savefig('cov_fe.pdf')


pl.show()
