""" Compare approaches to age group modeling
"""

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc

import data
import data_model
import fit_model
import graphics
import data_simulation
import pandas

import dismod3
reload(dismod3)

import book_graphics

reload(data)
reload(data_simulation)
reload(graphics)

def true_rate_function(a):
    return pl.exp(-.001*(a - 35)**2.  + .01*(a - 35))


def simulate_age_group_data(N=50, delta_true=150, pi_true=true_rate_function):
    """ generate simulated data
    """
    # start with a simple model with N rows of data
    model = data_simulation.simple_model(N)


    # record the true age-specific rates
    model.ages = pl.arange(0, 101, 1)
    model.pi_age_true = pi_true(model.ages)


    # choose age groups randomly
    age_width = pl.ones(N)*50
    age_mid = mc.runiform(age_width/2, 100-age_width/2, size=N)
    age_width[:10] = 10
    age_mid[:10] = pl.arange(5, 105, 10)
    #age_width[10:20] = 10
    #age_mid[10:20] = pl.arange(5, 105, 10)

    age_start = pl.array(age_mid - age_width/2, dtype=int)
    age_end = pl.array(age_mid + age_width/2, dtype=int)

    model.input_data['age_start'] = age_start
    model.input_data['age_end'] = age_end


    # choose effective sample size uniformly at random
    n = mc.runiform(100, 10000, size=N)
    model.input_data['effective_sample_size'] = n


    # integrate true age-specific rate across age groups to find true group rate
    age_weights = pl.ones_like(model.ages)
    sum_pi_wt = pl.cumsum(model.pi_age_true*age_weights)
    sum_wt = pl.cumsum(age_weights)
    p = (sum_pi_wt[age_end] - sum_pi_wt[age_start]) / (sum_wt[age_end] - sum_wt[age_start])


    # sample observed rate values from negative binomial distribution
    model.input_data['true'] = p
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true) / n

    print model.input_data.drop(['standard_error', 'upper_ci', 'lower_ci'], axis=1)
    return model


def fit_midpoint_model(model):
    """Midpoint model"""
    # Create age-group model
    ## Spline model to represent age-specific rate
    model.vars += dismod3.age_pattern.spline(name='mid', ages=model.ages,
                                             knots=pl.arange(0,101,20),
                                             smoothing=pl.inf,
                                             interpolation_method='linear')

    ## Midpoint model to represent age-group data
    model.vars += dismod3.age_group.midpoint_approx(name='mid', ages=model.ages,
                                                    mu_age=model.vars['mu_age'], 
                                                    age_start=model.input_data['age_start'], age_end=model.input_data['age_end'])

    ## Uniform prior on negative binomial rate model over-dispersion term
    model.vars += {'delta': mc.Uniform('delta_mid', 0., 1000., value=10.)}

    ## Negative binomial rate model
    model.vars += dismod3.rate_model.neg_binom(name='mid',
                                               pi=model.vars['mu_interval'],
                                               delta=model.vars['delta'],
                                               p=model.input_data['value'],
                                               n=model.input_data['effective_sample_size'])

    fit_model(model)


def fit_midpoint_covariate_model(model):
    """Midpoint/covariate model"""
    # Create age-group model
    ## Spline model to represent age-specific rate
    model.vars += dismod3.age_pattern.spline(name='midc', ages=model.ages,
                                             knots=pl.arange(0,101,20),
                                             smoothing=pl.inf,
                                             interpolation_method='linear')

    ## Midpoint model to represent age-group data
    model.vars += dismod3.age_group.midpoint_covariate_approx(name='midc', ages=model.ages,
                                                              mu_age=model.vars['mu_age'], 
                                                              age_start=model.input_data['age_start'], age_end=model.input_data['age_end'])

    ## Uniform prior on negative binomial rate model over-dispersion term
    model.vars += {'delta': mc.Uniform('delta_midc', 0., 1000., value=10.)}

    ## Negative binomial rate model
    model.vars += dismod3.rate_model.neg_binom(name='midc',
                                               pi=model.vars['mu_interval'],
                                               delta=model.vars['delta'],
                                               p=model.input_data['value'],
                                               n=model.input_data['effective_sample_size'])

    fit_model(model)


def fit_disaggregation_model(model):
    """Disaggregation approach"""
    ## Spline model to represent age-specific rate
    model.vars += dismod3.age_pattern.spline(name='dis', ages=model.ages,
                                             knots=pl.arange(0,101,20),
                                             smoothing=pl.inf,
                                             interpolation_method='linear')

    ## Disaggregate input data
    a = []
    p = []
    n = []
    for i in model.input_data.index:
        a_s, a_e = model.input_data.ix[i, 'age_start'], model.input_data.ix[i, 'age_end']
        a += range(a_s, a_e)
        p += [model.input_data.ix[i, 'value']] * (a_e - a_s)
        n += [float(model.input_data.ix[i, 'effective_sample_size']) / (a_e - a_s)] * (a_e - a_s)
    a = pl.array(a)
    p = pl.array(p)
    n = pl.array(n)

    model.vars['pi'] = mc.Lambda('pi_dis', lambda mu=model.vars['mu_age'], a=a: mu[a])
    model.vars['delta'] = mc.Uniform('delta_dis', 0., 1000., value=10.)


    ## Negative binomial rate model
    model.vars += dismod3.rate_model.neg_binom(name='dis',
                                               pi=model.vars['pi'],
                                               delta=model.vars['delta'],
                                               p=p,  # TODO: change this parameter name to "r" to match the book chapter
                                               n=n)

    fit_model(model)


def fit_age_standardizing_model(model):
    """Age-standardizing model"""
    ## Spline model to represent age-specific rate
    model.vars += dismod3.age_pattern.spline(name='std', ages=model.ages,
                                             knots=pl.arange(0,101,20),
                                             smoothing=pl.inf,
                                             interpolation_method='linear')

    model.vars += dismod3.age_group.age_standardize_approx(name='std', ages=model.ages,
                                                           age_weights=pl.ones_like(model.ages),
                                                           mu_age=model.vars['mu_age'], 
                                                           age_start=model.input_data['age_start'], age_end=model.input_data['age_end'])

    model.vars += {'delta': mc.Uniform('delta_std', 0., 1000., value=10.)}

    ## Negative binomial rate model
    model.vars += dismod3.rate_model.neg_binom(name='std',
                                               pi=model.vars['mu_interval'],
                                               delta=model.vars['delta'],
                                               p=model.input_data['value'],
                                               n=model.input_data['effective_sample_size'])

    fit_model(model)

fit_age_standardizing_model.fmt = '^-1w7'
fit_midpoint_model.fmt = 'o-1w5'
fit_midpoint_covariate_model.fmt = 'x-2k5'
fit_disaggregation_model.fmt = '*-1w13'

def fit_model(model):
    """ Fit model with MCMC, starting from MAP as initial value, and
    plot results"""

    print "Fitting vars:"
    print model.vars.describe()
    model.map = mc.MAP(model.vars)
    model.map.fit(method='fmin_powell', verbose=1)

    model.mcmc = mc.MCMC(model.vars)
    model.mcmc.use_step_method(mc.AdaptiveMetropolis, model.vars['gamma'])
    model.mcmc.sample(20000, 10000, 100)

    # Always check model convergence
    #mc.Matplot.plot(model.mcmc)
    model.vars.plot_acorr()
    model.vars.plot_trace()
    

    #model.vars.plot_age_pattern(model.input_data)  # TODO: determine if this should method of ModelVars or aux function in dismod3.graphics submodule
    #graphics.plot_one_type(model, model.vars, {}, 'p')
    #pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
    #pl.legend(fancybox=True, shadow=True, loc='upper left')



if __name__ == '__main__':
    m = {}
    for fit in [fit_midpoint_covariate_model, fit_age_standardizing_model, fit_midpoint_model, fit_disaggregation_model]:
        m[fit] = simulate_age_group_data(N=30, delta_true=5)
        m[fit].input_data = model.input_data.copy()
        fit(m[fit])

    pl.figure(**book_graphics.half_page_params)
    graphics.plot_data_bars(model.input_data)
    pl.plot(model.ages, model.pi_age_true, 'w-', linewidth=3)
    pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
    for fit in [fit_age_standardizing_model, fit_midpoint_model, fit_midpoint_covariate_model, fit_disaggregation_model]:
        pl.plot(model.ages[::10], m[fit].vars['mu_age'].stats()['mean'][::10], 'w-', linewidth=3)
        pl.plot(model.ages[::10], m[fit].vars['mu_age'].stats()['mean'][::10], marker=fit.fmt[0], linestyle=fit.fmt[1],
                mew=float(fit.fmt[2]), mec=fit.fmt[3], ms=float(fit.fmt[4:]), color='k', label=fit.__doc__)
    pl.legend(fancybox=True, shadow=True, loc='upper right')
    pl.xlabel('Age (Years)')
    pl.ylabel('Rate (Per PY)')
    pl.axis([-5, 105, 0., 1.5])
    pl.subplots_adjust(.1, .175, .98, .875, .275)
    pl.savefig('age_group_models.pdf')

    pl.show()
