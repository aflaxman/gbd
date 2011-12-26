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
    a = pl.arange(0, 100, 100./N)
    age_mid = mc.runiform(a, a+25, size=N)
    age_width = mc.runiform(0, 100, size=N)

    age_start = pl.array(age_mid - age_width/2, dtype=int).clip(0, 100)
    age_end = pl.array(age_mid + age_width/2, dtype=int).clip(0, 100)

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

    print model.input_data.drop(['standard_error', 'upper_ci', 'lower_ci'],1)
    return model


def fit_midpoint_model(model):
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


def fit_age_standardizing_model(model):
    model.parameters['p']['parameter_age_mesh'] = range(0, 101, 10)
    model.parameters['p']['smoothness'] = dict(amount='Slightly')

    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=25, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)
    graphics.plot_one_type(model, model.vars['p'], {}, 'p')
    pl.plot(model.ages, model.pi_age_true, 'r:', label='Truth')
    pl.legend(fancybox=True, shadow=True, loc='upper left')

    pl.show()

    return model


def fit_model(model):
    """ Fit model with MCMC, starting from MAP as initial value, and
    plot results"""

    print "Fitting vars:"
    print model.vars.describe()
    model.map = mc.MAP(model.vars)
    model.map.fit(method='fmin_powell', verbose=1)

    model.mcmc = mc.MCMC(model.vars)
    model.mcmc.use_step_method(mc.AdaptiveMetropolis, model.vars['gamma'])
    model.mcmc.sample(20000, 10000, 10)

    # Always check model convergence
    #mc.Matplot.plot(model.mcmc)
    model.vars.plot_acorr()
    model.vars.plot_trace()
    

    #model.vars.plot_age_pattern(model.input_data)  # TODO: determine if this should method of ModelVars or aux function in dismod3.graphics submodule
    graphics.plot_one_type(model, model.vars, {}, 'p')
    pl.plot(model.ages, model.pi_age_true, 'k--', label='Truth')
    pl.legend(fancybox=True, shadow=True, loc='upper left')

    pl.show()



if __name__ == '__main__':
    model = simulate_age_group_data()
    #fit_age_standardizing_model(model)
    fit_midpoint_model(model)
