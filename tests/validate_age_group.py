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
    return pl.exp(-.003*(a - 35)**2.  + .01*(a - 35))


def validate_age_group(model, replicate):
    # set random seed for reproducibility
    mc.np.random.seed(1234567+replicate)

    N = 30
    delta_true = 5.
    pi_true = true_rate_function
    m = simulate_age_group_data(N=N, delta_true=delta_true, pi_true=pi_true)
    
    if model == 'midpoint_covariate':
        fit_midpoint_covariate_model(m)
    elif model == 'age_standardizing':
        fit_age_standardizing_model(m)
    elif model == 'midpoint_model':
        fit_midpoint_model(m)
    elif model == 'disaggregation_model':
        fit_disaggregation_model(m)
    else:
        raise TypeError, 'Unknown model type: "%s"' % model


    # compare estimate to ground truth
    import data_simulation
    m.mu = pandas.DataFrame(dict(true=[pi_true(a) for a in range(101)],
                                 mu_pred=m.vars['mu_age'].stats()['mean'],
                                 sigma_pred=m.vars['mu_age'].stats()['standard deviation']))
    data_simulation.add_quality_metrics(m.mu)
    print '\nparam prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (m.mu['abs_err'].mean(),
                                                                         pl.median(pl.absolute(m.mu['rel_err'].dropna())),
                                                                         m.mu['covered?'].mean())
    print


    data_simulation.add_quality_metrics(m.mu)

    data_simulation.initialize_results(m)
    data_simulation.add_to_results(m, 'mu')
    data_simulation.finalize_results(m)

    return m


def simulate_age_group_data(N=50, delta_true=150, pi_true=true_rate_function):
    """ generate simulated data
    """
    # start with a simple model with N rows of data
    model = data_simulation.simple_model(N)


    # record the true age-specific rates
    model.ages = pl.arange(0, 101, 1)
    model.pi_age_true = pi_true(model.ages)


    # choose age groups randomly
    age_width = mc.runiform(1, 100, size=N)
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


# TODO: move all these models into the age_group_model.py and expose them as an option in the asr creating methods
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

def fit_model(model):
    """ Fit model with MCMC, starting from MAP as initial value, and
    plot results"""

    print "Fitting vars:"
    print model.vars.describe()

    import time
    start_time = time.time()
    model.map = mc.MAP(model.vars)
    model.map.fit(method='fmin_powell', verbose=1)

    model.mcmc = mc.MCMC(model.vars)
    model.mcmc.use_step_method(mc.AdaptiveMetropolis, model.vars['gamma'])
    model.mcmc.sample(20000, 10000, 100)
    model.mcmc.wall_time = time.time() - start_time
