""" Validate similarity prior
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc

pl.seterr('ignore')

import data
import data_model
import fit_model
import graphics
import data_simulation
import pandas

reload(data_simulation)
reload(graphics)

def quadratic(a):
    return .0001 * (a * (100. - a) + 100.)

def validate_similarity(N=15, delta_true=.15, pi_true=quadratic, heterogeneity='Very', bias=0., sigma_prior=.5):
    model = generate_data(N, delta_true, pi_true, heterogeneity, bias, sigma_prior)
    fit(model)
    return model

def generate_data(N, delta_true, pi_true, heterogeneity, bias, sigma_prior):
    a = pl.arange(0, 101, 1)
    pi_age_true = pi_true(a)

    model = data_simulation.simple_model(N)
    model.parameters['p']['parameter_age_mesh'] = range(0, 101, 10)
    model.parameters['p']['smoothness'] = dict(amount='Moderately')
    model.parameters['p']['heterogeneity'] = heterogeneity

    age_start = pl.array(mc.runiform(0, 100, size=N), dtype=int)
    age_end = pl.array(mc.runiform(age_start, 100, size=N), dtype=int)

    age_weights = pl.ones_like(a)
    sum_pi_wt = pl.cumsum(pi_age_true*age_weights)
    sum_wt = pl.cumsum(age_weights)
    p = (sum_pi_wt[age_end] - sum_pi_wt[age_start]) / (sum_wt[age_end] - sum_wt[age_start])

    # correct cases where age_start == age_end
    i = age_start == age_end
    if pl.any(i):
        p[i] = pi_age_true[age_start[i]]

    n = mc.runiform(10000, 100000, size=N)

    model.input_data['age_start'] = age_start
    model.input_data['age_end'] = age_end
    model.input_data['effective_sample_size'] = n
    model.input_data['true'] = p
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true*n*p) / n * pl.exp(bias)

    emp_priors = {}
    emp_priors['p', 'mu'] = pi_age_true
    emp_priors['p', 'sigma'] = sigma_prior*pi_age_true
    model.emp_priors = emp_priors

    model.a = a
    model.pi_age_true = pi_age_true
    model.delta_true = delta_true

    return model

def fit(model):
    emp_priors = model.emp_priors

    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, emp_priors['p', 'mu'], emp_priors['p', 'sigma'])
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=5000, burn=2000, thin=25, tune_interval=100)
    #model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=101, burn=0, thin=1, tune_interval=100)

    #graphics.plot_one_ppc(model.vars['p'], 'p')
    #graphics.plot_convergence_diag(model.vars)
    graphics.plot_one_type(model, model.vars['p'], emp_priors, 'p')
    pl.plot(model.a, model.pi_age_true, 'b--', linewidth=3, alpha=.5, label='Truth')
    pl.legend(fancybox=True, shadow=True, loc='upper left')
    pl.title('Heterogeneity %s'%model.parameters['p']['heterogeneity'])

    pl.show()

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    data_simulation.add_quality_metrics(model.input_data)

    model.delta = pandas.DataFrame(dict(true=[model.delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    data_simulation.add_quality_metrics(model.delta)

    print 'delta'
    print model.delta

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())

    model.mu = pandas.DataFrame(dict(true=model.pi_age_true,
                                     mu_pred=model.vars['p']['mu_age'].stats()['mean'],
                                     sigma_pred=model.vars['p']['mu_age'].stats()['standard deviation']))
    data_simulation.add_quality_metrics(model.mu)

    data_simulation.initialize_results(model)
    data_simulation.add_to_results(model, 'delta')
    data_simulation.add_to_results(model, 'mu')
    data_simulation.add_to_results(model, 'input_data')
    data_simulation.finalize_results(model)

    print model.results


if __name__ == '__main__':
    model = validate_similarity()
