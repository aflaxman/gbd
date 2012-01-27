""" Validate Age Integrating Model
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

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

reload(data_simulation)
reload(graphics)

def quadratic(a):
    return .0001 * (a * (100. - a) + 100.)

def validate_age_integrating_model_sim(N=500, delta_true=.15, pi_true=quadratic):
    ## generate simulated data
    a = pl.arange(0, 101, 1)
    pi_age_true = pi_true(a)

    model = data_simulation.simple_model(N)
    #model.parameters['p']['parameter_age_mesh'] = range(0, 101, 10)
    #model.parameters['p']['smoothness'] = dict(amount='Very')

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

    n = mc.runiform(100, 10000, size=N)

    model.input_data['age_start'] = age_start
    model.input_data['age_end'] = age_end
    model.input_data['effective_sample_size'] = n
    model.input_data['true'] = p
    model.input_data['value'] = mc.rnegative_binomial(n*p, delta_true*n*p) / n

    ## Then fit the model and compare the estimates to the truth
    model.vars = {}
    model.vars['p'] = data_model.data_model('p', model, 'p', 'all', 'total', 'all', None, None, None)
    model.map, model.mcmc = fit_model.fit_data_model(model.vars['p'], iter=10000, burn=5000, thin=25, tune_interval=100)

    graphics.plot_one_ppc(model.vars['p'], 'p')
    graphics.plot_convergence_diag(model.vars)
    graphics.plot_one_type(model, model.vars['p'], {}, 'p')
    pl.plot(a, pi_age_true, 'r:', label='Truth')
    pl.legend(fancybox=True, shadow=True, loc='upper left')

    pl.show()

    model.input_data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    model.input_data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']
    data_simulation.add_quality_metrics(model.input_data)

    model.delta = pandas.DataFrame(dict(true=[delta_true]))
    model.delta['mu_pred'] = pl.exp(model.vars['p']['eta'].trace()).mean()
    model.delta['sigma_pred'] = pl.exp(model.vars['p']['eta'].trace()).std()
    data_simulation.add_quality_metrics(model.delta)

    print 'delta'
    print model.delta

    print '\ndata prediction bias: %.5f, MARE: %.3f, coverage: %.2f' % (model.input_data['abs_err'].mean(),
                                                     pl.median(pl.absolute(model.input_data['rel_err'].dropna())),
                                                                       model.input_data['covered?'].mean())

    model.mu = pandas.DataFrame(dict(true=pi_age_true,
                                     mu_pred=model.vars['p']['mu_age'].stats()['mean'],
                                     sigma_pred=model.vars['p']['mu_age'].stats()['standard deviation']))
    data_simulation.add_quality_metrics(model.mu)

    model.results = dict(param=[], bias=[], mare=[], mae=[], pc=[])
    data_simulation.add_to_results(model, 'delta')
    data_simulation.add_to_results(model, 'mu')
    data_simulation.add_to_results(model, 'input_data')
    model.results = pandas.DataFrame(model.results, columns='param bias mae mare pc'.split())

    print model.results

    return model


if __name__ == '__main__':
    model = validate_age_integrating_model_sim()
