""" Validate Rate Model

Out-of-sample predictive validity approach

1. Load epilepsy model 32377 (or schiz or simulated data)
2. sample 75% of prevalence data
3. fit it without age pattern or covariates
4. measure predictive accuracy
5. repeat many times for all rate models

"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl
import pymc as mc
import dismod3
reload(dismod3)

pl.seterr('ignore')

def validate_rate_model(rate_type='neg_binom', data_type='epilepsy', replicate=0):
    # set random seed for reproducibility
    mc.np.random.seed(1234567 + replicate)
    
    # load data
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-32377/')

    data = model.get_data('p')

    #data = data.ix[:20, :]
    
    # replace data with synthetic data if requested
    if data_type == 'epilepsy':
        # no replacement needed
        pass

    elif data_type == 'schiz':
        import pandas as pd
        data = pd.read_csv('/homes/abie/gbd_dev/gbd/tests/schiz.csv')
    
    elif data_type == 'binom':
        N = 1.e6
        data['effective_sample_size'] = N
        mu = data['value'].mean()
        data['value'] = mc.rbinomial(N, mu, size=len(data.index)) / N

    elif data_type == 'poisson':
        N = 1.e6
        data['effective_sample_size'] = N
        mu = data['value'].mean()
        data['value'] = mc.rpoisson(N*mu, size=len(data.index)) / N

    elif data_type == 'normal':
        mu = data['value'].mean()
        sigma = .125*mu
        data['standard_error'] = sigma
        data['value'] = mc.rnormal(mu, sigma**-2, size=len(data.index))

    elif data_type == 'log_normal':
        mu = data['value'].mean()
        sigma = .25
        data['standard_error'] = sigma*mu
        data['value'] = pl.exp(mc.rnormal(pl.log(mu), sigma**-2, size=len(data.index)))

    else:
        raise TypeError, 'Unknown data type "%s"' % data_type

    # sample prevalence data
    i_test = mc.rbernoulli(.25, size=len(data.index))
    i_nan = pl.isnan(data['effective_sample_size'])
    
    data['lower_ci'] = pl.nan
    data['upper_ci'] = pl.nan
    data.ix[i_nan, 'effective_sample_size'] = 0.
    data['standard_error'] = pl.sqrt(data['value']*(1-data['value'])) / data['effective_sample_size']
    data.ix[pl.isnan(data['standard_error']), 'standard_error'] = pl.inf

    data['standard_error'][i_test] = pl.inf
    data['effective_sample_size'][i_test] = 0.

    data['value'] = pl.maximum(data['value'], 1.e-12)

    model.input_data = data


    # create model
    # TODO: set parameters in model.parameters['p'] dict
    # then have simple method to create age specific rate model
    #model.parameters['p'] = ...
    #model.vars += dismod3.ism.age_specific_rate(model, 'p')

    model.parameters['p']['parameter_age_mesh'] = [0,100]
    model.parameters['p']['heterogeneity'] = 'Very'
    model.vars['p'] = dismod3.data_model.data_model(
        'p', model, 'p',
        'all', 'total', 'all',
        None, None, None,
        rate_type=rate_type,
        interpolation_method='zero',
        include_covariates=False)
    
    # add upper bound on sigma in log normal model to help convergence
    #if rate_type == 'log_normal':
    #    model.vars['p']['sigma'].parents['upper'] = 1.5

    # add upper bound on sigma, zeta in offset log normal
    #if rate_type == 'offset_log_normal':
    #    model.vars['p']['sigma'].parents['upper'] = .1
    #    model.vars['p']['p_zeta'].value = 5.e-9
    #    model.vars['p']['p_zeta'].parents['upper'] = 1.e-8

    # fit model
    dismod3.fit.fit_asr(model, 'p', iter=20000, thin=10, burn=10000)
    #dismod3.fit.fit_asr(model, 'p', iter=100, thin=1, burn=0)

    # compare estimate to hold-out
    data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    data['lb_pred'] = model.vars['p']['p_pred'].stats()['95% HPD interval'][:,0]
    data['ub_pred'] = model.vars['p']['p_pred'].stats()['95% HPD interval'][:,1]

    import data_simulation
    model.test = data[i_test]
    data = model.test
    data['true'] = data['value']
    data_simulation.add_quality_metrics(data)

    data_simulation.initialize_results(model)
    data_simulation.add_to_results(model, 'test')
    data_simulation.finalize_results(model)


    return model

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Run a rate model validation.')
    parser.add_argument('rate_type',
                        help='rate distribution name')
    parser.add_argument('data_type',
                        help='data distribution name')
    parser.add_argument('replicate', type=int,
                        help='replicate number, for saving')
    args = parser.parse_args()

    validation_name = '%s/%s/%s-%s-%s.csv' % ('/home/j/Project/dismod', 'rate_model_validation', args.rate_type, args.data_type, args.replicate)
    print 'Running validation for:', validation_name
    model = validate_rate_model(args.rate_type, args.data_type, args.replicate)
    model.results.to_csv(validation_name)
                             
