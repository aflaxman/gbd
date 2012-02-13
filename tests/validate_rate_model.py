""" Validate Rate Model

Out-of-sample predictive validity approach

1. Load epilepsy model 32377
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

def validate_rate_model(rate_type='neg_binom', replicate=0):
    # set random seed for reproducibility
    mc.np.random.seed(123567 + replicate)
    
    # load epilepsy data
    model = dismod3.data.load('/home/j/Project/dismod/output/dm-32377/')
    data = model.get_data('p')

    # sample prevalence data
    i_test = mc.rbernoulli(.25, size=len(data.index))

    data['standard_error'][i_test] = pl.nan
    data['lower_ci'][i_test] = pl.nan
    data['upper_ci'][i_test] = pl.nan
    data['effective_sample_size'][i_test] = 0.
    
    model.input_data = data


    # create model
    # TODO: set parameters in model.parameters['p'] dict
    # then have simple method to create age specific rate model
    #model.parameters['p'] = ...
    #model.vars += dismod3.ism.age_specific_rate(model, 'p')

    model.parameters['p']['parameter_age_mesh'] = [0,100]
    model.vars['p'] = dismod3.data_model.data_model(
        'p', model, 'p',
        'all', 'total', 'all',
        None, None, None,
        rate_type=rate_type,
        interpolation_method='zero',
        include_covariates=False)
    
    # fit model
    dismod3.fit.fit_asr(model, 'p', iter=100, thin=1, burn=0)

    # compare estimate to hold-out
    data['mu_pred'] = model.vars['p']['p_pred'].stats()['mean']
    data['sigma_pred'] = model.vars['p']['p_pred'].stats()['standard deviation']

    import data_simulation
    data['true'] = data['value']
    data_simulation.add_quality_metrics(data)

    data_simulation.initialize_results(model)
    data_simulation.add_to_results(model, 'input_data')
    data_simulation.finalize_results(model)


    return model
