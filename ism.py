""" Stub module for planned refactoring of dismod3 model creation methods"""


# TODO: refactor data_model.data_model into ism.age_specific_rate
import data_model
import data
def age_specific_rate(model, data_type, reference_area='all', reference_sex='total', reference_year='all'):
    result = data.ModelVars()
    result[data_type] = data_model.data_model(data_type, model, data_type,
                                              reference_area, reference_sex, reference_year,
                                              None, None, None)
    return result


# TODO: refactor consistent_model.consistent_model into ism.consistent
import consistent_model
reload(consistent_model)
def consistent(model, reference_area='all', reference_sex='total', reference_year='all'):
    priors = {}
    return consistent_model.consistent_model(model, reference_area, reference_sex, reference_year, priors)


# TODO: refactor emp_priors into a class and document them
def emp_priors(dm, reference_area, reference_sex, reference_year):
    import dismod3.utils
    param_type = dict(i='incidence', p='prevalence', r='remission', f='excess-mortality', rr='relative-risk', pf='prevalence_x_excess-mortality', m_with='mortality')
    emp_priors = {}
    for t in 'i r pf p rr f'.split():
        key = dismod3.utils.gbd_key_for(param_type[t], reference_area, reference_year, reference_sex)
        mu = dm.get_mcmc('emp_prior_mean', key)
        sigma = dm.get_mcmc('emp_prior_std', key)
        
        if len(mu) == 101 and len(sigma) == 101:
            emp_priors[t, 'mu'] = mu
            emp_priors[t, 'sigma'] = sigma

    return emp_priors
