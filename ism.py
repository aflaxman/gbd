""" Stub module for planned refactoring of dismod3 model creation methods"""

import pymc as mc


# TODO: refactor data_model.data_model into ism.age_specific_rate
import data_model
reload(data_model)
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
def consistent(model, reference_area='all', reference_sex='total', reference_year='all', priors={}):
    """ dict priors can contain keys (t, 'mu') and (t, 'sigma') to
    tell the consistent model about the priors on levels for the
    age-specific rate of type t (these are arrays for mean and standard deviation a priori for mu_age[t]

    it can also contain dicts keyed by t alone to insert empirical priors on the fixed effects and random effects
    """
    # TODO: refactor the way priors are handled
    # current approach is much more complicated than necessary
    for t in 'i r pf p rr f'.split():
        if t in priors:
            model.parameters[t]['random_effects'].update(priors[t]['random_effects'])
            model.parameters[t]['fixed_effects'].update(priors[t]['fixed_effects'])

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

def effect_priors(model, type):
    """ Extract effect coeffs from model vars for rate type"""
    vars = model.vars[type]
    prior_vals = {}
    
    prior_vals['new_alpha'] = {}
    if 'alpha' in vars:
        for n, col in zip(vars['alpha'], vars['U'].columns):
            if isinstance(n, mc.Node):
                stats = n.stats()
                if stats:
                    #prior_vals['new_alpha'][col] = dict(dist='TruncatedNormal', mu=stats['mean'], sigma=stats['standard deviation'], lower=-5., upper=5.)
                    prior_vals['new_alpha'][col] = dict(dist='Constant', mu=stats['mean'], sigma=stats['standard deviation'])

        # uncomment below to save empirical prior on sigma_alpha, the dispersion of the random effects
        for n in vars['sigma_alpha']:
            stats = n.stats()
            prior_vals['new_alpha'][n.__name__] = dict(dist='TruncatedNormal', mu=stats['mean'], sigma=stats['standard deviation'], lower=.01, upper=.5)

    prior_vals['new_beta'] = {}
    if 'beta' in vars:
        for n, col in zip(vars['beta'], vars['X'].columns):
            stats = n.stats()
            if stats:
                #prior_vals['new_beta'][col] = dict(dist='normal', mu=stats['mean'], sigma=stats['standard deviation'], lower=-pl.inf, upper=pl.inf)
                prior_vals['new_beta'][col] = dict(dist='Constant', mu=stats['mean'])

    return prior_vals
