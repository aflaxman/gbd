""" Dismod-MR model creation methods"""

import pylab as pl
import pymc as mc
import networkx as nx
import pandas

import dismod3

import data
import rate_model
import age_pattern
import age_integrating_model
import covariate_model
import similarity_prior_model
import expert_prior_model
reload(expert_prior_model)
reload(similarity_prior_model)
reload(age_pattern)
reload(covariate_model)
reload(rate_model)

def age_specific_rate(model, data_type, reference_area='all', reference_sex='total', reference_year='all',
                      mu_age=None, mu_age_parent=None, sigma_age_parent=None, 
                      rate_type='neg_binom', lower_bound=None, interpolation_method='linear',
                      include_covariates=True, zero_re=False):
    # TODO: docstring,
    # expose (and document) interface for alternative rate_type as well as other options,
    # record reference values in the model
    """ Generate PyMC objects for model of epidemological age-interval data

    :Parameters:
      - `name` : str
      - `model` : data.ModelData
      - `data_type` : str, one of 'i', 'r', 'f', 'p', or 'pf'
      - `reference_area, reference_sex, reference_year` : the node of the model to fit consistently
      - `mu_age` : pymc.Node, will be used as the age pattern, set to None if not needed
      - `mu_age_parent` : pymc.Node, will be used as the age pattern of the parent of the root area, set to None if not needed
      - `sigma_age_parent` : pymc.Node, will be used as the standard deviation of the age pattern, set to None if not needed
      - `rate_type` : str, optional. One of 'beta_binom', 'binom', 'log_normal_model', 'neg_binom', 'neg_binom_lower_bound_model', 'neg_binom_model', 'normal_model', 'offest_log_normal', or 'poisson'
      - `lower_bound` : 
      - `interpolation_method` : str, optional, one of 'linear', 'nearest', 'zero', 'slinear', 'quadratic, or 'cubic'
      - `include_covariates` : boolean
      - `zero_re` : boolean, change one stoch from each set of siblings in area hierarchy to a 'sum to zero' deterministic
      - `TODO` : add to docstring about other options, and values allowed for them

    :Results:
      - Returns dict of PyMC objects, including 'pi', the covariate adjusted predicted values for each row of data

    """
    name = data_type
    import data
    result = data.ModelVars()
    
    if (mu_age_parent != None and pl.any(pl.isnan(mu_age_parent))) \
           or (sigma_age_parent != None and pl.any(pl.isnan(sigma_age_parent))):
        mu_age_parent = None
        sigma_age_parent = None
        print 'WARNING: nan found in parent mu/sigma.  Ignoring'

    ages = pl.array(model.parameters['ages'])
    data = model.get_data(data_type)
    if lower_bound:
        lb_data = model.get_data(lower_bound)
    parameters = model.parameters.get(data_type, {})
    area_hierarchy = model.hierarchy

    vars = dismod3.data.ModelVars()
    vars += dict(data=data)

    if 'parameter_age_mesh' in parameters:
        knots = pl.array(parameters['parameter_age_mesh'])
    else:
        knots = pl.arange(ages[0], ages[-1]+1, 5)

    smoothing_dict = {'No Prior':pl.inf, 'Slightly':1., 'Moderately': .5, 'Very': .05}
    if 'smoothness' in parameters:
        smoothing = smoothing_dict[parameters['smoothness']['amount']]
    else:
        smoothing = 0.

    if mu_age == None:
        vars.update(
            age_pattern.age_pattern(name, ages=ages, knots=knots, smoothing=smoothing, interpolation_method=interpolation_method)
            )
    else:
        vars.update(dict(mu_age=mu_age, ages=ages))

    vars.update(expert_prior_model.level_constraints(name, parameters, vars['mu_age'], ages))
    vars.update(expert_prior_model.derivative_constraints(name, parameters, vars['mu_age'], ages))

    if mu_age_parent != None:
        # setup a hierarchical prior on the simliarity between the
        # consistent estimate here and (inconsistent) estimate for its
        # parent in the areas hierarchy
        #weight_dict = {'Unusable': 10., 'Slightly': 10., 'Moderately': 1., 'Very': .1}
        #weight = weight_dict[parameters['heterogeneity']]
        vars.update(
            similarity_prior_model.similar('parent_similarity_%s'%name, vars['mu_age'], mu_age_parent, sigma_age_parent, 0.)
            )

        # also use this as the initial value for the age pattern, if it is not already specified
        if mu_age == None:
            if isinstance(mu_age_parent, mc.Node):  # TODO: test this code
                initial_mu = mu_age_parent.value
            else:
                initial_mu = mu_age_parent
                
            for i, k_i in enumerate(knots):
                vars['gamma'][i].value = (pl.log(initial_mu[k_i-ages[0]])).clip(-12,6)

    age_weights = pl.ones_like(vars['mu_age'].value) # TODO: use age pattern appropriate to the rate type
    if len(data) > 0:
        vars.update(
            age_integrating_model.age_standardize_approx(name, age_weights, vars['mu_age'], data['age_start'], data['age_end'], ages)
            )

        # uncomment the following to effectively remove alleffects
        #if 'random_effects' in parameters:
        #    for i in range(5):
        #        effect = 'sigma_alpha_%s_%d' % (name, i)
        #        parameters['random_effects'][effect] = dict(dist='TruncatedNormal', mu=.0001, sigma=.00001, lower=.00009, upper=.00011)
        #if 'fixed_effects' in parameters:
        #    for effect in ['x_sex', 'x_LDI_id_Updated_7July2011']:
        #        parameters['fixed_effects'][effect] = dict(dist='normal', mu=.0001, sigma=.00001)

        if include_covariates:
            vars.update(
                covariate_model.mean_covariate_model(name, vars['mu_interval'], data, parameters, model, reference_area, reference_sex, reference_year, zero_re=zero_re)
                )
        else:
            vars.update({'pi': vars['mu_interval']})

        ## ensure that all data has uncertainty quantified appropriately
        # first replace all missing se from ci
        missing_se = pl.isnan(data['standard_error']) | (data['standard_error'] < 0)
        data['standard_error'][missing_se] = (data['upper_ci'][missing_se] - data['lower_ci'][missing_se]) / (2*1.96)

        # then replace all missing ess with se
        missing_ess = pl.isnan(data['effective_sample_size'])
        data['effective_sample_size'][missing_ess] = data['value'][missing_ess]*(1-data['value'][missing_ess])/data['standard_error'][missing_ess]**2

        if rate_type == 'neg_binom':

            # warn and drop data that doesn't have effective sample size quantified, or is is non-positive
            missing_ess = pl.isnan(data['effective_sample_size']) | (data['effective_sample_size'] < 0)
            if sum(missing_ess) > 0:
                print 'WARNING: %d rows of %s data has invalid quantification of uncertainty.' % (sum(missing_ess), name)
                data['effective_sample_size'][missing_ess] = 0.0

            # warn and change data where ess is unreasonably huge
            large_ess = data['effective_sample_size'] >= 1.e10
            if sum(large_ess) > 0:
                print 'WARNING: %d rows of %s data have effective sample size exceeding 10 billion.' % (sum(large_ess), name)
                data['effective_sample_size'][large_ess] = 1.e10


            if 'heterogeneity' in parameters:
                lower_dict = {'Slightly': 9., 'Moderately': 3., 'Very': 1.}
                lower = lower_dict[parameters['heterogeneity']]
            else:
                lower = 1.

            # special case, treat pf data as poisson
            if data_type == 'pf':
                lower = 1.e12
            
            vars.update(
                covariate_model.dispersion_covariate_model(name, data, lower, lower*9.)
                )

            vars.update(
                rate_model.neg_binom_model(name, vars['pi'], vars['delta'], data['value'], data['effective_sample_size'])
                )
        elif rate_type == 'log_normal':

            # warn and drop data that doesn't have effective sample size quantified
            missing = pl.isnan(data['standard_error']) | (data['standard_error'] < 0)
            if sum(missing) > 0:
                print 'WARNING: %d rows of %s data has no quantification of uncertainty.' % (sum(missing), name)
                data['standard_error'][missing] = 1.e6

            # TODO: allow options for alternative priors for sigma
            vars['sigma'] = mc.Uniform('sigma_%s'%name, lower=.0001, upper=1., value=.01)
            #vars['sigma'] = mc.Exponential('sigma_%s'%name, beta=100., value=.01)
            vars.update(
                rate_model.log_normal_model(name, vars['pi'], vars['sigma'], data['value'], data['standard_error'])
                )
        elif rate_type == 'normal':

            # warn and drop data that doesn't have standard error quantified
            missing = pl.isnan(data['standard_error']) | (data['standard_error'] < 0)
            if sum(missing) > 0:
                print 'WARNING: %d rows of %s data has no quantification of uncertainty.' % (sum(missing), name)
                data['standard_error'][missing] = 1.e6

            vars['sigma'] = mc.Uniform('sigma_%s'%name, lower=.0001, upper=.1, value=.01)
            vars.update(
                rate_model.normal_model(name, vars['pi'], vars['sigma'], data['value'], data['standard_error'])
                )
        elif rate_type == 'binom':
            vars += rate_model.binom(name, vars['pi'], data['value'], data['effective_sample_size'])
        elif rate_type == 'beta_binom':
            vars += rate_model.beta_binom(name, vars['pi'], data['value'], data['effective_sample_size'])
        elif rate_type == 'poisson':
            missing_ess = pl.isnan(data['effective_sample_size']) | (data['effective_sample_size'] < 0)
            if sum(missing_ess) > 0:
                print 'WARNING: %d rows of %s data has invalid quantification of uncertainty.' % (sum(missing_ess), name)
                data['effective_sample_size'][missing_ess] = 0.0

            vars += rate_model.poisson(name, vars['pi'], data['value'], data['effective_sample_size'])
        elif rate_type == 'offset_log_normal':
            vars['sigma'] = mc.Uniform('sigma_%s'%name, lower=.0001, upper=10., value=.01)
            vars += rate_model.offset_log_normal(name, vars['pi'], vars['sigma'], data['value'], data['standard_error'])
        else:
            raise Exception, 'rate_model "%s" not implemented' % rate_type
    else:
        if include_covariates:
            vars.update(
                covariate_model.mean_covariate_model(name, [], data, parameters, model, reference_area, reference_sex, reference_year, zero_re=zero_re)
                )
    if include_covariates:
        vars.update(expert_prior_model.covariate_level_constraints(name, model, vars, ages))


    if lower_bound and len(lb_data) > 0:
        vars['lb'] = age_integrating_model.age_standardize_approx('lb_%s'%name, age_weights, vars['mu_age'], lb_data['age_start'], lb_data['age_end'], ages)

        if include_covariates:

            vars['lb'].update(
                covariate_model.mean_covariate_model('lb_%s'%name, vars['lb']['mu_interval'], lb_data, parameters, model, reference_area, reference_sex, reference_year, zero_re=zero_re)
                )
        else:
            vars['lb'].update({'pi': vars['lb']['mu_interval']})

        vars['lb'].update(
            covariate_model.dispersion_covariate_model('lb_%s'%name, lb_data, 1e12, 1e13)  # treat like poisson
            )

        ## ensure that all data has uncertainty quantified appropriately
        # first replace all missing se from ci
        missing_se = pl.isnan(lb_data['standard_error']) | (lb_data['standard_error'] <= 0)
        lb_data['standard_error'][missing_se] = (lb_data['upper_ci'][missing_se] - lb_data['lower_ci'][missing_se]) / (2*1.96)

        # then replace all missing ess with se
        missing_ess = pl.isnan(lb_data['effective_sample_size'])
        lb_data['effective_sample_size'][missing_ess] = lb_data['value'][missing_ess]*(1-lb_data['value'][missing_ess])/lb_data['standard_error'][missing_ess]**2

        # warn and drop lb_data that doesn't have effective sample size quantified
        missing_ess = pl.isnan(lb_data['effective_sample_size']) | (lb_data['effective_sample_size'] <= 0)
        if sum(missing_ess) > 0:
            print 'WARNING: %d rows of %s lower bound data has no quantification of uncertainty.' % (sum(missing_ess), name)
            lb_data['effective_sample_size'][missing_ess] = 1.0

        vars['lb'].update(
            rate_model.neg_binom_lower_bound_model('lb_%s'%name, vars['lb']['pi'], vars['lb']['delta'], lb_data['value'], lb_data['effective_sample_size'])
            )

    result[data_type] = vars
    return result
    
    
    
# TODO: refactor consistent_model.consistent_model into ism.consistent
import consistent_model
reload(consistent_model)
def consistent(model, reference_area='all', reference_sex='total', reference_year='all', priors={}):
    """ Create a consistent model
    
    :Parameters:
      - `model` : data.ModelData
      - `data_type` : str, one of 'i', 'r', 'f', 'p', or 'pf'
      - `root_area, root_sex, root_year` : the node of the model to fit consistently
      - `priors` : dictionary
    
    .. note::
      - dict priors can contain keys (t, 'mu') and (t, 'sigma') to tell the consistent model about the priors on levels for the age-specific rate of type t (these are arrays for mean and standard deviation a priori for mu_age[t]
      - it can also contain dicts keyed by t alone to insert empirical priors on the fixed effects and random effects

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
    """ Extract effect coefficients from model vars for rate type
    
    :Parameters:
      - `model` : data.ModelData
      - `type` : str, one of 'i', 'r', 'f', 'p', or 'pf'
    
    """
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
