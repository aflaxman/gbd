""" Data models"""

import pylab as pl
import pymc as mc
import networkx as nx
import pandas

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

def data_model(name, model, data_type, root_area, root_sex, root_year,
               mu_age, mu_age_parent, sigma_age_parent,
               rate_type='neg_binom',
               lower_bound=None):
    """ Generate PyMC objects for model of epidemological age-interval data

    Parameters
    ----------
    name : str
    model : data.ModelData
    data_type : str, i r f p pf
    root_area, root_sex, root_year : the node of the model to fit consistently
    mu_age : pymc.Node
      will be used as the age pattern, set to None if not needed
    mu_age_parent : pymc.Node
      will be used as the age pattern of the parent of the root area, set to None if not needed
    rate_type : str, optional
    
    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for each row of data
    """

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

    vars = dict(data=data)

    if 'parameter_age_mesh' in parameters:
        knots = pl.array(parameters['parameter_age_mesh'])
    else:
        knots = pl.arange(ages[0], ages[-1]+1, 5)

    sigma_dict = {'No Prior':pl.inf, 'Slightly':.1, 'Moderately': .01, 'Very': .001}
    if 'smoothness' in parameters:
        sigma = sigma_dict[parameters['smoothness']['amount']]
    else:
        sigma = 0.

    if mu_age == None:
        vars.update(
            age_pattern.pcgp(name, ages=ages, knots=knots, sigma=sigma)
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
            similarity_prior_model.similar('parent_similarity_%s'%name, vars['mu_age'], mu_age_parent, sigma_age_parent, .1)
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


        vars.update(
            covariate_model.mean_covariate_model(name, vars['mu_interval'], data, parameters, model, root_area, root_sex, root_year)
            )

        ## ensure that all data has uncertainty quantified appropriately
        # first replace all missing se from ci
        missing_se = pl.isnan(data['standard_error']) | (data['standard_error'] <= 0)
        data['standard_error'][missing_se] = (data['upper_ci'][missing_se] - data['lower_ci'][missing_se]) / (2*1.96)

        # then replace all missing ess with se
        missing_ess = pl.isnan(data['effective_sample_size'])
        data['effective_sample_size'][missing_ess] = data['value'][missing_ess]*(1-data['value'][missing_ess])/data['standard_error'][missing_ess]**2

        if rate_type == 'neg_binom':

            # warn and drop data that doesn't have effective sample size quantified, or is is non-positive
            missing_ess = pl.isnan(data['effective_sample_size']) | (data['effective_sample_size'] <= 0)
            if sum(missing_ess) > 0:
                print 'WARNING: %d rows of %s data has no quantification of uncertainty.' % (sum(missing_ess), name)
                data['effective_sample_size'][missing_ess] = 1.0

            # warn and change data where ess is unreasonably huge
            large_ess = data['effective_sample_size'] >= 1.e10
            if sum(large_ess) > 0:
                print 'WARNING: %d rows of %s data have effective sample size exceeding 10 billion.' % (sum(large_ess), name)
                data['effective_sample_size'][large_ess] = 1.e10


            if 'heterogeneity' in parameters:
                lower_dict = {'Slightly': .1, 'Moderately': .01, 'Very': .001}
                lower = lower_dict[parameters['heterogeneity']]
            else:
                lower = .1

                # special case, treat pf data more like poisson
                if data_type == 'pf':
                    lower = 10.

            # uncomment the following to make negative binomial effectively a poisson
            #lower=1.e6
            
            vars.update(
                covariate_model.dispersion_covariate_model(name, data, lower, lower*10)
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

            vars['sigma'] = mc.Uniform('sigma_%s'%name, lower=.0001, upper=.1, value=.01)
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
        else:
            raise Exception, 'rate_model "%s" not implemented' % rate_type
    else:
        vars.update(
            covariate_model.mean_covariate_model(name, [], data, parameters, model, root_area, root_sex, root_year)
            )
    vars.update(expert_prior_model.covariate_level_constraints(name, model, vars, ages))


    if lower_bound and len(lb_data) > 0:
        vars['lb'] = age_integrating_model.age_standardize_approx('lb_%s'%name, age_weights, vars['mu_age'], lb_data['age_start'], lb_data['age_end'], ages)

        vars['lb'].update(
            covariate_model.mean_covariate_model('lb_%s'%name, vars['lb']['mu_interval'], lb_data, parameters, model, root_area, root_sex, root_year)
            )

        vars['lb'].update(
            covariate_model.dispersion_covariate_model('lb_%s'%name, lb_data, 10, 100)  # treat like poisson
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

    return vars

