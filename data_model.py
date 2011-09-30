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
               mu_age, mu_age_parent,
               rate_type='neg_binom'):
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
    ages = pl.array(model.parameters['ages'])
    data = model.get_data(data_type)
    parameters = model.parameters.get(data_type, {})
    area_hierarchy = model.hierarchy

    vars = dict(data=data)

    if 'parameter_age_mesh' in parameters:
        knots = pl.array(parameters['parameter_age_mesh'])
    else:
        knots = pl.arange(ages[0], ages[-1]+1, 5)

    sigma_dict = {'No Prior':0., 'Slightly':.5, 'Moderately': .1, 'Very': .01}
    if 'smoothness' in parameters:
        sigma = sigma_dict[parameters['smoothness']['amount']]
    else:
        sigma = 0.

    if mu_age == None:
        vars.update(
            age_pattern.pcgp(name, ages=ages, knots=knots, sigma=sigma)
            )
    else:
        vars.update(dict(mu_age=mu_age))

    vars.update(expert_prior_model.level_constraints(name, parameters, vars['mu_age'], ages))
    vars.update(expert_prior_model.derivative_constraints(name, parameters, vars['mu_age'], ages))

    if mu_age_parent != None:
        # setup a hierarchical prior on the simliarity between the
        # consistent estimate here and (inconsistent) estimate for its
        # parent in the areas hierarchy
        if len(area_hierarchy.predecessors(root_area)) == 1:
            parent = area_hierarchy.predecessors(root_area)[0]
            vars.update(
                similarity_prior_model.similar('parent_similarity_%s'%name, vars['mu_age'], mu_age_parent, area_hierarchy[parent][root_area]['weight']) 
                )

        # also use this as the initial value for the age pattern, if it is not already specified
        if mu_age == None:
            # TODO: make this work when mu_age_parent is a stoch or an array (and test it)
            vars['gamma_bar'].value = pl.log(mu_age_parent.mean()).clip(-12,6)
            vars['gamma'].value = (pl.log(mu_age_parent[knots-ages[0]]) - vars['gamma_bar'].value).clip(-12,6)

    if len(data) > 0:
        age_weights = pl.ones_like(vars['mu_age'].value) # TODO: use age pattern appropriate to the rate type
        vars.update(
            age_integrating_model.age_standardize_approx(name, age_weights, vars['mu_age'], data['age_start'], data['age_end'], ages)
            #age_integrating_model.midpoint_approx(name, vars['mu_age'], data['age_start'], data['age_end'], ages)
            )

        vars.update(
            covariate_model.mean_covariate_model(name, vars['mu_interval'], data, model.output_template, area_hierarchy, root_area, root_sex, root_year)
            )

        ## ensure that all data has uncertainty quantified appropriately
        # first replace all missing se from ci
        missing_se = pl.isnan(data['standard_error']) | (data['standard_error'] <= 0)
        data['standard_error'][missing_se] = (data['upper_ci'][missing_se] - data['lower_ci'][missing_se]) / (2*1.96)

        # then replace all missing ess with se
        missing_ess = pl.isnan(data['effective_sample_size'])
        data['effective_sample_size'][missing_ess] = data['value'][missing_ess]*(1-data['value'][missing_ess])/data['standard_error'][missing_ess]**2

        if rate_type == 'neg_binom':
            vars.update(
                covariate_model.dispersion_covariate_model(name, data)
                )

            vars.update(
                rate_model.neg_binom_model(name, vars['pi'], vars['delta'], data['value'], data['effective_sample_size'])
                )
        elif rate_type == 'log_normal':
            vars['sigma'] = mc.Gamma('sigma_%s'%name, alpha=.1, beta=.1, value=.01)
            vars.update(
                rate_model.log_normal_model(name, vars['pi'], vars['sigma'], data['value'], data['standard_error'])
                )
        else:
            raise Exception, 'rate_model "%s" not implemented' % rate_model
    return vars

