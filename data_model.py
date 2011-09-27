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

def data_model(name, data, parameters, hierarchy, root, mu_age=None, mu_age_parent=None):
    """ Generate PyMC objects for model of epidemological age-interval data

    Parameters
    ----------
    name : str
    data : pandas.DataFrame
      data.columns must include value, sex, area, age_start, age_end, year_start,
      year_end, effective_sample_size, and each row will be included in the likelihood
    mu_age : pymc.Node, optional
      will be used as the age pattern if provided
    mu_age_parent : pymc.Node, optional
      will be used as the age pattern of the parent of the root area if provided
    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for each row of data
    """
    vars = dict(data=data)

    if 'parameter_age_mesh' in parameters:
        knots = pl.array(parameters['parameter_age_mesh'])
    else:
        knots = pl.arange(0,101,5)

    rho_dict = {'No Prior':1., 'Slightly':10., 'Moderately': 20, 'Very': 40}
    if 'smoothness' in parameters:
        rho = rho_dict[parameters['smoothness']['amount']]
    else:
        rho = 10.

    if mu_age == None:
        vars.update(
            age_pattern.pcgp(name, ages=pl.arange(101), knots=knots, rho=rho)
            )
    else:
        vars.update(dict(mu_age=mu_age))

    vars.update(expert_prior_model.level_constraints(name, parameters, vars['mu_age']))

    if mu_age_parent != None:
        # setup a hierarchical prior on the simliarity between the
        # consistent estimate here and (inconsistent) estimate for its
        # parent in the areas hierarchy
        if len(hierarchy.predecessors(root)) == 1:
            parent = hierarchy.predecessors(root)[0]
            vars.update(
                similarity_prior_model.similar('parent_similarity_%s'%name, vars['mu_age'], mu_age_parent, hierarchy[parent][root]['weight'])
                )

        # also use this as the initial value for the age pattern, if it is not already specified
        if mu_age == None:
            # TODO: make this work when mu_age_parent is a stoch or an array (and test it)
            vars['gamma_bar'].value = pl.log(mu_age_parent.mean()).clip(-20,20)
            vars['gamma'].value = (pl.log(mu_age_parent[knots]) - vars['gamma_bar'].value).clip(-20,20)

    age_weights = pl.ones_like(vars['mu_age'].value) # TODO: use age pattern appropriate to the rate type
    vars.update(
        age_integrating_model.age_standardize_approx(name, age_weights, vars['mu_age'], data['age_start'], data['age_end'])
        #age_integrating_model.midpoint_approx(name, vars['mu_age'], data['age_start'], data['age_end'])
        )

    vars.update(
        covariate_model.mean_covariate_model(name, vars['mu_interval'], data, hierarchy, root)
        )

    vars.update(
        covariate_model.dispersion_covariate_model(name, data)
        )

    missing_ess = pl.isnan(data['effective_sample_size'])
    data['effective_sample_size'][missing_ess] = data['value'][missing_ess]*(1-data['value'][missing_ess])/data['standard_error'][missing_ess]**2

    vars.update(
        rate_model.neg_binom_model(name, vars['pi'], vars['delta'], data['value'], data['effective_sample_size'])
        )

    return vars


def predict_for(output_template, hierarchy, root, area, sex, year, vars):
    """ Generate draws from posterior predicted distribution for a
    specific (area, sex, year)

    Parameters
    ----------
    output_template : pandas.DataFrame with covariate data for all leaf nodes in area hierarchy
    hierarchy : nx.DiGraph encoding hierarchical relationship of areas
    root : str, area for which this model was fit consistently
    area : str, area to predict for
    sex : str, sex to predict for
    year : str, year to predict for
    vars : dict, including entries for alpha, beta, mu_age, U, and X

    Results
    -------
    Returns array of draws from posterior predicted distribution
    """

    if isinstance(vars['alpha'], mc.Node):
        alpha_trace = vars['alpha'].trace()
    else:
        alpha_trace = pl.array([])

    if isinstance(vars['beta'], mc.Node):
        beta_trace = vars['beta'].trace()
    else:
        beta_trace = pl.array([])

    if len(alpha_trace) == 0 and len(beta_trace) == 0:
        return vars['mu_age'].trace()

    leaves = [n for n in nx.traversal.bfs_tree(hierarchy, area) if hierarchy.successors(n) == []]
    if len(leaves) == 0:
        # networkx returns an empty list when the bfs tree is a single node
        leaves = [area]

    covariate_shift = 0.
    total_population = 0.

    p_U = 2 + hierarchy.number_of_nodes()  # random effects for sex, time, area
    U_l = pandas.DataFrame(pl.zeros((1, p_U)), columns=['sex', 'time'] + hierarchy.nodes())
    U_l = U_l.filter(vars['U'].columns)
    
    output_template = output_template.groupby(['area', 'sex', 'year']).mean()
    covs = output_template.filter(vars['X'].columns)
    
    for l in leaves:
        log_shift_l = 0.

        # make U_l
        if len(alpha_trace) > 0:
            U_l.ix[0, :] = 0.
            U_l.ix[0, 'sex'] = float(sex == 'male')
            U_l.ix[0, 'time'] = year - 2000.
            for node in nx.shortest_path(hierarchy, root, l):
                if node in U_l.columns:
                    U_l.ix[0, node] = 1.
                else:
                    # TODO: include appropriate uncertainty for random effects that are not in model
                    pass

            log_shift_l += pl.dot(alpha_trace, pl.atleast_2d(U_l).T)
            
        # make X_l
        if len(beta_trace) > 0:
            X_l = covs.ix[l, sex, year]
            log_shift_l += pl.dot(beta_trace, pl.atleast_2d(X_l).T)

        shift_l = pl.exp(log_shift_l)
        covariate_shift += shift_l * output_template['pop'][l,sex,year]
        total_population += output_template['pop'][l,sex,year]
    covariate_shift /= total_population

    return vars['mu_age'].trace() * covariate_shift
    
