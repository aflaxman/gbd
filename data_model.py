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
reload(similarity_prior_model)
reload(age_pattern)
reload(covariate_model)

def data_model(name, data, hierarchy, node, mu_age=None):
    """ Generate PyMC objects for model of epidemological age-interval data

    Parameters
    ----------
    name : str
    data : pandas.DataFrame
      data.columns must include value, sex, area, age_start, age_end, year_start,
      year_end, effective_sample_size, and each row will be included in the likelihood
    
    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for each row of data
    """

    if mu_age == None:
        vars = {}

        vars.update(
            age_pattern.pcgp(name, ages=pl.arange(101), knots=pl.arange(0,101,5), rho=40.)
            )
    else:
        vars = dict(mu_age=mu_age)

    age_weights = pl.ones_like(vars['mu_age'].value) # TODO: use age pattern appropriate to the rate type
    vars.update(
        age_integrating_model.age_standardize_approx(name, age_weights, vars['mu_age'], data['age_start'], data['age_end'])
        #age_integrating_model.midpoint_approx(name, vars['mu_age'], data['age_start'], data['age_end'])
        )

    vars.update(
        covariate_model.mean_covariate_model(name, vars['mu_interval'], data, hierarchy, node)
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
    if not vars['alpha']:
        alpha_trace = pl.array([])
    else:
        alpha_trace = vars['alpha'].trace()

    if not vars['beta']:
        beta_trace = pl.array([])
    else:
        beta_trace = vars['beta'].trace()

    if len(alpha_trace) == 0 and len(beta_trace) == 0:
        return vars['mu_age'].trace()

    leaves = [n for n in nx.traversal.bfs_tree(hierarchy, area) if hierarchy.successors(n) == []] or [area]

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
                U_l.ix[0, node] = 1.

            log_shift_l += pl.dot(alpha_trace, U_l)
            
        # make X_l
        if len(beta_trace) > 0:
            X_l = covs.ix[l, sex, year]
            log_shift_l += pl.dot(beta_trace, X_l)

        shift_l = pl.exp(log_shift_l)
        covariate_shift += shift_l * output_template['pop'][l,sex,year]
        total_population += output_template['pop'][l,sex,year]
    covariate_shift /= total_population

    return (vars['mu_age'].trace().T*covariate_shift).T
    
