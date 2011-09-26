""" Data models"""

import pylab as pl
import pymc as mc

import data
import rate_model
import age_pattern
import age_integrating_model
import covariate_model
reload(age_pattern)
reload(covariate_model)

def data_model(name, data, hierarchy):
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

    vars = {}

    vars.update(
        age_pattern.pcgp(name, ages=pl.arange(101), knots=pl.arange(0,101,5), rho=40.)
        )

    age_weights = pl.ones_like(vars['mu_age'].value) # TODO: use age pattern appropriate to the rate type
    vars.update(
        age_integrating_model.age_standardize_approx(name, age_weights, vars['mu_age'], data['age_start'], data['age_end'])
        #age_integrating_model.midpoint_approx(name, vars['mu_age'], data['age_start'], data['age_end'])
        )

    vars.update(
        covariate_model.mean_covariate_model(name, vars['mu_interval'], data, hierarchy, 'all')
        )

    vars.update(
        covariate_model.dispersion_covariate_model(name, data)
        )

    vars.update(
        rate_model.neg_binom_model(name, vars['pi'], vars['delta'], data['value'], data['effective_sample_size'])
        )

    return vars
