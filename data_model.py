""" Data models"""

import pylab as pl
import pymc as mc

import data
import rate_model
import age_pattern
import age_integrating_model
import covariate_model
reload(covariate_model)

def data_model(name, data, hierarchy):
    """ Generate PyMC objects for model of epidemological age-interval data

    Parameters
    ----------
    name : str
    data : pandas.DataFrame
      data.columns must include value, age_start, age_end, year_start,
      year_end, effective_sample_size, and each row will be included in the likelihood
    
    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for each row of data
    """

    # make the "design matrices" X and U
    X = data.select(lambda col: col.startswith('x_'), axis=1)
    U = data.select(lambda col: col.startswith('u_'), axis=1)

    vars = {}

    vars.update(
        gamma_bar=mc.Uninformative('gamma_bar_%s'%name, value=0.),
        beta=mc.Uninformative('beta_%s'%name, value=pl.zeros_like(X.columns)),
        eta=mc.Normal('eta_%s'%name, mu=5., tau=1., value=5.)
        )

    if U:
        vars.update(
            zeta=mc.Uninformative('zeta_%s'%name, value=pl.zeros_like(U.columns))
            )
        vars.update(
            covariate_model.dispersion_covariate_model(name, vars['eta'], vars['zeta'], U),
            )
    else:
        vars.update(
            delta=mc.Lambda('delta_%s'%name, lambda eta=vars['eta']: 50. + pl.exp(eta))
            )


    vars.update(
        age_pattern.pcgp(name, vars['gamma_bar'], knots=pl.arange(0,101,5), rho=40.)
        )

    vars.update(
        age_integrating_model.midpoint_approx(name, vars['mu_age'], data['age_start'], data['age_end'])
        )

    vars.update(
        covariate_model.mean_covariate_model(name, vars['mu_interval'], vars['beta'], X, hierarchy)
        )

    vars.update(
        rate_model.neg_binom_model(name, vars['pi'], vars['delta'], data['value'], data['effective_sample_size'])
        )

    return vars
