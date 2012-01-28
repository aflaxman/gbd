""" Data models"""

import pylab as pl
import pymc as mc
import scipy.interpolate

import data
import rate_model
import age_pattern
import age_integrating_model
import covariate_model
import data_model
import similarity_prior_model
reload(age_pattern)
reload(covariate_model)
reload(data_model)

def consistent_model(model, root_area, root_sex, root_year, priors):
    """ Generate PyMC objects for consistent model of epidemological data

    Parameters
    ----------
    model : data.ModelData
    root_area, root_sex, root_year : str, node in hierarchy to fit consistently
    priors : dict, with keys for data types for lists of priors on age patterns
    
    Results
    -------
    Returns dict of dicts of PyMC objects, including 'i, p, r, f', the covariate
    adjusted predicted values for each row of data
    """
    rate = {}
    ages = model.parameters['ages']

    for t in 'irf':
        rate[t] = data_model.data_model(t, model, t,
                                        root_area, root_sex, root_year,
                                        mu_age=None,
                                        mu_age_parent=priors.get((t, 'mu')),
                                        sigma_age_parent=priors.get((t, 'sigma')))

        # set initial values from data
        if t in priors:
            if isinstance(priors[t], mc.Node):
                initial = priors[t].value
            else:
                initial = pl.array(priors[t])
        else:
            initial = rate[t]['mu_age'].value.copy()
            mean_data = model.get_data(t).groupby(['age_start', 'age_end']).mean().delevel()
            for i, row in mean_data.T.iteritems():
                start = row['age_start'] - rate[t]['ages'][0]
                end = row['age_end'] - rate[t]['ages'][0]
                initial[start:end] = row['value']

        for i,k in enumerate(rate[t]['knots']):
            rate[t]['gamma'][i].value = pl.log(initial[k - rate[t]['ages'][0]]+1.e-9)

    m_all = .01*pl.ones(101)
    mean_mortality = model.get_data('m_all').groupby(['age_start', 'age_end']).mean().delevel()

    if len(mean_mortality) == 0:
        print 'WARNING: all-cause mortality data not found, using m_all = .01'
    else:
        knots = []
        for i, row in mean_mortality.T.iteritems():
            knots.append(pl.clip((row['age_start'] + row['age_end'] + 1.) / 2., 0, 100))
            
            m_all[knots[-1]] = row['value']

        # extend knots as constant beyond endpoints
        knots = sorted(knots)
        m_all[0] = m_all[knots[0]]
        m_all[100] = m_all[knots[-1]]

        knots.insert(0, 0)
        knots.append(100)

        m_all = scipy.interpolate.interp1d(knots, m_all[knots], kind='linear')(pl.arange(101))
    m_all = m_all[ages]

    logit_C0 = mc.Uninformative('logit_C0', value=-10.)


    # use Runge-Kutta 4 ODE solver
    import dismod_ode

    N = len(m_all)
    num_step = 10  # double until it works
    ages = pl.array(ages, dtype=float)
    fun = dismod_ode.ode_function(num_step, ages, m_all)

    @mc.deterministic
    def mu_age_p(logit_C0=logit_C0,
                 i=rate['i']['mu_age'],
                 r=rate['r']['mu_age'],
                 f=rate['f']['mu_age']):

        # for acute conditions, it is silly to use ODE solver to
        # derive prevalence, and it can be approximated with a simple
        # transformation of incidence
        if r.min() > 5.99:
            return i / (r + m_all + f)
        
        C0 = mc.invlogit(logit_C0)

        x = pl.hstack((i, r, f, 1-C0, C0))
        y = fun.forward(0, x)

        susceptible = y[:N]
        condition = y[N:]

        p = condition / (susceptible + condition)
        p[pl.isnan(p)] = 0.
        return p

    p = data_model.data_model('p', model, 'p',
                              root_area, root_sex, root_year,
                              mu_age_p,
                              mu_age_parent=priors.get(('p', 'mu')),
                              sigma_age_parent=priors.get(('p', 'sigma')))

    @mc.deterministic
    def mu_age_pf(p=p['mu_age'], f=rate['f']['mu_age']):
        return p*f
    pf = data_model.data_model('pf', model, 'pf',
                               root_area, root_sex, root_year,
                               mu_age_pf,
                               mu_age_parent=priors.get(('pf', 'mu')),
                               sigma_age_parent=priors.get(('pf', 'sigma')),
                               lower_bound='csmr',
                               include_covariates=False)

    @mc.deterministic
    def mu_age_m(pf=pf['mu_age'], m_all=m_all):
        return (m_all - pf).clip(1.e-6, 1.e6)
    rate['m'] = data_model.data_model('m_wo', model, 'm_wo',
                              root_area, root_sex, root_year,
                              mu_age_m,
                              None, None,
                              include_covariates=False)

    @mc.deterministic
    def mu_age_rr(m=rate['m']['mu_age'], f=rate['f']['mu_age']):
        return (m+f) / m
    rr = data_model.data_model('rr', model, 'rr',
                               root_area, root_sex, root_year,
                               mu_age_rr,
                               mu_age_parent=priors.get(('rr', 'mu')),
                               sigma_age_parent=priors.get(('rr', 'sigma')),
                               rate_type='log_normal',
                               include_covariates=False)

    @mc.deterministic
    def mu_age_smr(m=rate['m']['mu_age'], f=rate['f']['mu_age'], m_all=m_all):
        return (m+f) / m_all
    smr = data_model.data_model('smr', model, 'smr',
                                root_area, root_sex, root_year,
                                mu_age_smr,
                                mu_age_parent=priors.get(('smr', 'mu')),
                                sigma_age_parent=priors.get(('smr', 'sigma')),
                                rate_type='log_normal',
                                include_covariates=False)

    @mc.deterministic
    def mu_age_m_with(m=rate['m']['mu_age'], f=rate['f']['mu_age']):
        return m+f
    m_with = data_model.data_model('m_with', model, 'm_with',
                                   root_area, root_sex, root_year,
                                   mu_age_m_with,
                                   mu_age_parent=priors.get(('m_with', 'mu')),
                                   sigma_age_parent=priors.get(('m_with', 'sigma')),
                                   include_covariates=False)
    
    # duration = E[time in bin C]
    @mc.deterministic
    def mu_age_X(r=rate['r']['mu_age'], m=rate['m']['mu_age'], f=rate['f']['mu_age']):
        hazard = r + m + f
        pr_not_exit = pl.exp(-hazard)
        X = pl.empty(len(hazard))
        X[-1] = 1 / hazard[-1]
        for i in reversed(range(len(X)-1)):
            X[i] = pr_not_exit[i] * (X[i+1] + 1) + 1 / hazard[i] * (1 - pr_not_exit[i]) - pr_not_exit[i]
        return X
    X = data_model.data_model('X', model, 'X',
                              root_area, root_sex, root_year,
                              mu_age_X,
                              mu_age_parent=priors.get(('X', 'mu')),
                              sigma_age_parent=priors.get(('X', 'sigma')),
                              rate_type='normal',
                              include_covariates=True)



    vars = rate
    vars.update(logit_C0=logit_C0, p=p, pf=pf, rr=rr, smr=smr, m_with=m_with, X=X)
    return vars
