""" Data models"""

import pylab as pl
import pymc as mc
import scipy.integrate

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
                                        mu_age=None, mu_age_parent=priors.get(t))

        # set initial values from data
        if t in priors:
            initial = pl.array(priors[t])
        else:
            initial = rate[t]['mu_age'].value.copy()
            mean_data = model.get_data(t).groupby(['age_start', 'age_end']).mean().delevel()
            for i, row in mean_data.T.iteritems():
                start = row['age_start'] - rate[t]['ages'][0]
                end = row['age_end'] - rate[t]['ages'][0]
                initial[start:end] = row['value']

        rate[t]['gamma_bar'].value = pl.log(initial.mean()+1.e-9)
        rate[t]['gamma'].value = pl.log(initial[rate[t]['knots'] - rate[t]['ages'][0]]+1.e-9) - pl.log(initial.mean()+1.e-9)

    m = .01*pl.ones_like(ages)
    mean_mortality = model.get_data('m').groupby(['age_start', 'age_end']).mean().delevel()
    for i, row in mean_mortality.T.iteritems():
        m[row['age_start']:row['age_end']] = row['value']

    logit_C0 = mc.Uninformative('logit_C0', value=-10.)

    # ODE functions for gradient and Jacobian
    def func(a, SC, i, r, f, m):
        a = int(a-ages[0])
        if a >= len(ages):
            return pl.array([0., 0.])
        return pl.array([-(i[a]+m[a])*SC[0] + r[a]*SC[1],
                          i[a]*SC[0] - (r[a]+m[a]+f[a])*SC[1]])
    def Dfun(a, SC, i, r, f, m):
        a = int(a-ages[0])
        if a >= len(ages):
            return pl.array([[0., 0.], [0., 0.]])
        return pl.array([[-(i[a]+m[a]), r[a]],
                         [i[a],  -(r[a]+m[a]+f[a])]])

    #rk = scipy.integrate.ode(func, Dfun).set_integrator('vode', method='adams', with_jacobian=True, rtol=.05)  # non-stiff
    rk = scipy.integrate.ode(func, Dfun).set_integrator('vode', method='bdf', with_jacobian=True)  # stiff
    #rk = scipy.integrate.ode(func, Dfun).set_integrator('vode', method='bdf', with_jacobian=True, rtol=1)  # stiff, but faster
    #rk = scipy.integrate.ode(func, Dfun).set_integrator('vode', method='adams', with_jacobian=True, nsteps=3, order=1, atol=.1)  # very non-stiff, much faster, and good results

    #rk = scipy.integrate.ode(func, Dfun).set_integrator('vode', method='adams', with_jacobian=True, nsteps=1, order=1, atol=.1)  # too non-stiff, gives errors
    #rk = scipy.integrate.ode(func).set_integrator('dopri5')  # doesn't work; why?

    @mc.deterministic
    def mu_age_p(logit_C0=logit_C0,
                 i=rate['i']['mu_age'], r=rate['r']['mu_age'], f=rate['f']['mu_age'], m=m,
                 ages=ages, rk=rk):
        C0 = mc.invlogit(logit_C0)
        SC = pl.zeros((len(ages),2))
        SC[0,:] = pl.array([1.-C0, C0])

        rk.set_f_params(i, r, f, m)
        rk.set_jac_params(i, r, f, m)
        rk.set_initial_value(SC[0,:], ages[0])
        while rk.t < ages[-1]:
            rk.integrate(rk.t+1.)
            SC[rk.t-ages[0],:] = rk.y

        return SC[:,1] / SC.sum(1)

    p = data_model.data_model('p', model, 'p',
                              root_area, root_sex, root_year,
                              mu_age_p, mu_age_parent=priors.get('p'))

    @mc.deterministic
    def mu_age_pf(p=p['mu_age'], f=rate['f']['mu_age']):
        return p*f
    pf = data_model.data_model('pf', model, 'pf',
                               root_area, root_sex, root_year,
                               mu_age_pf, mu_age_parent=priors.get('pf'))  # TODO: decide if including pf in priors is a good ideas, allow lower-bound data with data_type == csmr

    vars = rate
    vars.update(logit_C0=logit_C0, mu_age_p=mu_age_p, p=p, pf=pf)
    return vars


"""
    # iterative solution to difference equations to obtain bin sizes for all ages
    import scipy.linalg
    @mc.deterministic
    def mu_age_p(logit_C0=logit_C0, i=i['mu_age'], r=r['mu_age'], f=f['mu_age'], m_all_cause=m_all_cause, age_mesh=p_age_mesh):
        SC = pl.zeros([2, len(age_mesh)])
        p = pl.zeros(len(age_mesh))
        m = pl.zeros(len(age_mesh))
        
        SC[0,0] = 1 - mc.invlogit(logit_C0)
        SC[1,0] = 1. - SC[0,0]

        p[0] = SC[0,1] / (SC[0,0] + SC[0,1])
        m[0] = (m_all_cause[age_mesh[0]] - f[age_mesh[0]] * p[0]).clip(
            .1*m_all_cause[age_mesh[0]],
             1-1e-7)  # clip m[0] to avoid numerical instability

        for ii, a in enumerate(age_mesh[:-1]):
            A = pl.array([[-i[a]-m[ii], r[a]           ],
                          [ i[a]     , -r[a]-m[ii]-f[a]]]) * (age_mesh[ii+1] - age_mesh[ii])

            SC[:,ii+1] = pl.dot(scipy.linalg.expm(A), SC[:,ii])
            
            p[ii+1] = (SC[1,ii+1] / (SC[0,ii+1] + SC[1,ii+1])).clip(
                1e-7,
                1-1e-7)
            m[ii+1] = (m_all_cause[age_mesh[ii+1]] - f[age_mesh[ii+1]] * p[ii+1]).clip(
                .1*m_all_cause[age_mesh[ii+1]],
                 pl.inf)
        return scipy.interpolate.interp1d(age_mesh, p, kind='linear')(ages)
"""
