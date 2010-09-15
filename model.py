# -*- coding: utf-8 -*-
""" Functions to generate PyMC models of spatio-temporal health data
"""

import pylab as pl
import pymc as mc
import pymc.gp as gp

def count_covariates(data):
    return len([n for n in data.dtype.names if n.startswith('x')])

def fe(data):
    """ Fixed Effect model::
    
        Y_r,c,t = beta * X_r,c,t + e_r,c,t
        e_r,c,t ~ N(0, sigma^2)
    """
    # covariates
    K = count_covariates(data)
    X = pl.array([data['x%d'%i] for i in range(K)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    
    # predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        return pl.dot(beta, X)
    param_predicted = mu

    @mc.deterministic
    def predicted(mu=mu, sigma_e=sigma_e):
        return mc.rnormal(mu, sigma_e**-2.)

    # likelihood
    i_obs = pl.find(1 - pl.isnan(data.y))
    @mc.observed
    def obs(value=data.y, i_obs=i_obs, mu=mu, sigma_e=sigma_e):
        return mc.normal_like(value[i_obs], mu[i_obs], sigma_e**-2.)

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)

    # find good initial conditions with MAP approx
    print 'attempting to maximize likelihood'
    var_list = [mod_mc.beta, mod_mc.obs, mod_mc.sigma_e]
    mc.MAP(var_list).fit(method='fmin_powell', verbose=1)

    return mod_mc


def gp_re_a(data):
    """ Random Effect model, with gaussian process correlations in residuals that
    includes age::
    
        Y_r,c,t,a = beta * X_r,c,t + f_r(t) + g_r(a) + e_r,c,t,a

        f_r(t) ~ GP(0, C_rt)
        g_r(a) ~ GP(0, C_ra)

        C_t ~ Matern(2, sigma_f,0, tau_f,0)
        C_a ~ Matern(2, sigma_f,2, tau_f,1)
    """
    # covariates
    K = count_covariates(data)
    X = pl.array([data['x%d'%i] for i in range(K)])


    # semi-uninformative priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0., upper=1000., value=1.)

    sigma_f = mc.Gamma('sigma_f', alpha=.1, beta=.1, value=1.*pl.ones(2))
    tau_f = mc.Gamma('tau_f1', alpha=10., beta=.1, value=10.*pl.ones(2))


    # fixed-effect prediction
    @mc.deterministic
    def mu(X=X, beta=beta):
        """ mu_i,r,c,t,a = beta * X_i,r,c,t,a"""
        return pl.dot(beta, X)


    # setup gaussian processes random effects to model additional variation
    # make index dict to convert from region/country/age to array index
    regions = pl.unique(data.region)
    years = pl.unique(data.year)
    ages = pl.unique(data.age)

    r_index = dict([(r, i) for i, r in enumerate(regions)])
    t_index = dict([(t, i) for i, t in enumerate(years)])
    a_index = dict([(a, i) for i, a in enumerate(ages)])

    @mc.deterministic
    def C_t(grid=years, sigma_f=sigma_f, tau_f=tau_f, diff_degree=2.):
        return gp.matern.euclidean(grid, grid, amp=sigma_f[0], scale=tau_f[0], diff_degree=diff_degree)
    @mc.deterministic
    def C_a(grid=ages, sigma_f=sigma_f, tau_f=tau_f, diff_degree=2.):
        return gp.matern.euclidean(grid, grid, amp=sigma_f[1], scale=tau_f[1], diff_degree=diff_degree)

    f = [mc.MvNormalCov('f_%s'%r, pl.zeros_like(years), C_t, value=pl.zeros_like(years)) for r in regions]
    g = [mc.MvNormalCov('g_%s'%r, pl.zeros_like(ages), C_a, value=pl.zeros_like(ages)) for r in regions]

    # organize observations into country panels and calculate predicted value before sampling error
    country_param_pred = []
    country_tau = []
    for c in pl.unique(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]

        # find the index for this region, country, and for the relevant ages and times
        region = data.region[i_c[0]]
        r_index_c = r_index[region]
        t_index_c = [t_index[data.year[i]] for i in i_c]
        a_index_c = [a_index[data.age[i]] for i in i_c]
    
        # find predicted parameter value for all observations of country c
        @mc.deterministic(name='country_param_pred_%s'%c)
        def country_param_pred_c(i=i_c, mu=mu, f=f[r_index_c], g=g[r_index_c], a=a_index_c, t=t_index_c):
            """ country_param_pred_c[row] = parameter_predicted[row] * 1[row.country == c]"""
            country_param_pred_c = pl.zeros_like(data.y)
            country_param_pred_c[i] = mu[i] + f[t] + g[a]
            return country_param_pred_c
        country_param_pred.append(country_param_pred_c)

        # find predicted parameter precision for all observations of country c
        @mc.deterministic(name='country_tau_%s'%c)
        def country_tau_c(i=i_c, sigma_e=sigma_e, var_d_c=data.se[i_c]**2.):
            """ country_tau_c[row] = tau[row] * 1[row.country == c]"""
            country_tau_c = pl.zeros_like(data.y)
            country_tau_c[i] = 1 / (sigma_e**2. + var_d_c)
            return country_tau_c
        country_tau.append(country_tau_c)

    @mc.deterministic
    def param_predicted(country_param_pred=country_param_pred):
        return pl.sum(country_param_pred, axis=0)

    @mc.deterministic
    def tau(country_tau=country_tau):
        return pl.sum(country_tau, axis=0)

    @mc.deterministic
    def data_predicted(param_predicted=param_predicted, tau=tau):
        return mc.rnormal(param_predicted, tau)
    predicted = data_predicted

    i_obs = [i for i in range(len(data)) if not pl.isnan(data.y[i])]
    @mc.observed
    def obs(value=data.y, i=i_obs, param_predicted=param_predicted, tau=tau):
        return mc.normal_like(value[i], param_predicted[i], tau[i])

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    
    # TODO:  use covariance matrix to seed adaptive metropolis steps
    for r in range(len(regions)):
        mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.f[r], cov=pl.array(C_t.value*.01))
        mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.g[r], cov=pl.array(C_a.value*.01))

    # find good initial conditions with MAP approx
    for var_list in [[mod_mc.obs, mod_mc.beta]] + \
        [[mod_mc.obs, f_r] for f_r in mod_mc.f] + \
        [[mod_mc.obs, g_r] for g_r in mod_mc.g] + \
        [[mod_mc.obs, mod_mc.beta, mod_mc.sigma_e]] + \
        [[mod_mc.obs, mod_mc.beta] + mod_mc.f] + \
        [[mod_mc.obs, mod_mc.beta] + mod_mc.g] + \
        [[mod_mc.obs, mod_mc.beta] + mod_mc.f +  mod_mc.g]:
        print 'attempting to maximize likelihood of %s' % [v.__name__ for v in var_list]
        mc.MAP(var_list).fit(method='fmin_powell', verbose=1)

    return mod_mc


