# -*- coding: utf-8 -*-
""" Functions to generate PyMC models of spatio-temporal health data
"""

import pylab as pl
import pymc as mc
import pymc.gp as gp

def fe(data):
    """ Fixed Effect model::
    
        Y_r,c,t = beta * X_r,c,t + e_r,c,t
        e_r,c,t ~ N(0, sigma^2)
    """
    # covariates
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma = mc.Uniform('sigma', lower=0, upper=1000, value=1)
    
    # predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        return pl.dot(beta, X)

    @mc.deterministic
    def predicted(mu=mu, sigma=sigma):
        return mc.rnormal(mu, sigma**-2.)

    # likelihood
    i_obs = pl.find(1 - pl.isnan(data.y))
    @mc.observed
    def y(value=data.y, i_obs=i_obs, mu=mu, sigma=sigma):
        return mc.normal_like(value[i_obs], mu[i_obs], sigma**-2.)

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    return mod_mc


def re(data):
    """ Random Effect model::
    
        Y_r,c,t = (beta + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r,c,t[i] ~ N(0, sigma_i^2)
        e_r,c,t ~ N(0, sigma_e^2)
    """
    # covariates
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma = mc.Uniform('sigma', lower=0, upper=1000, value=pl.ones(K))
    
    # predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        return pl.dot(beta, X)

    @mc.deterministic
    def tau(X=X, sigma=sigma, sigma_e=sigma_e):
        """ y_i - mu_i ~ beta * N(0, sigma_k**2) + N(0, sigma_i**2)"""
        return 1 / (pl.dot(sigma**2., X**2.) + sigma_e**2.)

    @mc.deterministic
    def predicted(mu=mu, tau=tau):
        return mc.rnormal(mu, tau)

    # likelihood
    i_obs = pl.find(1 - pl.isnan(data.y))
    @mc.observed
    def y(value=data.y, i_obs=i_obs, mu=mu, tau=tau):
        return mc.normal_like(value[i_obs], mu[i_obs], tau[i_obs])

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    return mod_mc


def nested_re(data):
    """ Random Effect model, with country random effects nested in
    regions::
    
        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r[k] ~ N(0, sigma_k^2)
        u_r,c,t[k] ~ N(0, sigma_r,k^2)
        e_r,c,t ~ N(0, sigma_e^2)

    """
    # covariates, etc
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    regions = set(data.region)
    R = len(regions)
    region_id = dict(zip(sorted(regions), range(R)))
    region_i = [region_id[r] for r in data.region]


    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma = mc.Uniform('sigma', lower=0, upper=1000, value=pl.ones(K))
    sigma_r = mc.Uniform('sigma_r', lower=0, upper=1000, value=pl.ones((R, K)))


    # predictions
    @mc.stochastic
    def u_r(value=pl.zeros((R, K)), sigma=sigma):
        return mc.normal_like(value, 0, pl.resize(sigma**-2., (R,K)))

    @mc.deterministic
    def mu(X=X, beta=beta, u_r=u_r):
        """ mu_i,r,c,t = (beta + u_r) * X"""
        effect = beta + u_r[region_i, :]
        return pl.sum(effect.T * X, axis=0)

    @mc.deterministic
    def tau(X=X, sigma_r=sigma_r, sigma_e=sigma_e):
        """ y_i - mu_i ~ beta * (N(0, sigma_r,k**2)) + N(0, sigma_i**2)"""
        var_u_i = (sigma_r**2.)[region_i, :]
        variance = pl.sum(var_u_i.T * X**2., axis=0) + sigma_e**2.

        return 1. / variance

    @mc.deterministic
    def predicted(mu=mu, tau=tau):
        return mc.rnormal(mu, tau)


    # likelihood
    i_obs = pl.find(1 - pl.isnan(data.y))
    @mc.observed
    def y(value=data.y, i_obs=i_obs, mu=mu, tau=tau):
        return mc.normal_like(value[i_obs], mu[i_obs], tau[i_obs])

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.u_r)
    return mod_mc


def gp_re(data):
    """ Gaussian Process Random Effect Model, where variation that is
    not explained by fixed effects model is modeled with GP::
    
        Y_r,c,t = beta * X_r,c,t + f_c(t) + e_r,c,t
        f_c(t) ~ GP(0, C)
        C ~ Matern(2, sigma_f, tau_f)
    """
    # covariates, etc
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma_f = mc.Uniform('sigma_f', lower=0, upper=1000, value=1)
    tau_f = mc.Uniform('tau_f', lower=0, upper=1000, value=1)

    # predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        return pl.dot(beta, X)

    # gaussian process and residuals (implements the likelihood)
    # need to organize residuals in panels to measure GP likelihood
    sm = {}  # sm stands for submodel, and maybe should be changed to something more descriptive
    obs = []
    pred = []
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]
        M = gp.Mean(lambda x: pl.zeros(len(x)))
        C = gp.Covariance(gp.matern.euclidean, amp=sigma_f, scale=tau_f, diff_degree=2)
        sm_c = gp.GPSubmodel('sm_%s'%c, M, C, mesh=data.year[i_c], init_vals=pl.zeros_like(i_c))
        sm[c] = sm_c
    
        i_c = [i for i in range(len(data)) if data.country[i] == c and not pl.isnan(data.y[i])]
        @mc.observed(name='obs_%s'%c)
        def obs_c(value=data.y[i_c], t=data.year[i_c], i_c=i_c, f_c=sm_c.f, mu=mu, sigma_e=sigma_e):
            return mc.normal_like(value, mu[i_c] + f_c(t), sigma_e**-2.)
        obs.append(obs_c)

        i_c = [i for i in range(len(data)) if data.country[i] == c]
        @mc.deterministic(name='pred_%s'%c)
        def pred_c(t=data.year[i_c], i_c=i_c, f_c=sm_c.f, mu=mu, sigma_e=sigma_e):
            y_rep_c = pl.zeros_like(data.y)
            y_rep_c[i_c] = mc.rnormal(mu[i_c] + f_c(t), sigma_e**-2.)
            return y_rep_c
            
        pred.append(pred_c)

    @mc.deterministic
    def predicted(pred=pred):
        return pl.sum(pred, axis=0)

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    # TODO: determine if this is the appropriate step method
    for f in [sm_c.f for sm_c in sm.values()]:
        mod_mc.use_step_method(mc.NoStepper, f)

    return mod_mc

# alternative implementation
def gp_re2(data):
    """ Gaussian Process Random Effect Model, where variation that is
    not explained by fixed effects model is modeled with GP::
    
        Y_r,c,t = beta * X_r,c,t + f_c(t) + e_r,c,t
        f_c(t) ~ GP(0, C)
        C ~ Matern(2, sigma_f, tau_f)
        
    Possibly more efficient than gp_re implementation.
    """
    # covariates, etc
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma_f = mc.Uniform('sigma_f', lower=0, upper=1000, value=1)
    tau_f = mc.Uniform('tau_f', lower=0, upper=1000, value=1)

    # predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        return pl.dot(beta, X)

    # gaussian process implemented implicitly in observation likelihood
    # need to organize observations into panels to calculate likelihood including GP
    index_dict = {}
    obs = []
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c and not pl.isnan(data.y[i])]
        index_dict[c] = i_c

        @mc.observed(name='obs_%s'%c)
        def obs_c(value=data.y[i_c], t=data.year[i_c], i_c=i_c, sigma_f=sigma_f, tau_f=tau_f, mu=mu, sigma_e=sigma_e):
            C_c = gp.matern.euclidean(t, t, amp=sigma_f, scale=tau_f, diff_degree=2)
            return mc.mv_normal_cov_like(value, mu[i_c], C_c + sigma_e**2. * pl.eye(len(i_c)))
        obs.append(obs_c)

    @mc.deterministic
    def predicted(data=data, mu=mu, sigma_f=sigma_f, tau_f=tau_f, sigma_e=sigma_e):
        y_rep = pl.zeros_like(data.y)
        for c, i_c in index_dict.items():
            t = data.year[i_c]
            C_c = gp.matern.euclidean(t, t, amp=sigma_f, scale=tau_f, diff_degree=2)
            y_rep[i_c] = mc.rmv_normal_cov(mu[i_c], C_c + sigma_e**2. * pl.eye(len(i_c)))
        return y_rep

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    return mod_mc


def nested_gp_re(data):
    """ Random Effect model, with country random effects nested in
    regions and gaussian process correlations in residuals::
    
        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + f_r(t) + f_c(t) + e_r,c,t
        u_r[k] ~ N(0, sigma_k^2)
        u_r,c,t[k] ~ N(0, sigma_r,k^2)

        f_r(t) ~ GP(0, C_1)
        f_c(t) ~ GP(0, C_2)  
        C_1 ~ Matern(2, sigma^1_f, tau^1_f)
        C_2 ~ Matern(2, sigma^2_f, tau^2_f)
    """
    # covariates, etc
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    regions = set(data.region)
    R = len(regions)
    region_id = dict(zip(sorted(regions), range(R)))
    region_i = [region_id[r] for r in data.region]


    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma = mc.Uniform('sigma', lower=0, upper=1000, value=pl.ones(K))
    sigma_r = mc.Uniform('sigma_r', lower=0, upper=1000, value=pl.ones((R, K)))

    sigma_f1 = mc.Uniform('sigma_f1', lower=0, upper=1000, value=10)
    tau_f1 = mc.Uniform('tau_f1', lower=0, upper=1000, value=1)
    sigma_f2 = mc.Uniform('sigma_f2', lower=0, upper=1000, value=1)
    tau_f2 = mc.Uniform('tau_f2', lower=0, upper=1000, value=1)

    # predictions
    @mc.stochastic
    def u_r(value=pl.zeros((R, K)), sigma=sigma):
        return mc.normal_like(value, 0, pl.resize(sigma**-2., (R,K)))

    @mc.deterministic
    def mu(X=X, beta=beta, u_r=u_r):
        """ mu_i,r,c,t = (beta + u_r) * X"""
        effect = beta + u_r[region_i, :]
        return pl.sum(effect.T * X, axis=0)

    @mc.deterministic
    def tau(X=X, sigma_r=sigma_r, sigma_e=sigma_e):
        """ y_i - mu_i - f_c(t) ~ beta * (N(0, sigma_r,k**2)) + N(0, sigma_i**2)"""
        var_u_i = (sigma_r**2.)[region_i, :]
        variance = pl.sum(var_u_i.T * X**2., axis=0) + sigma_e**2.

        return 1. / variance

    # nested gaussian processes and residuals to implements the likelihood
    f = {}
    year_start = 1990  # FIXME: don't hard code year_start or year_end
    year_end = 2005
    t = pl.arange(year_start, year_end)
    for r in set(data.region):
        @mc.stochastic(name='f_%s'%r)
        def f_r(value=pl.zeros_like(t), t=t, sigma_f1=sigma_f1, tau_f1=tau_f1):
            C = gp.matern.euclidean(t, t, amp=sigma_f1, scale=tau_f1, diff_degree=2)
            return mc.mv_normal_cov_like(value, pl.zeros_like(value), C)
        f[r] = f_r

    # need to organize residuals in panels to measure GP likelihood
    sm = {}
    res = []
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]
        M = gp.Mean(lambda x: pl.zeros(len(x)))
        C = gp.Covariance(gp.matern.euclidean, amp=sigma_f2, scale=tau_f2, diff_degree=2)
        sm_c = gp.GPSubmodel('sm_%s'%c, M, C, mesh=data.year[i_c], init_vals=pl.zeros_like(i_c))
        sm[c] = sm_c
        f_r = f[data.region[i_c[0]]]  # find the latent gp var for the region which contains this country
    
        # data likelihood represented as a potential
        i_c = [i for i in range(len(data)) if data.country[i] == c and not pl.isnan(data.y[i])]
        @mc.potential(name='residual_%s'%c)
        def res_c(data=data, i_c=i_c, f_r=f_r, f_c=sm_c.f, mu=mu, tau=tau):
            return mc.normal_like(data.y[i_c] - mu[i_c] - f_r[data.year[i_c]-year_start] - f_c(data.year[i_c]), 0, tau[i_c])
        res.append(res_c)

    @mc.deterministic
    def predicted(data=data, mu=mu, f=f, sm=sm, tau=tau):
        y_rep = pl.zeros_like(data.y)
        for i in range(len(data)):  # TODO: speed up by vector operations instead of loop
            country = data.country[i]
            region = data.region[i]
            year = data.year[i]
            y_rep[i] = mc.rnormal(mu[i] + f[region][year-year_start] + sm[country].f(year), tau[i])
        return y_rep

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.u_r)
    # TODO: determine if this is the appropriate step method
    for f in [sm_c.f for sm_c in sm.values()]:
        mod_mc.use_step_method(mc.NoStepper, f)
    return mod_mc

# alternative implementation
def nested_gp_re2(data):
    """ Random Effect model, with country random effects nested in
    regions and gaussian process correlations in residuals::
    
        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + f_r(t) + f_c(t) + e_r,c,t
        u_r[k] ~ N(0, sigma_k^2)
        u_r,c,t[k] ~ N(0, sigma_r,k^2)

        f_r(t) ~ GP(0, C_1)
        f_c(t) ~ GP(0, C_2)  
        C_1 ~ Matern(2, sigma^1_f, tau^1_f)
        C_2 ~ Matern(2, sigma^2_f, tau^2_f)
    """
    # covariates, etc
    K = len(data.dtype)-4 # number of covariates
    X = pl.array([data['x%d'%i] for i in range(K)])

    regions = set(data.region)
    R = len(regions)
    region_id = dict(zip(sorted(regions), range(R)))
    region_i = [region_id[r] for r in data.region]

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)
    sigma = mc.Uniform('sigma', lower=0, upper=1000, value=pl.ones(K))
    sigma_r = mc.Uniform('sigma_r', lower=0, upper=1000, value=pl.ones((R, K)))

    sigma_f1 = mc.Uniform('sigma_f1', lower=0, upper=1000, value=10)
    tau_f1 = mc.Uniform('tau_f1', lower=0, upper=1000, value=1)
    sigma_f2 = mc.Uniform('sigma_f2', lower=0, upper=1000, value=1)
    tau_f2 = mc.Uniform('tau_f2', lower=0, upper=1000, value=1)

    # predictions
    @mc.stochastic
    def u_r(value=pl.zeros((R, K)), sigma=sigma):
        return mc.normal_like(value, 0, pl.resize(sigma**-2., (R,K)))

    @mc.deterministic
    def mu(X=X, beta=beta, u_r=u_r):
        """ mu_i,r,c,t = (beta + u_r) * X"""
        effect = beta + u_r[region_i, :]
        return pl.sum(effect.T * X, axis=0)

    @mc.deterministic
    def re_var(X=X, sigma_r=sigma_r, sigma_e=sigma_e):
        """ y_i - mu_i - f_c(t) ~ beta * (N(0, sigma_r,k**2)) + N(0, sigma_i**2)"""
        var_u_i = (sigma_r**2.)[region_i, :]
        variance = pl.sum(var_u_i.T * X**2., axis=0) + sigma_e**2.

        return variance

    # nested gaussian processes and residuals to implements the likelihood
    sm = {}
    for r in set(data.region):
        i_r = [i for i in range(len(data)) if data.region[i] == r]
        M = gp.Mean(lambda x: pl.zeros(len(x)))
        C = gp.Covariance(gp.matern.euclidean, amp=sigma_f1, scale=tau_f1, diff_degree=2)
        sm_r = gp.GPSubmodel('sm_%s'%r, M, C, mesh=pl.unique(data.year[i_r]))
        sm[r] = sm_r

    # organize data in panels to measure likelihood
    obs = []
    pred = []
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c and not pl.isnan(data.y[i])]
        f_r = sm[data.region[i_c[0]]].f  # find the latent gp var for the region which contains this country
        @mc.observed(name='obs_%s'%c)
        def obs_c(value=data.y[i_c], t=data.year[i_c], f_r=f_r, i_c=i_c, sigma_f2=sigma_f2, tau_f2=tau_f2, mu=mu, re_var=re_var):
            C_c = gp.matern.euclidean(t, t, amp=sigma_f2, scale=tau_f2, diff_degree=2)
            return mc.mv_normal_cov_like(value, mu[i_c] + f_r(t), C_c + pl.diagflat(re_var[i_c]))
        obs.append(obs_c)

        i_c = [i for i in range(len(data)) if data.country[i] == c]
        @mc.deterministic(name='pred_%s'%c)
        def pred_c(t=data.year[i_c], f_r=f_r, i_c=i_c, sigma_f2=sigma_f2, tau_f2=tau_f2, mu=mu, re_var=re_var):
            y_rep_c = pl.zeros_like(data.y)
            C_c = gp.matern.euclidean(t, t, amp=sigma_f2, scale=tau_f2, diff_degree=2)
            y_rep_c[i_c] = mc.rmv_normal_cov(mu[i_c] + f_r(t), C_c + pl.diagflat(re_var[i_c]))
            return y_rep_c
        pred.append(pred_c)

    @mc.deterministic
    def predicted(pred=pred):
        return pl.sum(pred, axis=0)

    # set up MCMC step methods
    mod_mc = mc.MCMC(vars())
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.beta)
    mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.u_r)
    # TODO: determine if this is the appropriate step method
    for f in [sm_c.f for sm_c in sm.values()]:
        mod_mc.use_step_method(mc.NoStepper, f)
    return mod_mc


def run_all_models(data, testing=False):
    """ Run models for testing and comparison
    """
    mc_dict = {}
    for mod in [fe, re, nested_re, gp_re, gp_re2, nested_gp_re, nested_gp_re2]:
        print "setting up model (%s)" % mod
        mod_mc = mod(data)

        print "sampling with MCMC"
        if testing == True:
            mod_mc.sample(iter=2)
        else:
            mod_mc.sample(iter=20000, burn=10000, thin=100, verbose=1)

        print "saving results"
        mc_dict[mod] = mod_mc

    return mc_dict

def test():
    """ Test that the models all run, data generation works, and graphics functions work
    """
    print "testing data generating module"
    import data
    data.test()

    print "testing models"
    data = pl.csv2rec('missing_noisy_test_data.csv')
    mc_dict = run_all_models(data, testing=True) # test fitting models

    print "testing graphics"
    import graphics
    graphics.test(data, mc_dict)

if __name__ == '__main__':
    test()
