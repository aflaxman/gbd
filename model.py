""" Functions to generate PyMC models of spatio-temporal health data
"""

import pylab as pl
import pymc as mc
import pymc.gp as gp

def fe():
    """ Fixed Effect model::
    
        Y_r,c,t = beta * X_r,c,t + e_r,c,t
        e_r,c,t ~ N(0, sigma^2)

    Example
    -------
    >>> import pymc, models
    >>> fe = models.fe()
    >>> pymc.MCMC(fe).sample(iter=20000, burn=10000, thin=10)
    >>> fe.beta.stats()
    """

    data = pl.csv2rec('data.csv')

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
    @mc.observed
    def y(value=data.y, mu=mu, sigma=sigma):
        return mc.normal_like(value, mu, sigma**-2.)

    return vars()


def re():
    """ Random Effect model::
    
        Y_r,c,t = (beta + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r,c,t[i] ~ N(0, sigma_i^2)
        e_r,c,t ~ N(0, sigma_e^2)
    """

    data = pl.csv2rec('data.csv')

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
    @mc.observed
    def y(value=data.y, mu=mu, tau=tau):
        return mc.normal_like(value, mu, tau)

    return vars()


def nested_re():
    """ Random Effect model, with country random effects nested in
    regions::
    
        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r[k] ~ N(0, sigma_k^2)
        u_r,c,t[k] ~ N(0, sigma_r,k^2)
        e_r,c,t ~ N(0, sigma_e^2)

    """

    data = pl.csv2rec('data.csv')


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
    @mc.observed
    def y(value=data.y, mu=mu, tau=tau):
        return mc.normal_like(value, mu, tau)

    return vars()


def gp_re():
    """ Gaussian Process Random Effect Model, where variation that is
    not explained by fixed effects model is modeled with GP::
    
        Y_r,c,t = beta * X_r,c,t + f_c(t) + e_r,c,t
        f_c(t) ~ GP(0, C)
        C ~ Matern(2, sigma_f, tau_f)
    """

    data = pl.csv2rec('data.csv')

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
    f = []
    res = []
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]
        M = gp.Mean(lambda x: pl.zeros(len(x)))
        C = gp.Covariance(gp.matern.euclidean, amp=sigma_f, scale=tau_f, diff_degree=2)
        f_c = gp.GPSubmodel('f_%s'%c, M, C, mesh=data.year[i_c])
        f.append(f_c)
    
        @mc.potential(name='residual_%s'%c)
        def res_c(data=data, i_c=i_c, f_ct=f_c.f_eval, mu=mu, sigma_e=sigma_e):
            return mc.normal_like(data.y[i_c] - mu[i_c] - f_ct, 0, sigma_e**-2.)
        res.append(res_c)

    return vars()

def nested_gp_re():
    """ Random Effect model, with country random effects nested in
    regions and gaussian process correlations in residuals::
    
        Y_r,c,t = (beta + u_r + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r[k] ~ N(0, sigma_k^2)
        u_r,c,t[k] ~ N(0, sigma_r,k^2)

        e_r,c,t ~ f_c(t)
        f_c(t) ~ GP(0, C)
        C ~ Matern(2, sigma_f, tau_f)
    """

    data = pl.csv2rec('data.csv')


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
    model_vars = dict(beta=beta, sigma_e=sigma_e, sigma=sigma, sigma_r=sigma_r)

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

    model_vars.update(u_r=u_r, mu=mu, tau=tau)

    # gaussian process and residuals (implements the likelihood)
    # need to organize residuals in panels to measure GP likelihood
    res = {}
    f = {}
    for c in set(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]
        M = gp.Mean(lambda x: pl.zeros(len(x)))
        C = gp.Covariance(gp.matern.euclidean, amp=1, scale=15, diff_degree=2)
        f_c = gp.GaussianProcess('f_%s'%c, gp.GPSubmodel('f_c', M, C, mesh=data.year[i_c]))
        f[c] = f_c
    
        @mc.potential(name='residual_%s'%c)
        def res_c(data=data, i_c=i_c, f_c=f_c, mu=mu, tau=tau):
            return mc.normal_like(data.y[i_c] - mu[i_c] - f_c(data.year[i_c]), 0, tau[i_c])
        res[c] = res_c

    #model_vars.update(res=res, f=f)
    return model_vars


def run_all_models():
    mc_dict = {}
    for mod in [fe, re, nested_re, gp_re]:
        mod_vars = mod()
        mod_mc = mc.MCMC(mod_vars)

        if mod == nested_re:
            mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_vars['u_r'])

        mod_mc.sample(iter=20000, burn=10000, thin=10, verbose=1)

        mc_dict[mod] = mod_mc

    return mc_dict
