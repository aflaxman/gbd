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
    K1 = count_covariates(data, 'x')
    X = pl.array([data['x%d'%i] for i in range(K1)])

    K2 = count_covariates(data, 'w')
    W = pl.array([data['w%d'%i] for i in range(K1)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K1))
    gamma = mc.Uninformative('gamma', value=pl.zeros(K2))
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
    
        Y_r,c,t,a = beta * X_r,c,t,a + f_r(t) + g_r(a) + h_c(t) + e_r,c,t,a

        f_r(t) ~ GP(0, C_rt)
        g_r(a) ~ GP(0, C_ra)
        h_c(t) ~ GP(0, C_ca)

        C_rt ~ Matern(2, sigma_f,0, tau_f,0)
        C_ra ~ Matern(2, sigma_f,1, tau_f,1)
        C_ct ~ Matern(2, sigma_f,2, tau_f,2)

        e_r,c,t,a ~ N(gamma * W_r,c,t,a, sigma_e)
    """
    # covariates
    K1 = count_covariates(data, 'x')
    X = pl.array([data['x%d'%i] for i in range(K1)])

    K2 = count_covariates(data, 'w')
    W = pl.array([data['w%d'%i] for i in range(K1)])

    # priors
    beta = mc.Uninformative('beta', value=pl.zeros(K1))
    gamma = mc.Uninformative('gamma', value=pl.zeros(K2))
    sigma_e = mc.Uniform('sigma_e', lower=0, upper=1000, value=1)

    # semi-uninformative priors
    sigma_f = mc.Gamma('sigma_f', alpha=.05, beta=.1, value=1.*pl.ones(3))
    tau_f = mc.Gamma('tau_f', alpha=10., beta=.1, value=10.*pl.ones(3))


    # fixed-effect predictions
    @mc.deterministic
    def mu(X=X, beta=beta):
        """ mu_i,r,c,t,a = beta * X_i,r,c,t,a"""
        return pl.dot(beta, X)
    @mc.deterministic
    def sigma_explained(W=W, gamma=gamma):
        """ sigma_explained_i,r,c,t,a = gamma * W_i,r,c,t,a"""
        return pl.dot(gamma, W)


    # setup gaussian processes random effects to model additional variation
    # make index dict to convert from region/country/age to array index
    regions = pl.unique(data.region)
    countries = pl.unique(data.country)
    years = pl.unique(data.year)
    ages = pl.unique(data.age)

    r_index = dict([(r, i) for i, r in enumerate(regions)])
    c_index = dict([(c, i) for i, c in enumerate(countries)])
    t_index = dict([(t, i) for i, t in enumerate(years)])
    a_index = dict([(a, i) for i, a in enumerate(ages)])

    C = []
    for i, grid in enumerate([years, ages, years]):
        @mc.deterministic(name='C_%d'%i)
        def C_i(i=i, grid=grid, sigma_f=sigma_f, tau_f=tau_f, diff_degree=2.):
            return gp.matern.euclidean(grid, grid, amp=sigma_f[i], scale=tau_f[i], diff_degree=diff_degree)
        C.append(C_i)

    f = [mc.MvNormalCov('f_%s'%r, pl.zeros_like(years), C[0], value=pl.zeros_like(years)) for r in regions]
    g = [mc.MvNormalCov('g_%s'%r, pl.zeros_like(ages), C[1], value=pl.zeros_like(ages)) for r in regions]
    h = [mc.MvNormalCov('h_%s'%c, pl.zeros_like(years), C[2], value=pl.zeros_like(years)) for c in countries]

    # organize observations into country panels and calculate predicted value before sampling error
    country_param_pred = []
    country_tau = []
    for c in pl.unique(data.country):
        i_c = [i for i in range(len(data)) if data.country[i] == c]

        # find the index for this region, country, and for the relevant ages and times
        region = data.region[i_c[0]]
        r_index_c = r_index[region]
        c_index_c = c_index[c]
        t_index_c = [t_index[data.year[i]] for i in i_c]
        a_index_c = [a_index[data.age[i]] for i in i_c]
    
        # find predicted parameter value for all observations of country c
        @mc.deterministic(name='country_param_pred_%s'%c)
        def country_param_pred_c(i=i_c, mu=mu, f=f[r_index_c], g=g[r_index_c], h=h[c_index_c], a=a_index_c, t=t_index_c):
            """ country_param_pred_c[row] = parameter_predicted[row] * 1[row.country == c]"""
            country_param_pred_c = pl.zeros_like(data.y)
            country_param_pred_c[i] = mu[i] + f[t] + g[a] + h[t]
            return country_param_pred_c
        country_param_pred.append(country_param_pred_c)

        # find predicted parameter precision for all observations of country c
        @mc.deterministic(name='country_tau_%s'%c)
        def country_tau_c(i=i_c, sigma_explained=sigma_explained, sigma_e=sigma_e, var_d_c=data.se[i_c]**2.):
            """ country_tau_c[row] = tau[row] * 1[row.country == c]"""
            country_tau_c = pl.zeros_like(data.y)
            country_tau_c[i] = 1 / (sigma_e**2. + sigma_explained[i]**2. + var_d_c)
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
    
    # use covariance matrix to seed adaptive metropolis steps
    for r in range(len(regions)):
        mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.f[r], cov=pl.array(C[0].value*.01))
        mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.g[r], cov=pl.array(C[1].value*.01))
    for c in range(len(countries)):
        mod_mc.use_step_method(mc.AdaptiveMetropolis, mod_mc.h[c], cov=pl.array(C[2].value*.01))

    # find good initial conditions with MAP approx
    try:
        for var_list in [[mod_mc.obs, mod_mc.beta]] + \
            [[mod_mc.obs, f_r] for f_r in mod_mc.f] + \
            [[mod_mc.obs, g_r] for g_r in mod_mc.g] + \
            [[mod_mc.obs, h_c] for h_c in mod_mc.h] + \
            [[mod_mc.obs, mod_mc.beta, mod_mc.sigma_e]] + \
            [[mod_mc.obs, mod_mc.beta] + mod_mc.f] + \
            [[mod_mc.obs, mod_mc.beta] + mod_mc.g] + \
            [[mod_mc.obs, mod_mc.beta] + mod_mc.f + mod_mc.g] + \
            [[mod_mc.obs, h_c] for h_c in mod_mc.h]:
            print 'attempting to maximize likelihood of %s' % [v.__name__ for v in var_list]
            mc.MAP(var_list).fit(method='fmin_powell', verbose=1)
    except mc.ZeroProbability, e:
        print 'Warning: Optimization became infeasible:\n', e
    except KeyboardInterrupt:
        print 'Warning: Likelihood maximization canceled'
                                    
    return mod_mc


# helper functions
def count_covariates(data, var_name='x'):
    return len([n for n in data.dtype.names if n.startswith(var_name)])

