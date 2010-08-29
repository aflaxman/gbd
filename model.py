""" Functions to generate PyMC models of spatio-temporal health data
"""


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
    import pylab as pl
    import pymc as mc

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

    return mc.Model(vars())

def re():
    """ Random Effect model::
    
        Y_r,c,t = (beta + u_r,c,t) * X_r,c,t + e_r,c,t
        u_r,c,t[i] ~ N(0, sigma_i^2)
        e_r,c,t ~ N(0, sigma_e^2)

    Example
    -------
    >>> import pymc, models
    >>> re = models.re()
    >>> pymc.MCMC(re).sample(iter=20000, burn=10000, thin=10, verbose=1)
    >>> re.beta.stats()
    """
    import pylab as pl
    import pymc as mc

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

    return mc.Model(vars())
