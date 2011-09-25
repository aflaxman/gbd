""" Several rate models"""

import pylab as pl
import pymc as mc


def neg_binom_model(name, pi, delta, p, n):
    """ Generate PyMC objects for a Negative Binomial model

    Parameters
    ----------
    pi : pymc.Node, expected values of rates
    delta : pymc.Node, dispersion parameters of rates
    p : array, observed values of rates
    n : array, effective sample sizes of rates

    Results
    -------
    Returns dict of PyMC objects, including 'p_obs' and 'p_pred'
    the observed stochastic likelihood and data predicted stochastic
    """
    @mc.observed(name='p_obs_%s'%name)
    def p_obs(value=p*n, pi=pi, delta=delta, n=n):
        return mc.negative_binomial_like(value, pi*n, delta)

    @mc.deterministic(name='p_pred_%s')
    def p_pred(pi=pi, delta=delta, n=n):
        return mc.rnegative_binomial(pi*n, delta) / pl.array(n, dtype=float)

    return dict(p_obs=p_obs, p_pred=p_pred)

def normal_model():
    pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
    sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

    @mc.potential
    def obs(pi=pi, sigma=sigma):
        return mc.normal_like(r, pi, 1./(s**2 + sigma**2))

    @mc.deterministic
    def pred(pi=pi, sigma=sigma):
        s_pred = pl.sqrt(pi*(1-pi)/n_pred)
        return mc.rnormal(pi, 1./(s_pred + sigma))

    ### @export 'normal-fit-and-store'
    mc.MCMC([pi, sigma, obs, pred]).sample(iter, burn, thin)

    results['Normal'] = dict(pi=pi.stats(), pred=pred.stats())


def log_normal_model():
    pi = mc.Uniform('pi', lower=0, upper=1, value=.5)
    sigma = mc.Uniform('sigma', lower=0, upper=10, value=.01)

    @mc.potential
    def obs(pi=pi, sigma=sigma):
        return mc.normal_like(pl.log(r), pl.log(pi), 1./((s/r)**2 + sigma**2))

    pred_s = pl.sqrt(r * (1-r) / n_pred)
    @mc.deterministic
    def pred(pi=pi, sigma=sigma):
        s_pred = pl.sqrt(pi*(1-pi)/n_pred)
        return pl.exp(mc.rnormal(pl.log(pi), 1./((s_pred/pi)**2 + sigma**2)))

    ### @export 'log-normal-fit-and-store'
    mc.MCMC([pi, sigma, obs, pred]).sample(iter, burn, thin)

    results['Lognormal'] = dict(pi=pi.stats(), pred=pred.stats())

