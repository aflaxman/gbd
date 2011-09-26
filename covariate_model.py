""" Covariate models"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

def mean_covariate_model(name, mu, data, hierarchy, root):
    """ Generate PyMC objects covariate adjusted version of mu

    Parameters
    ----------
    name : str
    mu : the unadjusted mean parameter for this node
    beta : effect coefficient prior
    data : pandas.DataFrame containing design matrix for covariates,
    as well as columns for year_start, year_end, and area
    hierarchy : nx.DiGraph encoding hierarchical structure and similarity weights

    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for the mu and X provided

    Notes
    -----
    This is used twice, first to predict the means of the observed
    data, and then to predict the level values for the node
    """
    n = len(data.index)

    # make U
    p_U = 2 + hierarchy.number_of_nodes()  # random effects for sex, time, area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=['sex', 'time'] + hierarchy.nodes())
    for i, row in data.T.iteritems():
        U.ix[i, 'sex'] = float(data.ix[i, 'sex'] == 'male')
        U.ix[i, 'time'] = .5 * (data.ix[i, 'year_start'] + data.ix[i, 'year_end']) - 2000.
        for node in nx.shortest_path(hierarchy, root, data.ix[i, 'area']):
            U.ix[i, node] = 1.
    U = U.select(lambda col: U[col].std() > 0, 1)  # drop blank columns

    # make tau_alpha and alpha
    if len(U.columns) > 0:
        tau_alpha = mc.InverseGamma(name='tau_alpha_%s', alpha=.1, beta=.1, value=pl.ones_like(U.columns))
        alpha = mc.Normal(name='alpha_%s'%name, mu=0, tau=tau_alpha, value=pl.zeros_like(U.columns))  
    else:
        tau_alpha = pl.array([])
        alpha = pl.array([])

    # TODO: consider faster ways to calculate dot(U, alpha), since the matrix is sparse and (half-)integral

    X = data.select(lambda col: col.startswith('x_'), axis=1)
    X = X.select(lambda col: X[col].std() > 0, 1)  # drop blank columns
    if len(X.columns) > 0:
        beta = mc.Uniform('beta', -5., 5., value=pl.zeros_like(X.columns))
    else:
        beta = pl.array([])

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=U, alpha=alpha, X=X, beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, tau_alpha=tau_alpha, alpha=alpha, X=X, beta=beta)

def dispersion_covariate_model(name, data):

    eta=mc.Normal('eta_%s'%name, mu=5., tau=1., value=5.)

    Z = data.select(lambda col: col.startswith('z_'), axis=1)
    Z = Z.select(lambda col: Z[col].std() > 0, 1)  # drop blank cols
    if len(Z.columns) > 0:
        zeta = mc.Uniform('zeta', -5, 5, value=pl.zeros_like(Z.columns))
    else:
        zeta = pl.array([])

    @mc.deterministic(name='delta_%s'%name)
    def delta(eta=eta, zeta=zeta, Z=Z):
        return (50. + pl.exp(eta)) * pl.exp(pl.dot(Z, zeta))

    return dict(eta=eta, Z=Z, zeta=zeta, delta=delta)
