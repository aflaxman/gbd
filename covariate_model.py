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

    # make U and alpha
    p_U = 2 + hierarchy.number_of_nodes()  # random effects for sex, time, area
    U = pandas.DataFrame(pl.zeros((n, p_U)), columns=['sex', 'time'] + hierarchy.nodes())
    for i, row in data.T.iteritems():
        U.ix[i, 'sex'] = float(data.ix[i, 'sex'] == 'male')
        U.ix[i, 'time'] = .5 * (data.ix[i, 'year_start'] + data.ix[i, 'year_end']) - 2000.
        for node in nx.shortest_path(hierarchy, root, data.ix[i, 'area']):
            U.ix[i, node] = 1.
    U = U.select(lambda col: U[col].std() > 0, 1)  # drop blank columns
    if len(U.columns) > 0:
        alpha = mc.Normal(name='alpha_%s'%name, mu=0, tau=1., value=pl.zeros_like(U.columns))  # TODO: put hyper-prior on tau
    else:
        alpha = pl.array([])

    # TODO: consider faster ways to calculate dot(U, alpha), since the matrix is sparse and (half-)integral

    X = data.select(lambda col: col.startswith('x_'), axis=1)
    if len(X.columns) > 0:
        beta = mc.Uninformative('beta', value=pl.zeros_like(X.columns))
    else:
        beta = pl.array([])

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, U=U, alpha=alpha, X=X, beta=beta):
        return mu * pl.exp(pl.dot(U, alpha) + pl.dot(X, beta))

    return dict(pi=pi, U=U, alpha=alpha, X=X, beta=beta)

""" Stuff for similarity potentials
    if X_parent != None:
        @mc.potential(name='pi_similarity_%s'%name)
        def pi_sim(X1=X, X2=X_parent, beta=beta, tau=weight**-2.):
            return mc.normal_like(pl.dot(X1-X2, beta), 0, tau)
        vars['pi_sim'] = pi_sim

    if mu_parent != None:
        @mc.potential(name='pi_similarity_%s'%name)
        def pi_sim(X1=X, X2=X_parent, beta=beta, tau=weight**-2.):
            return mc.normal_like(pl.dot(X1-X2, beta), 0, tau)
        vars['pi_sim'] = pi_sim

    return vars
"""
