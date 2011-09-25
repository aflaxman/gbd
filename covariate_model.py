""" Covariate models"""

import pylab as pl
import pymc as mc


def covariate_model(name, mu, X, beta):
    """ Generate PyMC objects covariate adjusted version of mu

    Parameters
    ----------
    name : str
    mu : the unadjusted mean parameter for this node
    X_obs : the covariate matrix for the data observations relevant to this node
    X_pred : the covariate matrix for the output predictions relevant to this node
    beta : the effect coefficients

    pi_pred_parent : adjusted mean for the node by parent in hierarchy, optional
    weight : similarity weight between this node and its parent, optional

    Results
    -------
    Returns dict of PyMC objects, including 'pi', the covariate
    adjusted predicted values for the mu and X provided

    Notes
    -----
    This is used twice, first to predict the means of the observed
    data, and then to predict the level values for the node
    """

    @mc.deterministic(name='pi_%s'%name)
    def pi(mu=mu, X=X, beta=beta):
        return mu * pl.exp(pl.dot(X, beta))
    return dict(pi=pi)
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
