""" Similarity prior model"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

def similar(name, mu_child, mu_parent, sigma_parent, sigma_difference, offset=1.e-9):
    """ Generate PyMC objects encoding a simliarity prior on mu_child
    to mu_parent

    Parameters
    ----------
    name : str
    mu_child : array or pymc.Node, the predicted values for the child node in the hierarchy
    mu_parent : array, the predicted values for the parent node in the hierarchy
    sigma_parent : array, the predicted standard devation for the parent node in the hierarchy
    sigma_difference : float, the prior on how similar between child is to parent
    offest : float, an offset to deal with log of zero

    Results
    -------
    Returns dict of PyMC objects, including parent_mu_age and parent_sim the similarity potential
    """
    if isinstance(mu_parent, mc.Node):
        tau = 1. / (len(mu_child.value) * sigma_difference**2)
    else:
        tau = 1. / (len(mu_child.value) * (((sigma_parent+offset)/(mu_parent+offset)).clip(0., sigma_difference)**2))

    @mc.potential(name='parent_similarity_%s'%name)
    def parent_similarity(mu_child=mu_child, mu_parent=mu_parent,
                          tau=tau):
        log_mu_child = pl.log(mu_child.clip(offset, pl.inf))
        log_mu_parent = pl.log(mu_parent.clip(offset, pl.inf))
        return mc.normal_like(log_mu_child, log_mu_parent, tau)

    return dict(parent_similarity=parent_similarity)
