""" Similarity prior model"""

import pylab as pl
import pymc as mc
import pandas
import networkx as nx

def similar(name, pi_child, pi_parent, sigma_difference, offset=1.e-9):
    """ Generate PyMC objects encoding a simliarity prior on pi_child
    to pi_parent

    Parameters
    ----------
    name : str
    pi : array or pymc.Node, the predicted values for the child node in the hierarchy
    pi_parent : array or pymc.Node, the predicted values for the parent node in the hierarchy
    sigma_difference : float, the prior on how similar between child is to parent

    Results
    -------
    Returns dict of PyMC objects, including 'pi_sim', the similarity potential
    """
    @mc.potential(name='pi_similarity_%s'%name)
    def pi_sim(pi_child=pi_child, pi_parent=pi_parent, tau=sigma_difference**-2.):
        return mc.normal_like(pl.log(pl.maximum(pi_child,offset)) - pl.log(pl.maximum(pi_parent,offset)), 0, tau)

    return dict(pi_sim=pi_sim)
