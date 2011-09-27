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
    tau = (sigma_difference*10.)**-2.  # try multiplying similarity by 10, since there are 100 points that are not independent
    @mc.potential(name='pi_similarity_%s'%name)
    def pi_sim(pi_child=pi_child, pi_parent=pi_parent, tau=tau):
        return mc.normal_like(pl.log(pi_child.clip(offset, pl.inf)), pl.log(pi_parent.clip(offset, pl.inf)), tau)

    return dict(pi_sim=pi_sim, pi_parent=pi_parent)
