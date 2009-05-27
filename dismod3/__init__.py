"""
This module provide the python shell interface to the generic disease
modeling system, as well as all of the parameters and functions
specific to the statistical modeling and plotting of generic disease.

see ``../docs/tutorial.rst`` for more details on the interface.
"""

from utils import get_disease_model, post_disease_model
from utils import gbd_regions, gbd_years, gbd_sexes, data_types, gbd_key_for

from plotting import plot_disease_model, sparkplot_disease_model, sparkplot_boxes, overlay_plot_disease_model
from bayesian_models.probabilistic_utils import *
