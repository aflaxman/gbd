"""
This module provide the python shell interface to the generic disease
modeling system, as well as all of the parameters and functions
specific to the statistical modeling and plotting of generic disease.

see ``../docs/tutorial.rst`` for more details on the interface.
"""

from settings import gbd_regions, gbd_years, gbd_sexes

from disease_json import get_job_queue, remove_from_job_queue, \
    try_posting_disease_model, load_disease_model, add_covariates_to_disease_model

import neg_binom_model, generic_disease_model, gbd_disease_model, plotting, utils
