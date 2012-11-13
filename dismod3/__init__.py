"""
This module provide the python shell interface to the generic disease
modeling system, as well as all of the parameters and functions
specific to the statistical modeling and plotting of generic disease.

see ``../docs/tutorial.rst`` for more details on the interface.
"""
import settings
reload(settings)

from settings import gbd_regions, gbd_years, gbd_sexes

from disease_json import get_job_queue, remove_from_job_queue, \
    try_posting_disease_model, load_disease_model, add_covariates_to_disease_model

import neg_binom_model, normal_model, log_normal_model, generic_disease_model, gbd_disease_model, regional_similarity_matrices
import plotting, table, utils

import age_integrating_model as age_group
import age_pattern
import covariate_model as covariates
import rate_model
import data
import data_model
import fit_model
import graphics

import ism
import fit

reload(data)
reload(data_model)
reload(age_pattern)
reload(covariates)
reload(rate_model)
reload(age_group)
reload(ism)
reload(fit)
reload(graphics)

import tests
