import pylab as pl
import numpy as np
import pymc as mc
from pymc import gp

from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import \
    trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE, MISSING

from disease_json import *
from model_utils import *

data_types = [ 'incidence data', 'prevalence data', 'remission data',
               'case-fatality data', 'all-cause mortality data', 'duration data',
               ]
gbd_regions = ['Australasia',
               'North America',
               'Asia, East',
               ]
gbd_years = [ 1995, 2005 ]
gbd_sexes = [ 'male', 'female']

def clean(str):
    """ Return a 'clean' version of a string, suitable for using as a hash
    string or a class attribute.
    """
    
    return str.strip().lower().replace(',', '').replace(' ', '_')

def gbd_key_for(type, region, year, sex):
    """ Make a human-readable string that can be used as a key for
    storing estimates for the given type/region/year/sex.
    """

    return '%s_%s_%s_%s' % (clean(type), clean(region), year, sex)
    

def fit(dm_id, probabilistic_model):
    """ Estimate disease parameters using a bayesian model

    Parameters
    ----------
    dm_id : int
      An id number for a disease model on the dismod server
    probabilistic_model : module, optional
      A python module that can do all the things a probabilistic model
      must do, Default is dismod3.generic_disease_model

    Returns
    -------
    dm : disease_json
      A thin wrapper around the json object returned by the dismod
      server

    Notes
    -----
    The probabilistic_model should be refactored to make it more
    OOPsy, and this function might need more features to become the
    workhorse of dismod analysis
    """
    dm = get_disease_model(dm_id)

    # filter out all data with type != data_type
    # dm.data = dm.filter_data(data_type=data_type)

    # store the probabilistic model code for future reference
    dm.set_model_source(probabilistic_model)
    dm.set_param_age_mesh([0.0, 0.5, 3.0, 10.0, 20.0, 30.0, 40.0,
                           50.0, 60.0, 70.0, 80.0, 90.0, 100.0])

    dm.set_estimate_age_mesh(range(MAX_AGE))

    probabilistic_model.initialize(dm)
    vars = probabilistic_model.setup(dm)
    map = probabilistic_model.map_fit(dm, vars)
    mcmc = probabilistic_model.mcmc_fit(dm, vars)

    url = post_disease_model(dm)
    print 'url for fit:\n\t%s' % url

    return {'vars': vars,
            'map': map,
            'mcmc': mcmc,
            'disease_model': dm}






