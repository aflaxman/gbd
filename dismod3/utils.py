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
               'case-fatality data', 'all-cause mortality data',
               'duration data',
               ]

output_data_types = ['Incidence', 'Remission', 'Case-fatality',
                     'Prevalence', 'Duration']

stoch_var_types = output_data_types + ['bins']

gbd_regions = [ u'Asia Pacific, High Income',
                u'Asia, Central',
                u'Asia, East',
                u'Asia, South',
                u'Asia, Southeast',
                u'Australasia',
                u'Caribbean',
                u'Europe, Central',
                u'Europe, Eastern',
                u'Europe, Western',
                u'Latin America, Andean',
                u'Latin America, Central',
                u'Latin America, Southern',
                u'Latin America, Tropical',
                u'North Africa/Middle East',
                u'North America, High Income',
                u'Sub-Saharan Africa, Central',
                u'Sub-Saharan Africa, East',
                u'Sub-Saharan Africa, Southern',
                u'Sub-Saharan Africa, West']



gbd_years = [ 1990, 2005 ]
gbd_sexes = [ 'male', 'female']


def gbd_keys(type_list=stoch_var_types,
             region_list=gbd_regions,
             year_list=gbd_years,
             sex_list=gbd_sexes):
    """ Make a list of gbd keys for the type, region, year, and sex
    specified

    Parameters
    ----------
    type_list : list, optional, subset of ['incidence', 'remission', 'case-fatality']
    region_list : list, optional, subset of 21 GBD regions
    year_list : list, optional, subset of ['1990', '2005']
    sex_list : list, optional, subset of ['male', 'female']

    Results
    -------
    A list of gbd keys corresponding to all combinations of list
    items.
    """
    key_list = []

    # special case: prevalence is controlled by incidence, remission,
    # and case-fatality
    if type_list == [ 'prevalence' ]:
        types = [clean(t) for t in output_data_types]
        
    for t in type_list:
        for r in region_list:
            for y in year_list:
                for s in sex_list:
                    key_list.append(gbd_key_for(t, r, y, s))
    return key_list

def clean(str):
    """ Return a 'clean' version of a string, suitable for using as a hash
    string or a class attribute.
    """
    str = str.strip()
    str = str.lower()
    str = str.replace(',', '')
    str = str.replace('/', '_')
    str = str.replace(' ', '_')
    return str

KEY_DELIM_CHAR = '+'

def gbd_key_for(type, region, year, sex):
    """ Make a human-readable string that can be used as a key for
    storing estimates for the given type/region/year/sex.
    """

    return KEY_DELIM_CHAR.join([clean(type), clean(region),
                                str(year), clean(sex)])
    

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






