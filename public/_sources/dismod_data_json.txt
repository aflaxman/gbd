Specification of the DisMod III Dataset JSON format
--------------------------

This document describes the fields in the JSON object of a DisMod
Dataset, which can be served and stored in the DisMod III Data Server
at:

* http://winthrop.gs.washington.edu:5432/new/dm/<dataset_id>
* http://winthrop.gs.washington.edu:5432/new/dm/new

::

    dismod_dataset = {
      'params' : params (required), see below
      'data' : data_list (required), see below
    }

    params = {
      'id' : int (required), unique id of this dataset,
      'region' : str, probably one of 21 GBD regions, or 'World'
      'year' : str, probably one of '1995', '2005'
      'param_age_mesh' : [ float, float, ... ] (required)
      'estimate_age_mesh' : [ float, float, ... ] (required)
      'sex' : str (required), one of 'male', 'female', 'total'
      'condition' : str (required)

      'units' : units_hash (required), see below
      'priors' : prior_hash (optional), see below
      'estimate_type' : str, optional, one of 'fit each region/year/sex individually', 'borrow strength within regions', 'borrow strength across regions'

      'initial_value' : value_hash (optional), see below
      'map' : value_hash (optional), see below
      'mcmc_median' : value_hash (optional), see below
      'mcmc_mean' : value_hash (optional), see below
      'mcmc_lower_ui' : value_hash (optional), see below
      'mcmc_upper_ui' : value_hash (optional), see below
    }

    units_hash = { data_type_1 : str (required),
                   data_type_2 : str (optional),
                   ...
                 }
    prior_hash = { data_type_1 : prior_str (required), see below
                   data_type_2 : prior_str (optional),
                   ...
                 }
    prior_str = a special string that specifies the priors for estimating data of this data_type
    value_hash = { data_type_1 : [ float, float, ... ] (required), list length equals length of estimate_age_mesh
                   data_type_2 : [ float, float, ... ] (optional),
                   ...
                 }
      
    data_list = [ data_1, data_2, ... ]
    data_i = { 'id' : int (required), unique id
               'condition' : string (required)
               'gbd_cause' : str (required)
               'data_type' : str (required), one of the following types
                             'incidence data', 'prevalence data', 'remission data',
                             'case-fatality data', 'all-cause mortality data', 'duration data'

               'region' : str (required)
               'gbd_region' : str (required)
               'country' : str (optional)

               'sex' : str (required), one of 'male', 'female', 'total'

               'age_start' : int (required)
               'age_end' : int (required)

               'age_weights' : [ float, float, ... ] (optional), length equals age_end - age_start + 1,
                               default/missing assume to be [ 1, ... ]

               'year_start' : int (required)
               'year_end' : int (required)

               'value' : float (required), -99 means missing
               'standard_error' : float (required), -99 means missing
               'radix' : float (required)

               'citation' : str (optional)
               additional keys, with corresponding strs (optional)
            }
