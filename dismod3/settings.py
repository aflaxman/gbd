
# server names, etc.
from server_settings import *

# disease model parameters
NEARLY_ZERO = 1.e-10
MAX_AGE = 101

MISSING = -99

PRIOR_SEP_STR = ','

KEY_DELIM_CHAR = '+'

data_types = ['prevalence data',
              'incidence data',
              'remission data',
              'case-fatality data',
              'duration data',
              'all-cause mortality data',
              ]

output_data_types = ['Prevalence',
                     'Incidence',
                     'Remission',
                     'Case-fatality',
                     'Relative-risk',
                     'Duration',
                     'Incidence x Duration']

stoch_var_types = output_data_types + ['bins']

gbd_regions = [u'Asia Pacific, High Income',
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
               u'Oceania',
               u'Sub-Saharan Africa, Central',
               u'Sub-Saharan Africa, East',
               u'Sub-Saharan Africa, Southern',
               u'Sub-Saharan Africa, West']

gbd_years = ['1990', '2005']

gbd_sexes = ['Male', 'Female']
