"""
Models are the object-oriented way to deal with the database.  Each
class in the model module corresponds to a database table, and each
instance of a class corresponds to a row

I'm putting doctests here until I figure out how to run them in the models
>>> probabilistic_utils.const_func([1,2,3], 17.0)
array([ 17., 17., 17.])

>>> M,C = probabilistic_utils.uninformative_prior_gp()
>>> M([0,10,100])
array([-10., -10., -10.])

>>> probabilistic_utils.gp_interpolate([0,100], [0,100], [0,100])
array([ 0., 100.])
"""

import probabilistic_utils
import django_utils

from disease import *
from region import *
#from study import *
#from table import *
from rate import *
from population import *
from age_specific_rate_function import *
from disease_model import *
#from jobs import *

