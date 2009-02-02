"""
Views handle all interactions between a browser and the database.  I'm
going to break them out into little files, one view to control each
model.  Maybe I will regret this decision.
"""

import plotting_utils

from population_view import *
from rate_view import *
from age_specific_rate_function_view import *
