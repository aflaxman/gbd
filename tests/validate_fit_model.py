""" Test fitting models
"""

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.', '..']

import pylab as pl

import data
import data_simulation
import data_model
import consistent_model
import fit_model
reload(consistent_model)
reload(data_model)

## create model

model = data.ModelData()
# generate simulated data
n = 50
sigma_true = .025
a = pl.arange(0, 100, 1)
pi_age_true = .0001 * (a * (100. - a) + 100.)
model.input_data = data_simulation.simulated_age_intervals('p', n, a, pi_age_true, sigma_true)
model.input_data['data_type'][-1] = 'r'  # make sure that there are multiple data types in the data set


def validate_fit_data_model():
    vars = data_model.data_model('validation', model, 'p',
                                 root_area='all', root_sex='total', root_year='all',
                                 mu_age=None, mu_age_parent=None)

    m = fit_model.fit_data_model(vars)

    return m


def validate_fit_consistent_model():
    # create model and priors
    vars = consistent_model.consistent_model(model, 'all', 'total', 'all', {})

    m = fit_model.fit_consistent_model(vars)

    return m


if __name__ == '__main__':
    m1 = validate_fit_data_model()
    m2 = validate_fit_consistent_model()
    
