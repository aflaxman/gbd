""" Test Data Model """

# matplotlib will open windows during testing unless you do the following
import matplotlib
matplotlib.use("AGG")

# add to path, to make importing possible
import sys
sys.path += ['.']

import pylab as pl
import pymc as mc

import data


class TestClass:
    def test_blank_model(self):
        d = data.Data()

        for field in 'data_type value area sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci age_weights'.split():
            assert field in d.input_data.columns, 'Input data CSV should have field "%s"' % field


if __name__ == '__main__':
    import nose
    nose.runmodule()
    
