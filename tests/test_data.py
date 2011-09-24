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
reload(data)

class TestClass:
    def test_blank_input_data(self):
        d = data.ModelData()

        for field in 'data_type value area sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci age_weights'.split():
            assert field in d.input_data.columns, 'Input data CSV should have field "%s"' % field

        for field in 'data_type area sex age_start age_end year_start year_end age_weights'.split():
            assert field in d.output_template.columns, 'Output template CSV should have field "%s"' % field

        for data_type in 'i p r f rr X'.split():
            assert data_type in d.parameters, 'Parameter dict should have entry for "%s"' % data_type

        assert d.areas_hierarchy.number_of_nodes() > 0, 'Areas hierarchy should be non-empty'

        assert len(d.areas_to_fit) > 0, 'Areas to fit should be non-empty'
        
    def test_from_gbd_json(self):
        d = data.ModelData.from_gbd_json('tests/dismoditis.json')

        assert len(d.input_data) == 17, 'dismoditis model has 17 data points'
        for field in 'data_type value area sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci age_weights'.split():
            assert field in d.input_data.columns, 'Input data CSV should have field "%s"' % field

        assert len(d.output_template) == 129168
        for field in 'data_type area sex age_start age_end year_start year_end age_weights'.split():
            assert field in d.output_template.columns, 'Output template CSV should have field "%s"' % field

        # TODO: test that expert priors were copied correctly

        assert d.areas_hierarchy.successors('asia_east') == ['MAC', 'PRK', 'TWN', 'HKG', 'CHN']
        assert d.areas_to_fit == 22

        
if __name__ == '__main__':
    import nose
    nose.runmodule()
    
