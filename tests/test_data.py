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
        """        
    def test_from_gbd_json(self):
        d = data.ModelData.from_gbd_json('tests/dismoditis.json')

        assert len(d.input_data) > 17, 'dismoditis model has more than 17 data points'
        for field in 'data_type value area sex age_start age_end year_start year_end standard_error effective_sample_size lower_ci upper_ci age_weights'.split():
            assert field in d.input_data.columns, 'Input data CSV should have field "%s"' % field
        assert len(d.input_data.filter(regex='x_').columns) == 1, 'should have added country-level covariates to input data'
        assert len(d.input_data['x_LDI_id_Updated_7July2011'].dropna().index) > 0

        assert len(d.output_template) > 100
        for field in 'data_type area sex age_start age_end year_start year_end age_weights'.split():
            assert field in d.output_template.columns, 'Output template CSV should have field "%s"' % field
        assert len(d.output_template.filter(regex='x_').columns) == 1, 'should have added country-level covariates to output template'
        assert len(d.output_template['x_LDI_id_Updated_7July2011'].dropna().index) > 0
        assert len(eval(d.output_template['age_weights'][0])) == d.output_template['age_end'][0] - d.output_template['age_start'][0]
        assert eval(d.output_template['age_weights'][1])[0] != eval(d.output_template['age_weights'][1])[1], 'output template population age weights should be changing'

        for data_type in 'i p r f rr X'.split():
            for prior in 'smoothness heterogeneity level_value level_bounds increasing decreasing'.split():
                assert prior in d.parameters[data_type], 'Parameters for %s should include prior on %s' % (data_type, prior)

        assert d.areas_hierarchy.successors('asia_east') == ['MAC', 'PRK', 'TWN', 'HKG', 'CHN']
        assert len(d.areas_to_fit) == 22
        """
    def test_save_and_load(self):
        d = data.ModelData.from_gbd_json('tests/dismoditis.json')

        # TODO: delete this dir if it exists
        d.save('/var/tmp/test_dm3/')

        # TODO: test that files really were created

        d2 = data.ModelData.load('/var/tmp/test_dm3/')
        assert d.input_data.shape == d2.input_data.shape, 'input data should be equal before and after save'
        assert pl.all(d.input_data['value'] == d2.input_data['value']), 'input data should be equal before and after save'

        assert d.output_template.shape == d2.output_template.shape, 'output template should be equal before and after save'
        assert pl.all(d.output_template['area'] == d2.output_template['area']), 'output template should be equal before and after save'

        assert d.parameters == d2.parameters, 'parameters should be equal before and after save'
        assert sorted(d.areas_hierarchy.edges()) == sorted(d2.areas_hierarchy.edges()), 'areas_hierarchy should be equal before and after save'
        assert d.areas_to_fit == d2.areas_to_fit, 'areas_to_fit should be equal before and after save'
if __name__ == '__main__':
    import nose
    nose.runmodule()
    
