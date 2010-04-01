from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
import urllib

from models import *

class DisModDataServerTestCase(TestCase):
    fixtures = ['dismod_data_server/fixtures',
                'population_data_server/fixtures',
                'covariate_data_server/fixtures']

    def create_users(self):
        """ Create users for functional testing of access control.

        It seems easier to create the users with a code block than as
        json fixtures, because the password is clearer before it is
        encrypted.
        """
        from django.contrib.auth.models import User
        user = User.objects.create_user('red', '', 'red')
        user = User.objects.create_user('green', '', 'green')
        user = User.objects.create_user('blue', '', 'blue')

    def assertPng(self, response):
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content[1:4], 'PNG')

    def assertSuccess(self, response):
        return self.assertEquals(response.status_code, 200)

    def assertNotFound(self, response):
        return self.assertEquals(response.status_code, 404)

    def setUp(self):
        self.dm = DiseaseModel.objects.latest('id')
        self.data = Data.objects.latest('id')

        self.create_users()

    # unit tests
    def test_str(self):
        """ Test all model string functions"""
        s = str(self.data)
        self.assertTrue(isinstance(s,str))

        s = self.data.get_absolute_url()
        self.assertTrue(isinstance(s,str))

        s = str(self.dm)
        self.assertTrue(isinstance(s,str))

        s = self.dm.get_absolute_url()
        self.assertTrue(isinstance(s,str))

    def test_calculate_and_cache_age_weights(self):
        """ Test that dismod data object can query to population data server"""
        self.assertFalse(self.data.params.has_key('age_weights'))

        age_weights = self.data.age_weights()
        self.assertTrue(self.data.params.has_key('age_weights'))

        # fixture has population skewed towards youth
        self.assertTrue(age_weights[0] > age_weights[1])

    def test_create_disease_model(self):
        """ Test creating a dismod model object from a dismod_data json string"""

        json_str = self.dm.to_json()
        dm2 = create_disease_model(json_str, User.objects.latest('id'))
        self.assertTrue(dm2.id != self.dm.id and
                        dm2.id == DiseaseModel.objects.latest('id').id)
        
    # functional tests
    #### Data Loading requirements

    def test_dismod_load_data_from_file(self):
        """ Make sure that a properly formatted data csv file can be loaded over the web"""

        c = Client()

        # first check that create requires a login
        url = reverse('gbd.dismod_data_server.views.data_upload')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_upload.html')

        response = c.post(url, {})
        self.assertTemplateUsed(response, 'data_upload.html')

        # now do it right, and make sure that data and datasets are added
        f = open("tests/data.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[DiseaseModel.objects.latest('id').id]))


    def test_dismod_load_data_ci(self):
        """ Make sure that data csv with confidence intervals is interpreted properly"""

        c = Client()
        c.login(username='red', password='red')
        url = reverse('gbd.dismod_data_server.views.data_upload')
        response = c.post(url, {'file': open("tests/data_ci.tsv")})


        # check that bounds on data from latest model are right
        d = Data.objects.latest('id')
        self.assertEqual(d.value, .25)
        self.assertEqual(float(d.params['lower_ci']), .2)
        self.assertEqual(float(d.params['upper_ci']), .3)

        from dismod3.disease_json import DiseaseJson
        dm = DiseaseJson(DiseaseModel.objects.latest('id').to_json())
        dm.calc_effective_sample_size(dm.data)
        d = dm.data[-1]
        
        self.assertEqual(dm.bounds_per_1(d), (.2, .3))

        se = (.3-.2)/(2*1.96)
        self.assertEqual(dm.se_per_1(d), se)

        n = .25*(1-.25)/se**2
        self.assertEqual(d['effective_sample_size'], n)

    def test_dismod_informative_error_for_badly_formed_data_file(self):
        """ Provide informative error if data file cannot be loaded"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        # data with required column, GBD Cause,  missing 
        f = open("tests/data_column_missing.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'GBD Cause')
        self.assertContains(response, 'is missing')

        # data with cell missing from line 2
        f = open("tests/data_cell_missing.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Error loading row 2:')

        # data with inconsistent gbd cause from line 2
        f = open("tests/data_inconsistent_gbd_cause.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 4:  could not understand entry for GBD Cause inconsistent')

        # data with wrong region from line 2
        f = open("tests/data_region.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Region')

        # data with wrong sex from line 2
        f = open("tests/data_sex.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Sex')

        # data with wrong sex from line 2
        f = open("tests/data_country_iso3_code.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Country ISO3 Code')

        # data with wrong ages from line 2
        f = open("tests/data_age.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Age Start greater than Age End')

        # data with wrong ages from line 2
        f = open("tests/data_year.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Year Start greater than Year End')

        # data with unrecognized parameter
        f = open("tests/data_unrecognized_parameter.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Parameter')

        # data with wrong parameter value
        f = open("tests/data_parameter_value.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Parameter Value less than 0')

        # data with wrong units
        f = open("tests/data_units.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Units less than 1')

        # data with wrong staudy id
        f = open("tests/data_study_id.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Study ID less than 0')

        # data with wrong coverage
        f = open("tests/data_coverage.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Coverage out of range [0, 1]')

        # data with wrong study size n for this year & sex
        f = open("tests/data_study_size_n_year_sex.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Study Size N For This Year and Sex less than or equal 0')

        # data with wrong ci from line 4
        f = open("tests/data_wrong_ci.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 4:  could not understand entry for Upper CI less than Parameter Value')

        # data with wrong standard error from line 2
        f = open("tests/data_standard_error.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Standard Error less than or equal 0')

        # data with wrong total_study_size_n from line 2
        f = open("tests/data_total_study_size_n.tsv")
        response = c.post(url, {'file':f})
        f.close()
        self.assertContains(response, 'Row 2:  could not understand entry for Total Study Size N less than or equal 0')

    def test_dismod_add_age_weights_to_data_file(self):
        """ Use the Population Data Server to get the age weights for a new piece of data"""

        # TODO: fix this test, which was broken when asynchronous
        # age_weight calculation was added to views.py

        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        f = open("tests/data_add_age_weights.tsv")
        response = c.post(url, {'file':f})
        f.close()

        id = DiseaseModel.objects.latest('id').id
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[id]))

        response = c.post(reverse('gbd.dismod_data_server.views.dismod_update_covariates', args=[id]))
        age_weights = Data.objects.latest('id').params.get('age_weights')
        # the fixture for Australia 2005 total population has a downward trend
        assert age_weights[0] > age_weights[1]
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_run', args=[DiseaseModel.objects.latest('id').id]))

    def test_dismod_set_covariates(self):
        """ Load the covariate selection panel for a new piece of data"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.dismod_set_covariates', args=[self.dm.id])
        c.login(username='red', password='red')

        self.dm.data.add(self.data)

        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_set_covariates.html')

    def test_dismod_add_additional_data_to_model_file(self):
        """ Test adding data from csv to existing model"""
        c = Client()
        c.login(username='red', password='red')

        self.data.cache_params()
        self.data.save()
        
        self.dm.data.add(self.data)
        
        url = reverse('gbd.dismod_data_server.views.data_upload', args=(self.dm.id,))

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_upload.html')

        f = open("tests/data.tsv")
        response = c.post(url, {'file':f})
        f.close()

        newest_data = Data.objects.latest('id')
        newest_dm = DiseaseModel.objects.latest('id')
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[newest_dm.id]))
        self.assertEqual(sorted([d.id for d in self.dm.data.all()] + [newest_data.id]),
                         sorted([d.id for d in newest_dm.data.all()]))

    def test_dismod_load_well_formed_data_csv(self):
        """ Make sure that a properly formatted data csv can be loaded over the web"""

        # TODO: fix this test, which was broken when asynchronous
        # age_weight calculation was added to views.py
        
        c = Client()

        # first check that create requires a login
        url = reverse('gbd.dismod_data_server.views.data_upload')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_upload.html')
        
        response = c.post(url, {'tab_separated_values': '', 'file': ''})
        self.assertTemplateUsed(response, 'data_upload.html')

        # now do it right, and make sure that data and datasets are added
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry ISo3 Code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tNorth America, High Income\tPrevalence\tTotal\tCAN\t15\t24\t2005\t2005\t.5\tper 1.0'})

        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[DiseaseModel.objects.latest('id').id]))
        #self.assertEqual([1.]*10, Data.objects.latest('id').params.get('age_weights'))

    def test_dismod_informative_error_for_badly_formed_data_csv(self):
        """ Provide informative error if data csv cannot be loaded"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        # csv with required column, GBD Cause,  missing 
        response = c.post(url, {'tab_separated_values': \
        'Region\tParameter\tSex\tCountry iso3 code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nAustralasia\tPrevalence\tTotal\tAUS\t15\t24\t2005\t2005\t.5\tper 1.0'})
        self.assertContains(response, 'GBD Cause')
        self.assertContains(response, 'is missing')

        # csv with cell missing from line 2
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry ISo3 code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tAustralasia\tPrevalence\tTotal\tAUS\t15\t24\t2005\t2005\t.5'})
        self.assertContains(response, 'Error loading row 2:')

        # csv with unrecognized parameter
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry iSo3 code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tAustralasia\tPrevalenceee\tTotal\tAUS\t15\t24\t2005\t2005\t.5\tper 1.0'})
        self.assertContains(response, 'Row 2:  could not understand entry for Parameter')

    def test_dismod_add_age_weights_to_data(self):
        """ Use the Population Data Server to get the age weights for a new piece of data"""

        # TODO: fix this test, which was broken when asynchronous
        # age_weight calculation was added to views.py

        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry iSO3 code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tAustralasia\tPrevalence\tTotal\tAUS\t15\t24\t2005\t2005\t.5\tper 1.0'})

        id = DiseaseModel.objects.latest('id').id
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[id]))

        response = c.post(reverse('gbd.dismod_data_server.views.dismod_update_covariates', args=[id]))
        age_weights = Data.objects.latest('id').params.get('age_weights')
        # the fixture for Australia 2005 total population has a downward trend
        assert age_weights[0] > age_weights[1]
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_run', args=[DiseaseModel.objects.latest('id').id]))

    def test_dismod_add_covariates_to_data(self):
        """ Use the Covariate Data Server to get the covariates for a new piece of data"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry iso3_code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tNorth America, High Income\tPrevalence\tTotal\tUSA\t15\t24\t2005\t2005\t.5\tper 1.0'})

        dm = DiseaseModel.objects.latest('id')
        dm.params.create(key='covariates', json=json.dumps({'Country_level':{'GDP': {'rate': {'value':1}}}}))
        url = reverse('gbd.dismod_data_server.views.dismod_update_covariates',
                      args=[dm.id])
        response = c.post(url)

        assert dm.data.latest('id').params.has_key('gdp'), \
            'should add GDP data from covariate data server'
        
    def test_dismod_add_additional_data_to_model(self):
        """ Test adding data from csv to existing model"""
        c = Client()
        c.login(username='red', password='red')

        self.data.cache_params()
        self.data.save()
        
        self.dm.data.add(self.data)
        
        url = reverse('gbd.dismod_data_server.views.data_upload', args=(self.dm.id,))

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_upload.html')

        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry iso3 code\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tUnits\nCannabis Dependence\tNorth America, High Income\tPrevalence\tTotal\tCAN\t15\t24\t2010\t2010\t.5\tper 1.0'})

        newest_data = Data.objects.latest('id')
        newest_dm = DiseaseModel.objects.latest('id')
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_summary', args=[newest_dm.id]))
        self.assertEqual(sorted([d.id for d in self.dm.data.all()] + [newest_data.id]),
                         sorted([d.id for d in newest_dm.data.all()]))
        
        
    #### Model Viewing requirements
    def test_data_show(self):
        """ Test displaying html version of a single data point"""
        c = Client()

        # first check that show requires login
        url = self.data.get_absolute_url()
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that show works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'data_show.html')

    def test_dismod_list(self):
        """ Test listing the existing disease models"""
        c = Client()

        # first check that show requires login
        url = reverse('gbd.dismod_data_server.views.dismod_list')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_list.html')

    def test_dismod_show(self):
        """ Test displaying html version of a disease model"""
        c = Client()

        # first check that show requires login
        url = self.dm.get_absolute_url()
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that show works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_show.html')

    def test_dismod_error_for_wrong_region_show(self):
        """ Test displaying non-existing region"""
        c = Client()
        c.login(username='red', password='red')

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region', args=[self.dm.id, 'theoryland'])
        response = c.get(url)
        self.assertNotFound(response)

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region', args=[self.dm.id, 'asia_east'])
        response = c.get(url)
        self.assertSuccess(response)

    def test_dismod_error_for_wrong_region_year_sex_show(self):
        """ Test non-existing region, year, or sex show"""
        c = Client()

        c.login(username='red', password='red')

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region_year_sex', args=[self.dm.id, 'asia_east', 2005, 'male'])
        response = c.get(url)
        self.assertSuccess(response)

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region_year_sex', args=[self.dm.id, 'theoryland', 2005, 'male'])
        response = c.get(url)
        self.assertNotFound(response)

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region_year_sex', args=[self.dm.id, 'asia_east', 1999, 'male'])
        response = c.get(url)
        self.assertNotFound(response)

        url = reverse('gbd.dismod_data_server.views.dismod_show_by_region_year_sex', args=[self.dm.id, 'asia_east', 2005, 'unspecified'])
        response = c.get(url)
        self.assertNotFound(response)

    def test_dismod_show_emp_priors(self):
        """ Test displaying empirical priors"""
        c = Client()
        c.login(username='red', password='red')
        url = reverse('gbd.dismod_data_server.views.dismod_show_emp_priors', args=[self.dm.id])
        response = c.get(url)
        self.assertSuccess(response)

        url = reverse('gbd.dismod_data_server.views.dismod_show_emp_priors', args=[self.dm.id, 'png'])
        response = c.get(url)
        self.assertSuccess(response)
        self.assertPng(response)

    def test_dismod_show_in_other_formats(self):
        """ Test displaying disease model as png, json, csv, etc"""
        c = Client()
        c.login(username='red', password='red')
        url = self.dm.get_absolute_url()

        # test png
        # (commented out, because it takes 30 seconds to run!)
        #response = c.get(url + '.png')
        #self.assertPng(response)

    def test_dismod_show_map(self):
        """ Test displaying map of disease model"""
        c = Client()
        c.login(username='red', password='red')
        url = reverse('gbd.dismod_data_server.views.dismod_show_map', args=[self.dm.id])

        response = c.post(url, {'year': '1990', 'sex': 'male', 'type': 'prevalence', 'map': 'data', 'data_count': 'Data Count Map', 'moment': 'median', 'age': 'all ages', 'weight': 'direct'})
        self.assertTemplateUsed(response, 'dismod_map.svg')

        response = c.post(url, {'year': '1997', 'sex': 'female', 'type': 'incidence', 'map': 'data_count', 'moment': 'median', 'age': 'all ages', 'weight': 'direct'})
        self.assertTemplateUsed(response, 'dismod_map.svg')

        #response = c.post(url, {'year': '2005', 'sex': 'total', 'type': 'remission', 'map': 'emp-prior', 'moment': 'sum', 'age': 'all ages', 'weight': 'weighted'})
        #self.assertTemplateUsed(response, 'dismod_map.svg')

    def test_dismod_compare_emp_priors(self):
        """ Test displaying map of disease model"""
        c = Client()
        c.login(username='red', password='red')
        url = reverse('gbd.dismod_data_server.views.dismod_comparison')
        response = c.get(url)

        url = reverse('gbd.dismod_data_server.views.dismod_compare', args=[self.dm.id, self.dm.id, 'alpha', 'png'])

        response = c.get(url)
        self.assertPng(response)

        url = reverse('gbd.dismod_data_server.views.dismod_compare', args=[self.dm.id, self.dm.id, 'overlay+prevalence+asia_southeast+1990+male', 'png'])

        response = c.get(url)
        self.assertPng(response)

    def test_dismod_sparkplot(self):
        """ Test sparkplot of disease model"""
        c = Client()

        # first check that sparkplot requires login
        url = '/dismod/show/spark_%d.png' % self.dm.id
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertPng(response)

    def test_dismod_overlay_plot(self):
        """ Test overlay plot of disease model"""
        c = Client()

        # first check that overlay plot requires login
        url = '/dismod/show/overlay_1_CHD+all+latin_america_southern+1995+male.png'
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s' % urllib.quote(url))

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertPng(response)

    def test_dismod_tile_plot(self):
        """ Test tile plot of disease model"""
        c = Client()

        # first check that overlay plot requires login
        url = '/dismod/show/tile_1_CHD+all+latin_america_southern+1995+male.png'
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s' % urllib.quote(url))

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertPng(response)

    def test_dismod_table_each_age(self):
        """ Test table of disease model"""
        c = Client()

        # first check that making tables requires login
        url = '/dismod/show/1.xls'
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s' % urllib.quote(url))

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertSuccess(response)

    def test_dismod_table_group_10(self):
        """ Test table of disease model"""
        c = Client()
        c.login(username='red', password='red')

        url = '/dismod/show/1.xls'
        response = c.get(url, dict(group_size=10))
        self.assertSuccess(response)

    def test_dismod_summary(self):
        """ Test the model summary view"""
        c = Client()

        # first check that overlay plot requires login
        url = reverse('gbd.dismod_data_server.views.dismod_summary', args=[self.dm.id])
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_summary.html')

    def test_dismod_export(self):
        """ Test the model export view"""
        c = Client()

        # first check that overlay plot requires login
        url = reverse('gbd.dismod_data_server.views.dismod_export', args=[self.dm.id])
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that it works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_export.html')
        
    
    #### Model Running requirements
    def test_get_model_json(self):
        """ Test getting a json encoding of the disease model"""
        c = Client()
        
        # first check that getting json requires login
        url = self.dm.get_absolute_url() + '.json'
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # now login, and check that you can get json
        c.login(username='red', password='red')
        response = c.get(url)
        r_json = json.loads(response.content)
        self.assertEqual(set(r_json.keys()), set(['params', 'data', 'id']))
        
    def test_post_model_json(self):
        """ Test posting a json encoding of the disease model"""
        c = Client()

        # first check that create requires a login
        url = reverse('gbd.dismod_data_server.views.dismod_upload')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_upload.html')

        # check that bad input is rejected
        # TODO: make this part of the test more extensive
        response = c.post(url, {'model_json': ''})
        self.assertTemplateUsed(response, 'dismod_upload.html')

        # check that if input is good and params.id equals a valid
        # model id, that model is updated
        p = DiseaseModelParameter(key='map', json=json.dumps({'prevalence': [0,0,0,0], 'incidence': [0,0,0,0]}))
        p.save()
        self.dm.params.add(p)
        response = c.post(url, {'model_json': self.dm.to_json()})
        self.assertRedirects(response, self.dm.get_absolute_url())
        dm = DiseaseModel.objects.get(id=self.dm.id)
        self.assertEqual(json.loads(dm.params.filter(key='map').latest('id').json)['prevalence'], [0,0,0,0])
        

        # now check that good input is accepted, and if params.id = -1
        # a new model is created
        initial_dm_cnt = DiseaseModel.objects.count()
        self.dm.id = -1
        response = c.post(url, {'model_json': self.dm.to_json()})
        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        self.assertEqual(DiseaseModel.objects.count(), initial_dm_cnt+1)
        

    def test_get_job_queue_list_and_remove(self):
        """ Test getting list of jobs waiting on queue to run"""
        c = Client()
        c.login(username='red', password='red')

        # make a model need to run
        p = DiseaseModelParameter(key='needs_to_run')
        p.save()
        self.dm.params.add(p)
        self.dm.save()

        # test GET list
        url = reverse('gbd.dismod_data_server.views.job_queue_list')
        response = c.get(url, {'format': 'json'})

        r_json = json.loads(response.content)
        self.assertEqual(r_json, [self.dm.params.filter(key='needs_to_run').latest('id').id])

        # test GET&POST remove
        self.assertTrue(self.dm.params.filter(key='needs_to_run').count() == 1)
        url = reverse('gbd.dismod_data_server.views.job_queue_remove')
        response = c.get(url)
        response = c.post(url, {'id': p.id})

        dm = DiseaseModel.objects.get(id=self.dm.id)
        self.assertTrue(self.dm.params.filter(key='needs_to_run').count() == 0)
        

    def test_job_queue_add(self):
        """ Test adding a job to job queue to run"""
        c = Client()
        c.login(username='red', password='red')

        self.assertTrue(self.dm.params.filter(key='needs_to_run').count() == 0)
        url = reverse('gbd.dismod_data_server.views.job_queue_add', args=[self.dm.id])
        response = c.post(url, {'estimate_type': 'fit each region/year/sex individually'})
        dm = DiseaseModel.objects.get(id=self.dm.id)
        self.assertTrue(dm.params.filter(key='needs_to_run').count()==1)
        p = dm.params.filter(key='needs_to_run').latest('id')
        self.assertEqual(json.loads(p.json).get('estimate_type'), 'fit each region/year/sex individually')

    def test_dismod_run(self):
        """ Test adding a job to job queue to run"""
        c = Client()
        c.login(username='red', password='red')

        url = reverse('gbd.dismod_data_server.views.dismod_run', args=[self.dm.id])
        response = c.get(url)
        self.assertTemplateUsed(response, 'dismod_run.html')

    def test_dismod_update_covariates(self):
        """ Test updating age weights for a model"""
        c = Client()
        c.login(username='red', password='red')

        url = reverse('gbd.dismod_data_server.views.dismod_update_covariates', args=[self.dm.id])
        response = c.post(url)
        self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_run', args=[self.dm.id]))
        

    # Model Adjusting Requirements
# This is currently handled by a Java Applet, but the test is for the previous implementation,
# an html form.  Perhaps we will go back to an html form in the future.
#     def test_dismod_adjust(self):
#         """ Test changing priors and ymax for a model"""
#         c = Client()

#         # first check that this requires login
#         url = reverse('gbd.dismod_data_server.views.dismod_adjust', args=[self.dm.id])
#         response = c.get(url)
#         self.assertRedirects(response, '/accounts/login/?next=%s'%url)

#         # now login, and check that you can get to adjust page
#         c.login(username='red', password='red')
#         response = c.get(url)

#         # post with no adjustments, this should not change the model
#         response = c.post(url, {})
#         self.assertSuccess(response)
        
#         # now make ymax adjustments, and check that they actually change ymax
#         response = c.post(url, {'ymax' : '0.1'})
        
#         dm = DiseaseModel.objects.get(id=self.dm.id)
#         self.assertRedirects(response, dm.get_absolute_url())
#         self.assertEqual(dm.params['ymax'], .1)
        
#         # now make prior adjustments, and check
#         response = c.post(url, {'prevalence_smoothness' : '10.0'})
        
#         new_dm = DiseaseModel.objects.latest('id')
#         self.assertRedirects(response, reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))
#         self.assertEqual(new_dm.params['priors']['prevalence+north_america_high_income+2005+male'], 'smooth 10.0, ')

    def test_dismod_preview_priors(self):
        """ Test generating png to preview priors"""
        c = Client()
        c.login(username='red', password='red')
        url = reverse('gbd.dismod_data_server.views.dismod_preview_priors', args=[self.dm.id])

        # test get
        response = c.get(url)
        self.assertSuccess(response)
        self.assertPng(response)

        # test post
        response = c.post(url, {'JSON': json.dumps({'smoothing': {'incidence': 'pretty smooth', 'prevalence': 'hello, world'},
                                                    'parameter_age_mesh': [0, 10, 100],
                                                    'y_maximum': 1.0,
                                                    'note': '',
                                                    
                                                    })})
        self.assertSuccess(response)
        self.assertPng(response)
