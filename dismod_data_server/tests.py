from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
from models import *

class DisModDataServerTestCase(TestCase):
    fixtures = ['dismod_data_server/fixtures',
                'population_data_server/fixtures']

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
        self.assertTrue(age_weights[0] > age_weights[-1])

    def test_create_disease_model(self):
        """ Test creating a dismod model object from a dismod_data json string"""

        json_str = self.dm.to_json()
        dm2 = create_disease_model(json_str)
        self.assertTrue(dm2.id != self.dm.id and
                        dm2.id == DiseaseModel.objects.latest('id').id)
        
    # functional tests
    #### Data Loading requirements

    def test_dismod_load_well_formed_csv(self):
        """ Make sure that a properly formatted data csv can be loaded over the web"""
        c = Client()

        # first check that create requires a login
        url = reverse('gbd.dismod_data_server.views.data_upload')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_upload.html')
        
        response = c.post(url, {'tab_separated_values': ''})
        self.assertTemplateUsed(response, 'data_upload.html')

        # now do it right, and make sure that data and datasets are added
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalence\tTotal\tCanada\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})

        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        self.assertEqual([1.]*10, Data.objects.latest('id').params.get('age_weights'))

    def test_dismod_informative_error_for_badly_formed_csv(self):
        """ Provide informative error if data csv cannot be loaded"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        # csv with required column, GBD Cause,  missing 
        response = c.post(url, {'tab_separated_values': \
        'Region\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nWorld\tPrevalence\tTotal\tAustralia\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})
        self.assertContains(response, 'GBD Cause')
        self.assertContains(response, 'is missing')

        # csv with cell missing from line 2
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalence\tTotal\tAustralia\t15\t24\t2005\t2005\t.5\t.1\tper 1.0'})
        self.assertContains(response, 'Error loading row 2:')

        # csv with unrecognized parameter
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalenceee\tTotal\tAustralia\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})
        self.assertContains(response, 'Row 2:  could not understand entry for Parameter')

    def test_dismod_add_age_weights(self):
        """ Use the Population Data Server to get the age weights for a new piece of data"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalence\tTotal\tAustralia\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})

        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        age_weights = Data.objects.latest('id').params.get('age_weights')
        # the fixture for Australia 2005 total population has a downward trend
        assert age_weights[0] > age_weights[-1]

    def test_dismod_add_covariates(self):
        """ Use the Covariate Data Server to get the covariates for a new piece of data"""
        c = Client()
        url = reverse('gbd.dismod_data_server.views.data_upload')
        c.login(username='red', password='red')

        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalence\tTotal\tAustralia\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})

        assert Data.objects.latest('id').params.has_key('gdp'), \
            'should add GDP data from covariate data server (not yet implemented)'
        
        
        
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

    def test_dismod_show_in_other_formats(self):
        """ Test displaying disease model as png, json, csv, etc"""
        c = Client()
        c.login(username='red', password='red')
        url = self.dm.get_absolute_url()

        # test png
        response = c.get(url + '.png')
        self.assertPng(response)


    
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
        self.assertEqual(set(r_json.keys()), set(['params', 'data']))
        
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

        # now check that good input is accepted
        initial_dm_cnt = DiseaseModel.objects.count()
        response = c.post(url, {'model_json': self.dm.to_json()})
        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        self.assertEqual(DiseaseModel.objects.count(), initial_dm_cnt+1)

