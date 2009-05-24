from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
from models import *

class PopulationDataServerTestCase(TestCase):
    fixtures = ['population_data_server/fixtures']

    def assertPng(self, response):
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content[1:4], 'PNG')

    def assertSuccess(self, response):
        return self.assertEquals(response.status_code, 200)

    def setUp(self):
        self.pop = Population.objects.all()[0]

    # unit tests
    def test_str(self):
        """ Test all model string functions"""
        s = str(self.pop)
        self.assertTrue(isinstance(s,str))

        s = self.pop.get_absolute_url()
        self.assertTrue(isinstance(s,str))

    def test_gp(self):
        """ Test Gaussian Process interpolation"""
        M, C = self.pop.gaussian_process()
        self.assertEqual(M(0), 1.)

    # functional tests
    def test_population_show(self):
        """ Test plotting population curve"""
        c = Client()

        # first check that show requires login
        url = self.pop.get_absolute_url()
        response = c.get(url)
        self.assertPng(response)
        
