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

    def test_interpolate(self):
        """ Test interpolation"""
        f = self.pop.interpolate([0,10,20])
        self.assertEqual(f[0], 1.5)

    # functional tests
    def test_population_show(self):
        """ Test plotting population curve"""
        c = Client()

        url = self.pop.get_absolute_url()
        response = c.get(url)
        self.assertPng(response)

    def test_population_show_in_other_formats(self):
        """ Test getting population curve as json, csv, etc"""
        c = Client()

        # test png
        url = self.pop.get_absolute_url()
        response = c.get(url + '.png')
        self.assertPng(response)

        # test json
        response = c.get(url + '.json')
        r_json = json.loads(response.content)
        self.assertEqual(set(r_json.keys()), set(['age', 'population']))

        # test csv
        response = c.get(url + '.csv')
        self.assertEqual(response.content.split('\r\n')[0], 'Age (years),Population (thousands)')
