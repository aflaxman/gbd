from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
from models import *

class CovariateDataServerTestCase(TestCase):
    fixtures = ['covariate_data_server/fixtures']

    def assertPng(self, response):
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content[1:4], 'PNG')

    def assertSuccess(self, response):
        return self.assertEquals(response.status_code, 200)

    def setUp(self):
        self.cov = Covariate.objects.all()[0]

    # unit tests
    def test_str(self):
        """ Test all model string functions"""
        s = str(self.cov)
        self.assertTrue(isinstance(s,str))

        s = self.cov.get_absolute_url()
        self.assertTrue(isinstance(s,str))

