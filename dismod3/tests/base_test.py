from django.test import TestCase
from django.test.client import Client

class DisModTestCase(TestCase):
    fixtures = ['dismod3/fixtures/populations',
                'dismod3/fixtures/regions',
                'dismod3/fixtures/diseases',
                'dismod3/fixtures/rates',
                'dismod3/fixtures/age_specific_rate_functions',
                'dismod3/fixtures/disease_models'
                ]

    def assertPng(self, response):
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content[1:4], 'PNG')  # is there a better way to test that the response is a png?
