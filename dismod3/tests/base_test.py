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

    def create_users(self):
        """
        it seems easier to create the users with a code block than as
        json fixtures, because the password is clearer before it is
        encrypted
        """
        from django.contrib.auth.models import User
        user = User.objects.create_user('red', '', 'red')
        user = User.objects.create_user('green', '', 'green')
        user = User.objects.create_user('blue', '', 'blue')
        

    def assertPng(self, response):
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.content[1:4], 'PNG')  # is there a better way to test that the response is a png?
