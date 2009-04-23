from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
from models import *

class DisModTestCase(TestCase):
    fixtures = ['new_dm3/fixtures']

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


class DataTestCase(DisModTestCase):
    def setUp(self):
        from models import Data
        self.data = Data.objects.all()[0]
        self.create_users()

    # unit tests
    def test_str(self):
        s = str(self.data)
        self.assertTrue(isinstance(s,str))

        s = self.data.get_absolute_url()
        self.assertTrue(isinstance(s,str))

    # functional tests
    def test_data_show(self):
        c = Client()

        # first check that show requires login
        response = c.get(self.data.get_absolute_url())
        self.assertRedirects(response, 'http://testserver/accounts/login/?next=/new/data/1')

        # then check that show works after login
        c.login(username='red', password='red')
        response = c.get(self.data.get_absolute_url())
        self.assertTemplateUsed(response, 'data_show.html')

        
    def test_data_create(self):
        c = Client()

        # first check that create requires a login

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(reverse("new_dm3.views.data_new"))
        self.assertTemplateUsed(response, 'data_new.html')
        
        response = c.post(reverse("new_dm3.views.data_new"), {'tab_separated_values': ''})
        self.assertTemplateUsed(response, 'data_new.html')

        # now do it right, and make sure that data and datasets are added
        response = c.post(reverse("new_dm3.views.data_new"), {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tEstimate Year Start\tEstimate Year End\tParameter Value\tLower Value\tUpper Value\tUnits\tType of Bounds\nCannabis Dependence\tWorld\tPrevalence\tTotal\tCanada\t15\t24\t2005\t2005\t.5\t.4\t.6\t1.0\t95% CI'})

        self.assertRedirects(response, Data.objects.all()[Data.objects.count()-1].get_edit_url())
