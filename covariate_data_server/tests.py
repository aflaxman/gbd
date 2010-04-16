from django.test import TestCase
from django.test.client import Client

from django.core.urlresolvers import reverse
from models import *

class CovariateDataServerTestCase(TestCase):
    fixtures = ['covariate_data_server/fixtures']

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
        self.ctype = CovariateType.objects.all()[0]
        self.cov = Covariate.objects.all()[0]

        self.create_users()

    # unit tests
    def test_str(self):
        """ Test all model string functions"""
        s = str(self.cov)
        self.assertTrue(isinstance(s,str))

        s = self.cov.get_absolute_url()
        self.assertTrue(isinstance(s,str))

    # functional tests
    def test_show(self):
        c = Client()

        url = self.cov.get_absolute_url()
        response = c.get(url)
        self.assertPng(response)

        url = self.ctype.get_absolute_url()
        response = c.get(url)
        
    def test_upload(self):
        c = Client()

        url = reverse('gbd.covariate_data_server.views.covariate_upload')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)
        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'covariate_upload.html')

        response = c.post(url, {})
        self.assertTemplateUsed(response, 'covariate_upload.html')

        # now do it right, and make sure that data and datasets are added
        from StringIO import StringIO
        f = StringIO(',iso3,year,LDI_id,LDI_usd\n1,ABW,1950,1533.743774,1105.747437\n1,ABW,1951,1533.843774,1105.87437\n')
        f.name = 'LDI.csv'
        response = c.post(url, {'file':f, 'type': 'LDI_id'})
        self.assertRedirects(response, reverse('gbd.covariate_data_server.views.covariate_type_show', args=[CovariateType.objects.latest('id').id]))
        
        self.assertEqual(CovariateType.objects.filter(slug='LDI_id').count(), 1)
        self.assertEqual(Covariate.objects.filter(type__slug='LDI_id', sex='male').count(), 2)
        self.assertEqual(Covariate.objects.filter(type__slug='LDI_id', sex='female').count(), 2)
        self.assertEqual(Covariate.objects.filter(type__slug='LDI_id', sex='total').count(), 2)
