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

    def assertSuccess(self, response):
        return self.assertEquals(response.status_code, 200)



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
        url = self.data.get_absolute_url()
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then check that show works after login
        c.login(username='red', password='red')
        response = c.get(url)
        self.assertTemplateUsed(response, 'data_show.html')

    def test_data_create(self):
        c = Client()

        # first check that create requires a login
        url = reverse('new_dm3.views.data_new')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'data_new.html')
        
        response = c.post(url, {'tab_separated_values': ''})
        self.assertTemplateUsed(response, 'data_new.html')

        # now do it right, and make sure that data and datasets are added
        response = c.post(url, {'tab_separated_values': \
        'GBD Cause\tRegion\tParameter\tSex\tCountry\tAge Start\tAge End\tYear Start\tYear End\tParameter Value\tStandard Error\tUnits\tType of Bound\nCannabis Dependence\tWorld\tPrevalence\tTotal\tCanada\t15\t24\t2005\t2005\t.5\t.1\tper 1.0\t95% CI'})

        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        self.assertEqual([1.]*10, Data.objects.latest('id').params.get('age_weights'))
        



class DiseaseModelTestCase(DisModTestCase):
    def setUp(self):
        from models import DiseaseModel
        self.dm = DiseaseModel.objects.all()[0]
        self.create_users()

    # unit tests
    def test_str(self):
        s = str(self.dm)
        self.assertTrue(isinstance(s,str))

        s = self.dm.get_absolute_url()
        self.assertTrue(isinstance(s,str))

#     def test_fit(self):
#         import dismod3
#         dismod3.fit(1)
    
    # functional tests
    def test_disease_model_show(self):
        c = Client()

        # first check that create requires a login
        url = self.dm.get_absolute_url()
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # then login and do functional tests
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'disease_model_show.html')
        
        # try getting a json version of this disease model, as well
        response = c.get(reverse("new_dm3.views.disease_model_show", args=(self.dm.id,'json')))
        # FIXME: this brittle test is now wrong
        #self.assertEqual(response.content.replace('\n','').replace(' ',''),
        #                 ('{"data": [], "params": %s}' % json.dumps(self.dm.params)).replace('\n','').replace(' ',''))

        # and finally, try getting a png version
        d = Data.objects.all()[0]
        d.cache_params()
        d.save()
        response = c.get(reverse("new_dm3.views.disease_model_show", args=(self.dm.id,'png')))
        self.assertPng(response)

        
    def test_disease_model_new(self):
        c = Client()

        # first check that create requires a login
        url = reverse('new_dm3.views.disease_model_new')
        response = c.get(url)
        self.assertRedirects(response, '/accounts/login/?next=%s'%url)

        # login and do it again
        c.login(username='red', password='red')

        response = c.get(url)
        self.assertTemplateUsed(response, 'disease_model_new.html')
        
        response = c.post(url, {'model_json': ''})
        self.assertTemplateUsed(response, 'disease_model_new.html')

        # now do it right, and make sure that data and datasets are added
        response = c.post(url, {'model_json': json.dumps({'params': {'condition': 'test disease', 'sex': 'male', 'region': 'World', 'year': 1999}, 'data': [{'id': 1}]})})

        self.assertRedirects(response, DiseaseModel.objects.latest('id').get_absolute_url())
        self.assertEqual([1], [d.id for d in DiseaseModel.objects.latest('id').data.all()])
        
