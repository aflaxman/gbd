from base_test import DisModTestCase
from django.test.client import Client

import simplejson as json

from dismod3.models import Population

class PopulationTestCase(DisModTestCase):
    def setUp(self):
        self.pop = Population.objects.get(country='Australia', year=2009, sex='male')
        self.create_users()

    # unit tests
    def test_strs(self):
        s = str(self.pop)
        self.assertTrue(isinstance(s,str))

        s = self.pop.get_absolute_url()
        self.assertTrue(isinstance(s,str))

        s = self.pop.get_png_url()
        self.assertTrue(isinstance(s,str))

    def test_url(self):
        url = self.pop.get_absolute_url()
        self.assertTrue(isinstance(url,str))

    def testGP(self):
        import numpy as np
        M, C = self.pop.gaussian_process()
        self.assertEquals(len(M([0,10,100])), 3)
    

    # functional tests
    def test_population_fixtures(self):
        pop = Population.objects.get(country='Australia', year=2009, sex='male')
        self.assertEquals(pop, self.pop)

    def test_plot(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get(self.pop.get_png_url())
        self.assertEquals(response.status_code, 200)

    def test_show(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get(self.pop.get_absolute_url())
        self.assertTemplateUsed(response, 'population/show.html')

    def test_redirect(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get(self.pop.get_absolute_url() + '/png')
        self.assertRedirects(response, 'population/%d.png' % self.pop.id)

        response = c.get(self.pop.get_absolute_url() + '/csv')
        self.assertRedirects(response, 'population/%d.csv' % self.pop.id)

        response = c.get(self.pop.get_absolute_url() + '/json')
        self.assertRedirects(response, 'population/%d.json' % self.pop.id)

        response = c.get(self.pop.get_absolute_url() + '/no-nexistent_command')
        self.assertEqual(response.status_code, 404)
        
