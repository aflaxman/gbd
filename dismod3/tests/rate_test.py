from base_test import DisModTestCase
from django.test.client import Client
import simplejson as json

from dismod3.models import *
from dismod3.views import *

class RateTestCase(DisModTestCase):
    def setUp(self):
        from dismod3.models import Rate
        self.rate = Rate.objects.all()[0]
        self.create_users()

    # unit tests
    def test_str(self):
        s = str(self.rate)
        self.assertTrue(isinstance(s,str))

        s = self.rate.get_absolute_url()
        self.assertTrue(isinstance(s,str))


    def test_urls(self):
        url = self.rate.get_absolute_url()
        self.assertTrue(isinstance(url,str))

        url = self.rate.get_edit_url()
        self.assertTrue(isinstance(url,str))

    def test_params(self):
        self.assertTrue(isinstance(self.rate.params, dict))

    # functional tests
    def test_plot(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get('/rate/plot/disease_1-rate_prevalence.png')
        self.assertEquals(response.status_code, 200)

        response = c.get('/rate/plot/disease_1-rate_prevalence-region_1.pdf')
        self.assertEquals(response.status_code, 200)

        response = c.get('/rate/plot/rate_prevalence-region_1-sex_total.svg')
        self.assertEquals(response.status_code, 200)

    def test_show(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get(self.rate.get_absolute_url())
        self.assertTemplateUsed(response, 'rate/show.html')

    def test_index(self):
        c = Client()
        c.login(username='red', password='red')
        response = c.get('/rate/')
        self.assertTemplateUsed(response, 'rate/index.html')
        
        response = c.post('/rate/', {'tab_separated_values': ''})
        self.assertTemplateUsed(response, 'rate/index.html')

        # now do it right, and make sure that asrfs are added
        asrf_cnt = AgeSpecificRateFunction.objects.count()
        response = c.post('/rate/', {'tab_separated_values': \
        'Disease\tRegion\tRate Type\tSex\tCountry\tAge Start\tAge End\tEstimate Year Start\tEstimate Year End\tRate\tNumber of Subjects\tStandard Error\nCannabis Dependence\tWorld\tPrevalence\tTotal\tCanada\t15\t24\t2005\t2005\t.5\t1000\t.01'})

        id_str = view_utils.objects_to_id_str(AgeSpecificRateFunction.objects.all()[asrf_cnt:])
        redirect_url = '/age_specific_rate_function/%s' % id_str
        self.assertRedirects(response, redirect_url)
        assert asrf_cnt < AgeSpecificRateFunction.objects.count()
