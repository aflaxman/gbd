from base_test import DisModTestCase
from django.test.client import Client
import simplejson as json

class RateTestCase(DisModTestCase):
    def setUp(self):
        from dismod3.models import Rate
        self.rate = Rate.objects.all()[0]

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
        response = c.get('/rate/plot/disease_1-rate_prevalence.png')
        self.assertEquals(response.status_code, 200)

        response = c.get('/rate/plot/disease_1-rate_prevalence-region_1.pdf')
        self.assertEquals(response.status_code, 200)

        response = c.get('/rate/plot/rate_prevalence-region_1-sex_total.svg')
        self.assertEquals(response.status_code, 200)

    def test_show(self):
        c = Client()
        response = c.get(self.rate.get_absolute_url())
        self.assertTemplateUsed(response, 'rate/show.html')

    def test_index(self):
        c = Client()
        response = c.get('/rate/')
        self.assertTemplateUsed(response, 'rate/index.html')
        
        response = c.post('/rate/', {'tab_separated_values': ''})
        self.assertTemplateUsed(response, 'rate/index.html')

        response = c.post('/rate/', {'tab_separated_values': \
        """Disease\tRegion\tRate Type\tSex\tCountry\tAge Start\tAge End\tEstimate Year Start\tEstimate Year End\tRate\tNumber of Subjects\tStandard Error
        of the multiline system
        """})
        self.assertTemplateUsed(response, 'age_specific_rate_function/show.html')
