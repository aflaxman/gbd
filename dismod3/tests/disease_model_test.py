from base_test import DisModTestCase
from django.test.client import Client
import simplejson as json

from dismod3.models import *
from dismod3.views import *

class DiseaseModelTestCase(DisModTestCase):
    def setUp(self):
        self.disease = Disease.objects.get(name='Cannabis Use')
        self.region = Region.objects.get(name='Australasia')
        self.rate = Rate.objects.all()[0]
        self.asrf = AgeSpecificRateFunction.objects.all()[1]
        self.dm = DiseaseModel.objects.all()[0]

    ############
    # unit tests
    #
    def test_save(self):
        dm = DiseaseModel(name='Test Model')
        dm.save()

    def test_strs(self):
        s = str(self.dm)
        self.assertNotEqual(s, '')

        s = self.dm.get_absolute_url()
        self.assertNotEqual(s, '')

        s = self.dm.get_asrf_id_str()
        self.assertNotEqual(s, '')

    def test_fit(self):
        #from dismod3.bayesian_models import fit_rate_function
        #fit_rate_function.mcmc_fit(self.asrf, speed='testing fast')
        pass

    ##################
    # functional tests
    #
    def test_index_view(self):
        pass
        #c = Client()

        #response = c.get('/disease_model/')
        #self.assertTemplateUsed(response, 'age_specific_rate_function/index.html')

    def test_show_view(self):
        """
        the show view supports arguments about subplot size and also
        about the axis, for zooming
        """

        c = Client()
        response = c.get('/disease_model/%d' % self.dm.id)
        self.assertTemplateUsed(response, 'disease_model/show.html')
