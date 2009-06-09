from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import simplejson as json

import fields

class PopulationAdmin(admin.ModelAdmin):
    list_display  = ('id', 'region', 'country', 'sex', 'year',)
    list_filter   = ['sex', 'region', ]
    search_fields = ['country', 'sex', 'year',]

class Population(models.Model):
    """
    most of these fields are self-explanatory, but
    data_json is not.  this is an place to store
    any relevant information in a dictionary. the
    dictionary is expected to contain:
        data['mesh'] = age-points where the population size has been measured
        data['vals'] = values at each mesh point
    """
    region = models.ForeignKey('Region')
    country = models.CharField(max_length=200)
    iso_code = models.CharField(max_length=10)
    year = models.IntegerField()
    sex = fields.SexField()
    data_json = models.TextField(blank=True)

    def _get_data(self):
        return json.loads(self.data_json)
    data = property(_get_data)

    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3' 

    def __unicode__(self):
        return "%s, %s, %s" % (self.country, self.year, self.get_sex_display(),)

    def get_absolute_url(self):
        return reverse("dismod3.views.population_show", args=(self.id,))

    def get_png_url(self):
        return reverse("dismod3.views.population_plot", args=(self.id,'png'))        

    def gaussian_process(self):
        """
        return a PyMC Gaussian Process mean and covariance to interpolate
        the population-by-age mesh/value data
        """
        from pymc import gp
        from dismod3.bayesian_models import probabilistic_utils

        M, C = probabilistic_utils.uninformative_prior_gp(c=0.,  diff_degree=2., amp=10., scale=200.)
        gp.observe(M, C, self.data['mesh'], self.data['vals'], 0.0)
    
        return M, C

