from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import simplejson as json


SEX_CHOICES = [
    ('male', _('Male')),
    ('female', _('Female')),
    ('total', _('Total')),
]

class SexField(models.CharField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('max_length', 10)
        kwargs.setdefault('choices', SEX_CHOICES)

        super(SexField, self).__init__(*args, **kwargs)

    def get_internal_type(self):
        return "CharField"


class PopulationAdmin(admin.ModelAdmin):
    list_display  = ('id', 'region', 'sex', 'year',)
    list_filter   = ['sex', 'region', ]
    search_fields = ['region', 'sex', 'year',]

class Population(models.Model):
    """ Model for Population Data
    Parameters
    ----------
    region : str
    year : int
    sex : str
    params_json : str
        this is an place to store any relevant information in a
        dictionary. the dictionary is expected to contain::

            data['mesh'] = age-points where the population size has been measured
            data['vals'] = values at each mesh point
            super(Population, self).__init__(*args, **kwargs)
    """
    region = models.CharField(max_length=200)
    year = models.IntegerField()
    sex = SexField()
    params_json = models.TextField(default=json.dumps({}))

    def __init__(self, *args, **kwargs):
        super(Population, self).__init__(*args, **kwargs)
        try:
            self.params = json.loads(self.params_json)
        except ValueError:
            self.params = {}

    def cache_params(self):
        """
        store the params dict as json text
        
        this must be called before population.save()
        to preserve any changes to params dict

        do it this way, instead of automatically in the save method to
        permit direct json editing in the admin interface
        """
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s' % (self.region, self.year, self.get_sex_display(),)

    def get_absolute_url(self):
        return reverse('gbd.population_data_server.views.population_show', args=(self.id,))

    def gaussian_process(self):
        """ return a PyMC Gaussian Process mean and covariance to interpolate
        the population-by-age mesh/value data
        """
        from pymc import gp
        from dismod3.bayesian_models import probabilistic_utils

        M, C = probabilistic_utils.uninformative_prior_gp(c=0.,  diff_degree=2., amp=10., scale=200.)
        gp.observe(M, C, self.params['mesh'], self.params['vals'], 0.0)
    
        return M, C


