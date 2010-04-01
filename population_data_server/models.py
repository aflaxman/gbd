from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import numpy as np
from pymc import gp
import simplejson as json

import gbd.fields
from dismod3.utils import MAX_AGE

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

            params['mesh'] = midpoints of age intervals for population size
            params['vals'] = average value for each age interval
            params['interval_start'] = left endpoint of age intervals (optional)
            params['interval_length'] = width of age intervals (optional)
    """
    region = models.CharField(max_length=200)
    year = models.IntegerField()
    sex = gbd.fields.SexField()
    params_json = models.TextField(default=json.dumps({}))

    def __init__(self, *args, **kwargs):
        super(Population, self).__init__(*args, **kwargs)
        try:
            self.params = json.loads(self.params_json)
        except ValueError:
            self.params = {}

    def cache_params(self):
        """ Store the params dict as json text

        Notes
        -----
        This must be called before population.save() to preserve any
        changes to params dict

        I do it this way, instead of automatically in the save method
        to permit direct json editing in the admin interface
        """
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s' % (self.region, self.year, self.get_sex_display(),)

    def get_absolute_url(self):
        return reverse('gbd.population_data_server.views.population_show', args=(self.id,))

    def interpolate(self, age_range):
        #M,C = self.gaussian_process()
        #return M(a)
        from dismod3.utils import interpolate
        self.params['mesh'][0] = 0.0
        wts = interpolate(self.params['mesh'] + [ MAX_AGE ], self.params['vals'] + [ 0. ], age_range)
        return wts


#     def gaussian_process(self):
#         """ return a PyMC Gaussian Process mean and covariance to interpolate
#         the population-by-age mesh/value data
#         """
#         # TODO: make this evaluate the function on arange(MAX_AGE) and store the results in the db for better performance
#         M, C = uninformative_prior_gp(c=0.,  diff_degree=2., amp=10., scale=200.)
#         gp.observe(M, C, self.params['mesh'] + [ MAX_AGE ], self.params['vals'] + [ 0. ], 0.0)
    
#         return M, C

# def const_func(x, c):
#     """ A constant function, f(x) = c

#     To be used as a non-informative prior on a Gaussian process.

#     Example
#     -------
#     >>> const_func([1,2,3], 17.0)
#     array([ 17., 17., 17.])
#     """
#     return np.zeros(np.shape(x)) + c

# def uninformative_prior_gp(c=-10.,  diff_degree=2., amp=100., scale=200.):
#     """ Uninformative Mean and Covariance Priors
#     Parameters
#     ----------
#     c : float, the prior mean
#     diff_degree : float, the prior on differentiability (2 = twice differentiable?)
#     amp : float, the prior on the amplitude of the Gaussian Process
#     scale : float, the prior on the scale of the Gaussian Process

#     Results
#     -------
#     M, C : mean and covariance objects
#       this constitutes an uninformative prior on a Gaussian Process
#       with a euclidean Matern covariance function
#     """
#     M = gp.Mean(const_func, c=c)
#     C = gp.Covariance(gp.matern.euclidean, diff_degree=diff_degree,
#                       amp=amp, scale=scale)

#     return M,C
