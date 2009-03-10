from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import simplejson as json

import fields

class RateAdmin(admin.ModelAdmin):
    list_display = ('id', 'disease', 'rate_type', 'region',
                    'country', 'epoch_start', 'epoch_end', 'sex',
                    'age_start', 'age_end', 'numerator',
                    'denominator', 'rate')
    list_filter = ['disease', 'region', 'rate_type', 'sex',]
    search_fields = ['country', 'id',]

class Rate(models.Model):
    """
    Model for a rate for a particular disease, region, rate_type, sex, and age range

    I've attempted to put all of the essential stuff in the model, but any additional
    information should be stored in a dictionary that is saved in the params_json field.

        rate_type will take values 'incidence data', 'prevalence data', ...
    """
    # information which is included implicitly in a row of Theo's table
    disease = models.ForeignKey('Disease')
    region = models.ForeignKey('Region')
    rate_type = fields.RateTypeField()
    sex = fields.SexField()
    country = models.CharField(max_length=200)
    age_start = models.IntegerField()
    age_end = models.IntegerField()
    epoch_start = models.IntegerField()
    epoch_end = models.IntegerField()
    numerator = models.IntegerField()
    denominator = models.IntegerField()
    rate = models.FloatField(default=0.)
    params_json = models.TextField(default=json.dumps({}))

    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3' 

    def __init__(self, *args, **kwargs):
        super(Rate, self).__init__(*args, **kwargs)
        try:
            self.params = json.loads(self.params_json)
        except ValueError:
            self.params = {}

    def save(self, force_insert=False, force_update=False):
        # store the params dict as json text
        self.params_json = json.dumps(self.params)

        # calculate the rate, for sorting purposes
        self.rate = float(self.numerator)/float(self.denominator)        
        super(Rate,self).save(force_insert, force_update)

    def __unicode__(self):
        return "%s, %s, %d-%d, %s, %d-%d, %d/%d" % (self.disease, self.country, self.epoch_start, self.epoch_end,
                                                    self.sex, self.age_start, self.age_end,
                                                    self.numerator, self.denominator)

    def get_absolute_url(self):
        return reverse("dismod3.views.rate_show", args=(self.id,))

    def get_edit_url(self):
        return "/admin/dismod3/rate/%i" % self.id

    def pretty_rate(self):
        ci = self.ci()
        return "%.3f (%.3f, %.3f)" % (self.rate, max(0,ci[0]), ci[1])


    def ci(self):
        """
        estimate [p0,p1] such that Pr[B(r.denominator, p0 (1) ) < (>) r.numerator] <= 0.025 (>= 0.975)
        """
        import numpy as np
        
        a = float(self.numerator)
        b = float(self.denominator)
        var = (a * b) / ((a + b)**2. * (a + b + 1.))
        std = np.sqrt(var)
        return (a/b - 2*std), (a/b + 2*std)


    def population(self):
        """
        return a population-by-age vector for the country and year(s)
        that this rate was sampled from.

        The returned vector has length self.age_end - self.age_start,
        and entry i is the proportion of the population with age =
        (self.age_start + i)

        Cache this vector in self.params to avoid repeatedly making
        the database queries required to compute it.
        """
        from dismod3.bayesian_models import probabilistic_utils

        if self.params.has_key('population'):
            return self.params['population']
        else:
            self.params['population'] = list(probabilistic_utils.population_during(self))
            self.save()
            return self.params['population']
