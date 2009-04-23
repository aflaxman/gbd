from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import simplejson as json

import dismod3.models.fields as fields

class DataAdmin(admin.ModelAdmin):
    list_display = ('id', 'condition', 'data_type', 'gbd_region',
                    'region', 'sex', 'time_start', 'time_end',
                    'age_start', 'age_end', 'value',
                    'standard_error')
    list_filter = ['condition', 'gbd_region', 'sex',]
    search_fields = ['region', 'id',]

class Data(models.Model):
    """
    Model for generic piece of dismod data, which is semi-structured,
    and will usually include a condition of interest, data_type,
    region, sex, time range, age range, value, and standard error

    any additional information should be stored in a dictionary that
    is saved in the params_json field.

    data_type will take values 'incidence data', 'prevalence data', ...
    """
    
    # information which is included implicitly in a row of Theo's table
    condition = models.CharField(max_length=200)
    data_type = fields.RateTypeField()

    gbd_region = models.CharField(max_length=200)
    region = models.CharField(max_length=200)

    sex = fields.SexField()

    time_start = models.IntegerField()
    time_end = models.IntegerField()
    age_start = models.IntegerField()
    age_end = models.IntegerField()

    value = models.FloatField()
    standard_error = models.FloatField()

    params_json = models.TextField(default=json.dumps({}))

    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3' 

    def __init__(self, *args, **kwargs):
        super(Data, self).__init__(*args, **kwargs)
        try:
            self.params = json.loads(self.params_json)
        except ValueError:
            self.params = {}

    def cache_params(self):
        """
        store the params dict as json text
        this must be called before data.save()
        to preserve any changes to params dict
        """
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s, %s, %s, %s' \
               % (self.condition, self.region, self.time_str(),
                  self.sex, self.age_str(),
                  self.value_str())

    def get_absolute_url(self):
        return reverse("new_dm3.views.data_show", args=(self.id,))

    def get_edit_url(self):
        return '/admin/dismod3/data/%i/' % self.id

    def age_str(self):
        if self.age_end == fields.MISSING:
            return '%d+' % self.age_start
        elif self.age_start == self.age_end:
            return '%d' % self.age_start
        else:
            return '%d-%d' % (self.age_start, self.age_end)

    def time_str(self):
        if self.time_start == self.time_end:
            return '%d' % self.time_start
        else:
            return '%d-%d' % (self.time_start, self.time_end)

    def value_str(self):
        return '%3f (%3f, %3f)' % (self.value,
                                   max(0.,self.value - 2. * self.standard_error),
                                   self.value + 2. * self.standard_error)

    def age_weights(self):
        """
        return a population-by-age weight vector for the country and year(s)
        that this data comes from.

        The returned vector has length self.age_end + 1 - self.age_start,
        and entry i is the proportion of the population with age =
        (self.age_start + i)

        Cache this vector in self.params to avoid repeatedly making
        the database queries required to compute it.
        """
        import numpy as np

        if self.params.has_key('age_weights'):
            return self.params['age_weights']

        # calculate age_weights, and cache it for future reference
        if self.age_end == fields.MISSING:
            self.age_end = MAX_AGE-1

        if rate.age_end < rate.age_start:
            raise ValueError('rate %d has age_end < age_start' % rate.id)
        elif rage.age_end == rate.age_start:
            # don't need to look in database for a single year
            pop_vals = [ 1. ]
        else:
            a = range(self.age_start, self.age_end+1)
            total = np.zeros(len(a))

            relevant_populations = Data.objects.filter(data_type='population', region=self.region, sex=self.sex,
                                                       year__gte=self.time_start, year__lte=self.time_end)
            if relevant_populations.count() == 0:
                print "WARNING: Population for %s-%d-%d-%s not found, using uniform distribution instead of age-weighted distribution (rate_id=%d)" \
                      % (rate.region, rate.time_start, rate.time_end, rate.sex, rate.pk)
                pop_vals = np.ones(len(a))
            else:
                for population in relevant_populations:
                    a0 = population.age_start - self.age_start
                    a1 = population.age_end + 1 - self.age_start
                    total[min(0,a0):max(len(a),a1)] += population.value / (a1 - a0)

                pop_vals = total/(rate.time_end + 1. - rate.time_start)

        self.params['population'] = list(pop_vals)
        self.save()

        return self.params['population']
