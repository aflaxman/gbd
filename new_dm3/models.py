from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import simplejson as json

import dismod3.models.fields as fields

class DataAdmin(admin.ModelAdmin):
    list_display = ('id', 'condition', 'data_type', 'gbd_region',
                    'region', 'sex', 'year_start', 'year_end',
                    'age_start', 'age_end', 'value',
                    'standard_error')
    list_filter = ['condition', 'gbd_region', 'sex',]
    search_fields = ['region', 'id',]

class Data(models.Model):
    """
    Model for generic piece of dismod data, which is semi-structured,
    and will usually include a condition of interest, data_type,
    region, sex, year range, age range, value, and standard error

    any additional information should be stored in a dictionary that
    is saved in the params_json field.

    data_type will take values 'incidence data', 'prevalence data', ...
    """
    
    # information which is included implicitly in a row of Theo's table
    condition = models.CharField(max_length=200)
    data_type = models.CharField(max_length=200)

    gbd_region = models.CharField(max_length=200)
    region = models.CharField(max_length=200)

    sex = fields.SexField()

    year_start = models.IntegerField()
    year_end = models.IntegerField()
    age_start = models.IntegerField()
    age_end = models.IntegerField()

    value = models.FloatField()
    standard_error = models.FloatField()

    params_json = models.TextField(default=json.dumps({}))

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

        do it this way, instead of automatically in the save method to
        permit direct json editing in the admin interface
        """
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s, %s, %s, %s' \
               % (self.condition, self.region, self.year_str(),
                  self.get_sex_display(), self.age_str(),
                  self.value_str())

    def get_absolute_url(self):
        return reverse("new_dm3.views.data_show", args=(self.id,))

    def get_edit_url(self):
        return '/admin/new_dm3/data/%i/' % self.id

    def age_str(self):
        if self.age_end == fields.MISSING:
            return '%d+' % self.age_start
        elif self.age_start == self.age_end:
            return '%d' % self.age_start
        else:
            return '%d-%d' % (self.age_start, self.age_end)

    def year_str(self):
        if self.year_start == self.year_end:
            return '%d' % self.year_start
        else:
            return '%d-%d' % (self.year_start, self.year_end)

    def value_str(self):
        if self.standard_error == 0.:
            return '%.2f' % self.value
        else:
            return '%.2f (%.2f, %.2f)' % (self.value,
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
        if not self.params.has_key('age_weights'):
            self.calculate_age_weights()

        return self.params['age_weights']

    def calculate_age_weights(self):
        """
        Calculate and cache age_weights vector in self.params to avoid repeatedly making
        the database queries required to compute it.
        """
        import numpy as np
        from dismod3.models import Population

        # deal with 'missing' data about end of age interval
        if self.age_end == fields.MISSING:
            self.age_end = MAX_AGE-1

        # sanity check
        if self.age_end < self.age_start:
            raise ValueError('Data %d has age_end < age_start' % self.id)
        elif self.age_end == self.age_start:
            # don't need to look in database for a single year
            pop_vals = [ 1. ]
        else:
            a = range(self.age_start, self.age_end+1)
            relevant_populations = Population.objects.filter(country=self.region, sex=self.sex,
                                                             year__gte=self.year_start, year__lte=self.year_end)
            if relevant_populations.count() == 0:
                print "WARNING: Population for %s-%d-%d-%s not found, using uniform distribution instead of age-weighted distribution (Data_id=%d)" \
                      % (self.region, self.year_start, self.year_end, self.sex, self.pk)
                pop_vals = np.ones(len(a))
            else:
                total = np.zeros(len(a))
                for population in relevant_populations:
                    M,C = population.gaussian_process()
                    total += M(a)

                pop_vals = np.maximum(NEARLY_ZERO, total)
                pop_vals /= sum(pop_vals)

        self.params['age_weights'] = list(pop_vals)
        self.cache_params()
        self.save()

        return self.params['age_weights']




class DiseaseModelAdmin(admin.ModelAdmin):
    list_display = ('id', 'condition',
                    'region', 'sex', 'year',)
    list_filter = ['condition', 'region', 'sex', 'year']
    search_fields = ['region', 'id',]

class DiseaseModel(models.Model):
    """
    Model for a collection of dismod data, together with priors and
    any other relevant parameters
    
    also may include the results of fitting the model
    """

    condition = models.CharField(max_length=200)
    region = models.CharField(max_length=200)
    sex = fields.SexField()
    year = models.CharField(max_length=200)

    data = models.ManyToManyField(Data)

    params_json = models.TextField(default=json.dumps({}))

    def __init__(self, *args, **kwargs):
        super(DiseaseModel, self).__init__(*args, **kwargs)
        try:
            self.params = json.loads(self.params_json)
        except ValueError:
            self.params = {}

    def cache_params(self):
        """
        store the params dict as json text
        
        this must be called before data.save()
        to preserve any changes to params dict

        do it this way, instead of automatically in the save method to
        permit direct json editing in the admin interface
        """
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s, %s' \
               % (self.condition, self.region, self.get_sex_display(), self.year,)

    def get_absolute_url(self):
        return reverse("new_dm3.views.disease_model_show", args=(self.id,))

    def get_edit_url(self):
        return '/admin/dismod3/diseasemodel/%i/' % self.id
