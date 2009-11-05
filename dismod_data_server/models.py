from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User

import simplejson as json

import gbd.fields
import gbd.settings
import dismod3

def debug(string):
    """ Print string, or output it in the appropriate way for the
    environment (i.e. don't output it at all on production server).
    """
    if gbd.settings.DEBUG_TO_STDOUT:
        import sys
        print string
        sys.stdout.flush()

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
    region, sex, year range, age range, value, and standard error.

    Any additional information should be stored in a dictionary that
    is saved in the params_json field.

    For more details, see ``dismod_data_json.html``.
    """
    
    condition = models.CharField(max_length=200)
    data_type = gbd.fields.DataTypeField()

    gbd_region = models.CharField(max_length=200)
    region = models.CharField(max_length=200)

    sex = gbd.fields.SexField()

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
            debug('WARNING: could not load params_json for DiseaseModel %d' % self.id)
            self.params = {}

    def cache_params(self):
        """ Store the params dict as json text

        Notes
        -----
        This must be called before data.save() to preserve any
        changes to params dict.

        I do it this way, instead of automatically in the save method
        to permit direct json editing in the admin interface.

        Example
        -------
        >>> d = Data(condition='ats_use', data_type='prevalence', \
                     region='Canada', gbd_region='North America, High Income', \
                     sex='male', year_start=1990, year_end=1990, \
                     age_start=15, age_end=15, \
                     value=.01, standard_error=.005)
        >>> d.save()
        >>> assert d.params == {}, 'params should not yet be cached'
        >>> d.cache_params()
        >>> d.save()
        >>> assert d.params != {}, 'params should now be cached'
        """

        self.params.update(id=self.id,
                           condition=self.condition,
                           data_type=self.data_type,
                           gbd_region=self.gbd_region,
                           region=self.region,
                           sex=self.sex,
                           year_start=self.year_start,
                           year_end=self.year_end,
                           age_start=self.age_start,
                           age_end=self.age_end,
                           value=self.value,
                           standard_error=self.standard_error)
        self.params_json = json.dumps(self.params)

    def __unicode__(self):
        return '%s, %s, %s, %s, %s, %s, %s' \
               % (self.condition, self.data_type, self.region, self.year_str(),
                  self.get_sex_display(), self.age_str(),
                  self.value_str())

    def get_absolute_url(self):
        return reverse('gbd.dismod_data_server.views.data_show', args=[self.id])

    def age_str(self):
        """ Return a pretty string describing the age range of this data
        point.
        """
        
        if self.age_end == dismod3.MISSING:
            return '%d+' % self.age_start
        elif self.age_start == self.age_end:
            return '%d' % self.age_start
        else:
            return '%d-%d' % (self.age_start, self.age_end)

    def year_str(self):
        """ Return a pretty string describing the time range of this data
        point.
        """
        
        if self.year_start == self.year_end:
            return '%d' % self.year_start
        else:
            return '%d-%d' % (self.year_start, self.year_end)

    def value_str(self):
        """ Return a pretty string describing the value and uncertainty of this data
        point.
        """
        
        if self.standard_error == 0. or self.standard_error == dismod3.MISSING:
            return '%f' % self.value
        else:
            return '%f (%f, %f)' % (self.value,
                                    max(0., self.value - 1.96 * self.standard_error),
                                    self.value + 1.96 * self.standard_error)

    def age_weights(self):
        """ Return a population-by-age weight vector for the country and
        year(s) that this data comes from.
        
        The returned vector has length self.age_end + 1 - self.age_start,
        and entry i is the proportion of the population with age =
        (self.age_start + i).
        
        Cache this vector in self.params to avoid repeatedly making
        the database queries required to compute it.
        """
        
        if not self.params.has_key('age_weights'):
            self.calculate_age_weights()

        return self.params['age_weights']

    def to_dict(self):
        return self.params

    def calculate_age_weights(self):
        """ Calculate and cache age_weights vector in self.params to avoid
        repeatedly making the database queries required to compute it.
        """

        import numpy as np
        from population_data_server.models import Population

        # deal with 'missing' data at end of age interval
        if self.age_end == dismod3.MISSING:
            self.age_end = dismod3.MAX_AGE - 1

        # sanity check
        if self.age_end < self.age_start:
            raise ValueError('Data %d has age_end < age_start' % self.id)
        elif self.age_end == self.age_start:
            # don't need to look in database for a single year
            pop_vals = [ 1. ]
        else:
            a = range(self.age_start, self.age_end+1)
            relevant_populations = Population.objects.filter(region=self.region,
                                                             sex=self.sex,
                                                             year__gte=self.year_start,
                                                             year__lte=self.year_end)
            if relevant_populations.count() == 0:
                debug(("WARNING: Population for %s-%d-%d-%s not found, "
                       + "using uniform distribution instead of age-weighted "
                       + "distribution (Data_id=%d)" )
                      % (self.region, self.year_start, self.year_end, self.sex, self.id))
                total = np.ones(len(a))
            else:
                total = np.zeros(len(a))
                for population in relevant_populations:
                    M,C = population.gaussian_process()
                    total += M(a)

            pop_vals = np.maximum(dismod3.NEARLY_ZERO, total)
            pop_vals /= sum(pop_vals)

        self.params['age_weights'] = list(pop_vals)
        self.cache_params()
        self.save()

        return self.params['age_weights']

    def relevant_to(self, type, region, year, sex):
        """ Determine if this data is relevant to the requested
        type, region, year, and sex"""
        return dismod3.relevant_to(self.to_dict(), type, region, year, sex)


class DiseaseModelAdmin(admin.ModelAdmin):
    list_display = ('id', 'condition',
                    'region', 'sex', 'year',)
    list_filter = ['condition', 'region', 'sex', 'year']
    search_fields = ['region', 'id',]

def create_disease_model(dismod_dataset_json, creator):
    """ Turn a dismod_dataset json into a honest DiseaseModel object and
    save it in the database.
    """

    model_dict = json.loads(dismod_dataset_json)

    params = model_dict['params']
    args = {}
    args['region'] = params['region']
    args['year'] = params['year']
    args['sex'] = params['sex']
    args['condition'] = params['condition']
    args['creator'] = creator

    dm = DiseaseModel.objects.create(**args)
    for d_data in model_dict['data']:
        dm.data.add(d_data['id'])

    for key in params:
        p, flag = dm.params.get_or_create(key=key)
        p.json = json.dumps(params[key])
    return dm

class DiseaseModelParameter(models.Model):
    """ Any sort of semi-structured data that is associated with a
    disease model.

    Used for holding priors, initial values, model fits, etc.
    """
    region = models.CharField(max_length=200, blank=True)
    sex = gbd.fields.SexField(blank=True)
    year = models.CharField(max_length=200, blank=True)
    type = gbd.fields.DataTypeField(blank=True)

    key = models.CharField(max_length=200)
    json = models.TextField(default=json.dumps({}))

    def __unicode__(self):
        if self.region and self.sex and self.year and self.type:
            return '%d: %s (%s, %s, %s, %s)' \
                   % (self.id, self.key, self.region, self.get_sex_display(), self.year, self.get_type_display())
        else:
            return '%d: %s' % (self.id, self.key)
    
class DiseaseModel(models.Model):
    """ Model for a collection of dismod data, together with priors and
    all other relevant parameters.
    
    The params hash can also include the parameter estimates produced
    by fitting the model.
    """

    condition = models.CharField(max_length=200)
    region = models.CharField(max_length=200)
    sex = gbd.fields.SexField()
    year = models.CharField(max_length=200)

    data = models.ManyToManyField(Data)
    params = models.ManyToManyField(DiseaseModelParameter)

    creator = models.ForeignKey(User)
    # TODO: add notes textfield, add data created, modified fields

    def __unicode__(self):
        return '%s, %s, %s, %s' \
               % (self.condition, self.region, self.get_sex_display(), self.year,)

    def get_absolute_url(self):
        return reverse('gbd.dismod_data_server.views.dismod_show', args=(self.id,))

    def to_json(self, filter_args={}):
        """ Return a dismod_dataset json corresponding to this model object

        See ``dismod_data_json.html`` for details.
        """
        param_dict = {}

        for p in self.params.filter(**filter_args):
            if p.type and p.region and p.sex and p.year:
                if not param_dict.has_key(p.key):
                    param_dict[p.key] = {}
                param_dict[p.key][dismod3.gbd_key_for(p.type,p.region,p.year,p.sex)] = json.loads(p.json)
            else:
                param_dict[p.key] = json.loads(p.json)
        # include params for all regions as well, if params were filtered above
        if len(filter_args) > 0:
            for p in self.params.filter(region=''):
                if param_dict.has_key(p.key):
                    continue
                param_dict[p.key] = json.loads(p.json)

        param_dict.update(id=self.id,
                          condition=self.condition,
                          sex=self.sex,
                          region=self.region,
                          year=self.year)

        return json.dumps({'params': param_dict,
                           'data': [d.params for d in self.data.all()],
                           'id': self.id})
#                          sort_keys=True, indent=2)
