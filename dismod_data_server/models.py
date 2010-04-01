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

            if self.data_type.startswith('prevalence') or \
                   self.data_type.startswith('incidence'): # use population structure for age weights
                # SPECIAL CASE: sex == 'all' will be applied to males and females, so use pop structure for total
                if self.sex == 'all':
                    sex = 'total'
                else:
                    sex = self.sex
                relevant_populations = Population.objects.filter(region=self.region,
                                                                 sex=sex,
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
                        total += population.interpolate(a)
            else: # don't use population structure for age weights
                total = np.ones(len(a))
                
            pop_vals = np.maximum(dismod3.NEARLY_ZERO, total)
            pop_vals /= sum(pop_vals)

        self.params['age_weights'] = list(pop_vals)
        self.cache_params()
        self.save()

        return self.params['age_weights']

    def calculate_covariate(self, covariate_type):
        """ Calculate and cache specified covariate in self.params to avoid
        repeatedly making the database queries required to compute it.
        """

        import numpy as np
        from covariate_data_server.models import Covariate
        from gbd.dismod3.utils import clean

        # TODO: allow a way for one db query to calculate covariates for many data points
        covariates = Covariate.objects.filter(
            type__slug=covariate_type,
            sex=self.sex,
            country_year__in=['%s-%d' % (self.region, y) for y in [gbd.fields.ALL_YEARS] + range(self.year_start,self.year_end+1)])
        if len(covariates) == 0:
            debug(("WARNING: Covariate %s not found for %s %s-%s, "
                   + "(Data_id=%d)" )
                  % (covariate_type, self.sex, self.region, self.year_str(), self.id))

        else:
            self.params[clean(covariate_type)] = np.mean([c.value for c in covariates])
            self.cache_params()
            self.save()
            debug('updated %s %s %s-%s, (Data_id=%d)' % (covariate_type, self.sex, self.region, self.year_str(), self.id))
    def relevant_to(self, type, region, year, sex):
        """ Determine if this data is relevant to the requested
        type, region, year, and sex"""
        return dismod3.relevant_to(self.to_dict(), type, region, year, sex)


class DiseaseModelAdmin(admin.ModelAdmin):
    list_display = ('id', 'condition',
                    'region', 'sex', 'year',)
    list_filter = ['condition', 'region', 'sex', 'year']
    search_fields = ['region', 'id',]

class DiseaseModelParameterAdmin(admin.ModelAdmin):
    list_display = ('id', 'region', 'sex', 'year', 'type', 'key')
    list_filter = ('key', 'sex', 'year', 'type', 'region', )
    search_fields = ['key', 'id',]

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
        if params[key]:
            p, flag = dm.params.get_or_create(key=key)
            p.json = json.dumps(params[key])
            p.save()
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

        filter_args : dict

        Example
        -------
        >> dm = DiseaseModel.objects.get(id=1)
        >> dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'}))

        Notes
        -----
        This {'region': 'none'} business is a tricky optimization, so
        that the DiseaseModel itself is small, and the large amounts
        of generated data are stored in DiseaseModelParameter objects.

        {'region': 'none'} says don't merge in any of the region specific
        diseasemodelparameters when you are converting the disease
        model to json.
        """
        param_dict = {}
        for p in self.params.filter(**filter_args):
            if p.type and p.region and p.sex and p.year:
                if not param_dict.has_key(p.key):
                    param_dict[p.key] = {}
                param_dict[p.key][dismod3.gbd_key_for(p.type,p.region,p.year,p.sex)] = json.loads(p.json)
            else:
                try:
                    param_dict[p.key] = json.loads(p.json)
                except ValueError:
                    # skip bad json, it sometimes happens, for unknown reasons (HTTP glitches?)
                    pass
        # include params for all regions as well, if params were filtered above
        if len(filter_args) > 0:
            for p in self.params.filter(region=''):
                if param_dict.has_key(p.key):
                    continue
                try:
                    param_dict[p.key] = json.loads(p.json)
                except ValueError:
                    # skip bad json, it sometimes happens, for unknown reasons (HTTP glitches?)
                    pass

        param_dict.update(id=self.id,
                          condition=self.condition,
                          sex=self.sex,
                          region=self.region,
                          year=self.year)

        return json.dumps({'params': param_dict,
                           'data': [d.params for d in self.data.all()],
                           'id': self.id})

    def country_level_covariates(self):
        from gbd.covariate_data_server.models import CovariateType, Covariate
        cov_dict = {}
        for ct in CovariateType.objects.all():
            cov_dict[ct.slug] = {
                'rate': dict(value=0, default=0),
                'error': dict(value=0, default=0),
                'value': dict(value='0', default='0'),  # value must be a string be a string
                'range': [0, 10^6],
                'category': ['', ''],
                'defaults': dict([[c.iso3, c.value] for c in ct.covariate_set.all()])
                }
        return cov_dict

    def study_level_covariates(self):
        data_list = [d.params for d in self.data.all()]
        
        all_keys = set()

        for d in data_list:
            all_keys |= set(d.keys())

        required_keys = ['GBD Cause', 'Parameter', 'GBD Region', 'Country ISO3 Code',
                         'Sex', 'Year Start', 'Year End', 'Age Start', 'Age End',
                         'Parameter Value', 'Standard Error', 'Units', ]

        redundant_keys = ['_row', 'age_weights', 'id', 'value', 'condition', 'data_type', 'region']

        from dismod3.utils import clean
        from numpy import inf
        additional_keys = sorted(all_keys - set([clean(k) for k in required_keys] + redundant_keys))

        cov_dict = {}
        for k in  additional_keys:
            x_vals = set()
            x_min = 10000.
            x_max = -10000.
            for x in [d.get(k) or 0. for d in data_list]:
                try:
                    x = float(x)
                    x_min = min(x_min, x)
                    x_max = max(x_max, x)
                except ValueError:
                    x_vals.add(x)

            if x_min == 10000. and x_max == -10000.:
                x_min = 0.
                x_max = 0.
            if len(x_vals) <= 1 and x_min == x_max:
                continue

            # for now, only allow numerical covariates
            if x_min == x_max:
                continue
            cov_dict[k] = dict(rate=dict(value=0, default=0),
                               error=dict(value=0, default=0),
                               value=dict(value='0', default='0'),
                               range=[x_min, x_max],
                               category=sorted(x_vals)
                               )

        if len(cov_dict) == 0:
            cov_dict['none'] = {
                'rate': dict(value=0, default=0),
                'error': dict(value=0, default=0),
                'value': dict(value='', default='0.'),  # value must be a string
                'range': [0, 1],
                'category': ['0', '.5', '1']
                }
            
        return cov_dict
