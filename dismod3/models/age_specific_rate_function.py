from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse
from django import forms

import copy
import pylab as pl
import numpy as np
import simplejson as json

import fields
import django_utils
from dismod3.models import Region, Disease, Rate

default_fit = {'age_mesh': [0.0, 0.5, 3.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]}
def default_fit_json():
    return json.dumps(default_fit)

class AgeSpecificRateFunction(models.Model):
    disease = models.ForeignKey('Disease')
    region = models.ForeignKey('Region')
    rate_type = fields.RateTypeField()
    sex = fields.SexField()
    notes = models.TextField(blank=True)
    rates = models.ManyToManyField('Rate', blank=True)
    num_rates = models.IntegerField(default=0)
    fit_json = models.TextField(default=default_fit_json())

    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3'

    def __init__(self, *args, **kwargs):
        super(AgeSpecificRateFunction, self).__init__(*args, **kwargs)
        try:
            self.fit = json.loads(self.fit_json)
        except ValueError:
            self.fit = copy.copy(default_fit)

    def save(self, force_insert=False, force_update=False):
        # store the fit dict as json text
        self.fit_json = json.dumps(self.fit)

        # denormalize the number of rates, for sorting-by-number
        try:
            self.num_rates = self.rates.count()
        except ValueError: # instance needs to have a primary key before rates.count() will work
            self.num_rates = 0
        
        super(AgeSpecificRateFunction,self).save(force_insert, force_update)

    def clone(self, notes=''):
        new_asrf = django_utils.copy_model_instance(self)
        if notes == '':
            notes = 'Copy of Age Specific Rate Function %s' % self.id
        new_asrf.notes = notes
        new_asrf.fit = copy.copy(default_fit)
        new_asrf.fit['ancestor_ids'] = [self.id] + self.fit.get('ancestor_ids', [])
        new_asrf.save()
        for rate in self.rates.all():
            new_asrf.rates.add(rate)
        return new_asrf
        
        
    def __unicode__(self):
        return '%s; %s; %s; %s' % (self.disease, self.get_rate_type_display(), self.region, self.get_sex_display())

    def get_absolute_url(self):
        return reverse('dismod3.views.age_specific_rate_function_show', args=(self.id,))

    def get_edit_url(self):
        return "/admin/dismod3/agespecificratefunction/%i" % self.id

    def relevant_rates(self):
        """
        construct a list of all Rates with disease, gbd_region, rate_type, and sex
        that match those of self
        """
        return Rate.objects.filter(disease=self.disease, region=self.region,
                                   rate_type=self.rate_type, sex=self.sex)

    def similar_asrfs(self, limit=25):
        """
        return a list of age specific rate functions most similar to this one
        """
        asrfs = AgeSpecificRateFunction.objects
        if self.fit.has_key('ancestor_ids'):
            asrfs = asrfs.filter(id__in=self.fit['ancestor_ids'])
        else:
            asrfs = asrfs.filter(disease=self.disease, rate_type=self.rate_type)

        asrfs = asrfs.exclude(id=self.id)
        asrfs = asrfs.order_by('region')

        return asrfs
    
class ASRFAdmin(admin.ModelAdmin):
    list_display = ('id', 'disease', 'region', 'rate_type',
                    'num_rates', 'sex',)
    list_filter = ['rate_type','disease', 'sex', 'region',]
    search_fields = ['region', 'disease',]

class ASRFForm(forms.ModelForm):
    class Meta:
        model = AgeSpecificRateFunction
        exclude = ('fit_json','rates','num_rates')

def create_multiple(disease, region=None, rate_type='all', sex='all', notes=''):
    params = {'disease': disease, 'notes': notes}

    if region == None:
        regions = Region.objects.all()
    elif not isinstance(region, list):
        regions = [region]
    else:
        regions = region # already a list
    
    if rate_type == 'all':
        rate_types = fields.all_options(fields.RATE_TYPE_CHOICES)[:-1]
    else:
        rate_types = [rate_type]

    if sex == 'all':
        sexes = fields.all_options(fields.SEX_CHOICES)
    else:
        sexes = [sex]

    asrfs = []
    for params['region'] in regions:
        for params['rate_type'] in rate_types:
            for params['sex'] in sexes:
                rf = AgeSpecificRateFunction.objects.create(**params)
                rf.rates = rf.relevant_rates()
                # following two lines should not be necessary---they are handled in model __init__()
                # rf.num_rates = rf.rates.count()
                # rf.fit_json = default_fit_json()
                rf.save()
                asrfs.append(rf)
    

    return asrfs
