from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse
from datetime import datetime
import numpy as np
import gbd.fields
from gbd.dismod3.utils import clean
from gbd.dismod3.neg_binom_model import countries_for

class CovariateType(models.Model):
    slug = models.CharField(max_length=50)
    uploader = models.CharField(max_length=30)
    upload_time = models.DateTimeField(default=datetime.now)
    source = models.TextField()
    last_modified_time = models.DateTimeField(default=datetime.now)
    description = models.TextField()
    year_start = models.IntegerField()
    year_end = models.IntegerField()
    region_only = models.BooleanField()

    def __unicode__(self):
        return self.slug

    def get_absolute_url(self):
        return reverse('gbd.covariate_data_server.views.covariate_type_show', args=[self.id])

class CovariateTypeAdmin(admin.ModelAdmin):
    list_display = ('id', 'slug', 'uploader', 'upload_time', 'source', 'last_modified_time', 'description')

class CovariateAdmin(admin.ModelAdmin):
    list_display  = ('id', 'type', 'iso3', 'region', 'sex', 'age', 'year', 'value')
    list_filter   = ['sex', 'type',]
    search_fields = ['type', 'iso3', 'year',]

class Covariate(models.Model):
    """ Model for Covariate Data
    """
    type = models.ForeignKey(CovariateType)
    iso3 = models.CharField(max_length=3)
    year = models.IntegerField()
    sex = gbd.fields.SexField()
    country_year = models.CharField(max_length=8)
    value = models.FloatField()
    age = models.CharField(default='all',max_length=3)
    region = models.CharField(default='',max_length=28)

    def save(self, *args, **kwargs):
        self.country_year = '%s-%d' % (self.iso3, self.year)
        super(Covariate, self).save(*args, **kwargs)

    def __unicode__(self):
        return '%s: %s, %s, %s' % (self.type, self.iso3, self.year, self.get_sex_display(),)

    def get_absolute_url(self):
        return reverse('gbd.covariate_data_server.views.covariate_show', args=[self.type, self.iso3, self.sex, 'png'])
    

    def to_dict(self):
        return dict(type=self.type, iso3=self.iso3, year=self.year, sex=self.sex, value=self.value)







