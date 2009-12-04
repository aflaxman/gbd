from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse

import gbd.fields

class CovariateType(models.Model):
    slug = models.CharField(max_length=50)
    description = models.TextField()

    def __unicode__(self):
        return self.slug


class CovariateAdmin(admin.ModelAdmin):
    list_display  = ('id', 'type', 'iso3', 'sex', 'year', 'value')
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

    def save(self, *args, **kwargs):
        self.country_year = '%s-%d' % (self.iso3, self.year)
        super(Covariate, self).save(*args, **kwargs)

    def __unicode__(self):
        return '%s: %s, %s, %s' % (self.type, self.iso3, self.year, self.get_sex_display(),)

    def get_absolute_url(self):
        return reverse('gbd.covariate_data_server.views.covariate_show', args=(self.id,))

    def to_dict(self):
        return dict(type=self.type, iso3=self.iso3, year=self.year, sex=self.sex, value=self.value)
