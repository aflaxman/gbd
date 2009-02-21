from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse
from django import forms

import fields
from dismod3.models import AgeSpecificRateFunction

class DiseaseModel(models.Model):
    disease = models.ForeignKey('Disease')
    region = models.ForeignKey('Region')
    sex = fields.SexField()

    notes = models.TextField(blank=True)

    rates = models.ManyToManyField('AgeSpecificRateFunction')
    params_json = models.TextField(blank=True)
    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3'

    def rate_in(self, rate_type):
        filtered_rates = self.rates.filter(rate_type=rate_type)
        if filtered_rates.count() > 0:
            return filtered_rates[0]
        else:
            rf = AgeSpecificRateFunction(disease=self.disease, region=self.region,
                                         sex=self.sex, rate_type=rate_type,
                                         notes='From Disease Model %d' % self.id)
            rf.save()
            self.rates.add(rf)
            return rf
    def i_in(self):
        self.i = self.rate_in('incidence data')
        return self.i
    def p_in(self):
        self.p = self.rate_in('prevalence data')
        return self.p
    def r_in(self):
        self.r = self.rate_in('remission data')
        return self.r
    def f_in(self):
        self.f = self.rate_in('case fatality data')
        return self.f

    def __unicode__(self):
        return '%s; %s; %s' % (self.id, self.disease, self.get_sex_display())
    
    def get_absolute_url(self):
        return reverse("dismod3.views.disease_model_show", args=(self.id,))

    def get_asrf_id_str(self):
        return '_'.join([str(r.id) for r in self.rates.all()])

    def clone(self):
        dm_copy = copy_model_instance(self)
        dm_copy.notes = "Copy of " + self.id
        dm_copy.save()
        for dset in self.rates.all():
            dm_copy.rates.add(dset)
        return dm_copy

    def parse_age_range(self):
        """
        the text field age_range should have the form
        'a0-a1' where a0 and a1 are integers.
        #TODO: include a skip value
        """
        (a0_str, a1_str) = self.age_range.split('-')
        return int(a0_str), int(a1_str)

    def min_age(self):
        return self.parse_age_range()[0]

    def max_age(self):
        return self.parse_age_range()[1]

    def max_rate(self):
        return max([dset.max_rate() for dset in self.datasets.all()])
