from django.db import models
from django.contrib import admin
from django.utils.translation import ugettext as _
from django.core.urlresolvers import reverse
from django import forms
import fields


class DiseaseModel(models.Model):
    name = models.CharField(max_length=200)
    notes = models.TextField(blank=True)
    rates = models.ManyToManyField('AgeSpecificRateFunction')
    params_json = models.TextField(blank=True)
    class Meta:
        # needed to make db work with models directory
        # instead of models.py file
        app_label = 'dismod3'

    def rate_in(self, rate_type):
        filtered_rates = self.rates.filter(rate_type__contains=rate_type, num_rate__ge=1)
        if filtered_rates.count() > 0:
            return filtered_rates[0]
        else:
            return None
    def i_in(self):
        return self.rate_in('incidence')
    def p_in(self):
        return self.rate_in('prevalence')
    def r_in(self):
        return self.rate_in('remission')
    def f_in(self):
        return self.rate_in('case')

    def __unicode__(self):
        return '%s - %s' % (self.name, self.notes)
    
    def get_absolute_url(self):
        return reverse("dismod3.views.disease_model_show", args=(self.id,))

    def get_asrf_id_str(self):
        return '_'.join([str(r.id) for r in self.rates.all()])

    def clone(self):
        dm_copy = copy_model_instance(self)
        dm_copy.name = "Copy of " + self.name
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
