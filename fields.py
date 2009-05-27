from django.db import models
from django.utils.translation import ugettext as _

SEX_CHOICES = [
    ('male', _('Male')),
    ('female', _('Female')),
    ('total', _('Total')),
]

class SexField(models.CharField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('max_length', 10)
        kwargs.setdefault('choices', SEX_CHOICES)

        super(SexField, self).__init__(*args, **kwargs)

    def get_internal_type(self):
        return "CharField"

standardize_sex = {
    'male': 'male',
    'Male': 'male',
    'Male ': 'male',
    'm': 'male',
    'M': 'male',

    'female': 'female',
    'Female': 'female',
    'Female ': 'female',
    'f': 'female',
    'F': 'female',

    'total': 'total',
    'Total': 'total',
    'Total ': 'total',
    't': 'total',
    'T': 'total',
    '': 'total',
    }
    

DATA_TYPE_CHOICES = [
    ('incidence data', _('Incidence')),
    ('prevalence data', _('Prevalence')),
    ('remission data', _('Remission')),
    ('case-fatality data', _('Case-fatality')),
    ('duration data', _('Duration')),
    ('all-cause mortality data', _('All-cause mortality')),
]

class DataTypeField(models.CharField):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('max_length', 100)
        kwargs.setdefault('choices', DATA_TYPE_CHOICES)

        super(DataTypeField, self).__init__(*args, **kwargs)

    def get_internal_type(self):
        return "CharField"

standardize_data_type = {
    'incidence data': 'incidence data',
    'incidence': 'incidence data',
    'Incidence': 'incidence data',
    'Incidence ': 'incidence data',
    'i': 'incidence data',
    'I': 'incidence data',

    'prevalence data': 'prevalence data',
    'prevalence': 'prevalence data',
    'Prevalence': 'prevalence data',
    'Prevalence ': 'prevalence data',
    'p': 'prevalence data',
    'P': 'prevalence data',

    'remission data': 'remission data',
    'remission': 'remission data',
    'Remission': 'remission data',
    'Remission ': 'remission data',
    'r': 'remission data',
    'R': 'remission data',

    'case fatality data': 'case-fatality data',
    'Case Fatality': 'case-fatality data',
    'Case fatality': 'case-fatality data',
    'Case Fatality ': 'case-fatality data',
    'case-fatality': 'case-fatality data',
    'cf': 'case-fatality data',
    'CF': 'case-fatality data',
    'c': 'case-fatality data',
    'C': 'case-fatality data',
    }    
