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
    'P': 'total',
    '': 'total',
    }
    

DATA_TYPE_CHOICES = [
    ('incidence data', _('Incidence Rate')),
    ('prevalence data', _('Prevalence Rate')),
    ('remission data', _('Remission Rate')),
    ('excess-mortality data', _('Excess Mortality Rate')),
    ('mrr data', _('Mortality Rate Ratio')),
    ('smr data', _('Standardized Mortality Ratio')),
    ('mortality data', _('With-condition Mortality Rate')),
    ('duration data', _('Case Duration')),
    ('all-cause mortality data', _('All-cause Mortality Rate')),
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
    'remission rate': 'remission data',
    'Remission Rate': 'remission data',
    'r': 'remission data',
    'R': 'remission data',

    'Excess Mortality': 'excess-mortality data',
    'excess mortality data': 'excess-mortality data',
    'case fatality data': 'excess-mortality data',
    'Case Fatality': 'excess-mortality data',
    'Case fatality': 'excess-mortality data',
    'Case Fatality ': 'excess-mortality data',
    'case-fatality': 'excess-mortality data',

    'Relative Mortality Ratio': 'relative-risk data',
    'mrr data': 'relative-risk data',
    'relative-risk': 'relative-risk data',
    'relative risk': 'relative-risk data',
    'RR': 'relative-risk data',
    'Mortality: RR': 'relative-risk data',
    'relative-risk data': 'relative-risk data',
    'relative risk data': 'relative-risk data',

    'Standardized Mortality Ratio': 'smr data',
    'SMR': 'smr data',
    'Specific Mortality Risk': 'smr data',
    'smr': 'smr data',

    'mortality data': 'mortality data',
    'Mortality': 'mortality data',
    'mortality': 'mortality data',
    'crude mortality': 'mortality data',
    'with-condition mortality rate': 'mortality data',
    'with-condition mortality': 'mortality data',
    'With-Condition Mortality': 'mortality data',
    'With Condition Mortality': 'mortality data',
    
    'Duration': 'duration data',
    'duration': 'duration data',
    'duration data': 'duration data',
    }    
