from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

import pymc.gp as gp
import numpy as np
import pylab as pl

from dismod3.models import fields
from dismod3.views import view_utils

from models import *

def max_min_str(num_list):
    a = min(num_list)
    b = max(num_list)
    if a == b:
        return '%d' % a
    else:
        return '%d-%d' % (a,b)

def clean(str):
    return str.strip().lower().replace(' ', '_')

class NewDataForm(forms.Form):
    required_data_fields = ['GBD Cause', 'Region', 'Parameter', 'Sex', 'Country',
                            'Age Start', 'Age End', 'Year Start', 'Year End',
                            'Parameter Value', 'Standard Error', 'Units', ]

    tab_separated_values = forms.CharField(required=True, widget=forms.widgets.Textarea(), help_text='See "data input spec" for details')

    def clean_tab_separated_values(self):
        # TODO: make a special form derived from CharField that does this split with the csv package
        tab_separated_values = self.cleaned_data['tab_separated_values']
        lines = tab_separated_values.split('\n')

        col_names = [clean(col) for col in lines.pop(0).split('\t')]
        #import pdb; pdb.set_trace()

        # check that required fields appear
        for field in NewDataForm.required_data_fields:
            if not clean(field) in col_names:
                raise forms.ValidationError('Column "%s" is missing' % field)

        data_list = []
        for ii, line in enumerate(lines):
            # skip blank lines
            if line == '':
                continue

            cells = line.split('\t')
            
            # ensure that something appears for each column
            if len(cells) != len(col_names):
                raise forms.ValidationError('Error loading row %d:  missing fields detected' % (ii+2))

            # make an associative array from the row data
            data = {}
            for key, val in zip(col_names, cells):
                data[clean(key)] = val.strip()

            data_list.append(data)

        # ensure that certain cells are the right format
        for r in data_list:
            r['parameter'] = fields.standardize_data_type[r['parameter']]
            r['sex'] = fields.standardize_sex[r['sex']]
            r['age_start'] = int(r['age_start'])
            r['age_end'] = int(r['age_end'] or fields.MISSING)
            r['year_start'] = int(r['year_start'])
            r['year_end'] = int(r['year_end'])
            r['parameter_value'] = float(r['parameter_value'])
            r['standard_error'] = float(r['standard_error'])
            # TODO: catch ValueError and KeyError, and raise informative error instead, forms.ValidationError('useful msg here')
            # and write tests for this, too
        return data_list


@login_required
def data_new(request):
    if request.method == 'GET':     # no form data is associated with page, yet
        form = NewDataForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = NewDataForm(request.POST) # A form bound to the POST data

        if form.is_valid(): # All validation rules pass, create new data based on the form contents
            data_table = form.cleaned_data['tab_separated_values']

            # make rates from rate_list
            data_list = []
            for d in data_table:
                # add a rate, save it on a list
                args = {}
                args['condition'] = d['gbd_cause']
                args['gbd_region'] = d['region']
                args['region'] = d['country']
                args['data_type'] = d['parameter']
                args['sex'] = d['sex']
                args['age_start'] = d['age_start']
                args['age_end'] = d['age_end']
                args['year_start'] = d['year_start']
                args['year_end'] = d['year_end']

                args['value'] = d['parameter_value']
                args['standard_error'] = d['standard_error']
                
                args['defaults'] = {'params_json': json.dumps(d)}

                d, is_new = Data.objects.get_or_create(**args)
                d.calculate_age_weights()
                data_list.append(d)
                
            # collect this data together into a new model
            args = {}
            args['condition'] = ', '.join(set([d.condition for d in data_list]))
            args['sex'] = ', '.join(set([d.sex for d in data_list]))
            args['region'] = '; '.join(set([d.region for d in data_list]))
            args['year'] = max_min_str([d.year_start for d in data_list] + [d.year_end for d in data_list])

            #import pdb; pdb.set_trace()
            dm = DiseaseModel.objects.create(**args)
            for d in data_list:
                dm.data.add(d)
            dm.cache_params()
            dm.save()
            
            return HttpResponseRedirect(dm.get_absolute_url()) # Redirect after POST

    return render_to_response('data_new.html', {'form': form})


@login_required
def data_show(request, id):
    data = get_object_or_404(Data, pk=id)
    data.view_list = [[_('Condition'), data.condition],
                      [_('Data Type'), data.data_type],
                      [_('Sex'), data.get_sex_display()],
                      [_('GBD Region'), data.gbd_region],
                      [_('Region'), data.region],
                      [_('Age'), data.age_str()],
                      [_('Year'), data.year_str()],
                      [_('Value'), data.value_str()],
                      ]
    return render_to_response('data_show.html', view_utils.template_params(data))


@login_required
def disease_model_show(request, id, format='html'):
    dm = get_object_or_404(DiseaseModel, id=id)

    if format == 'html':
        return render_to_response('disease_model_show.html', view_utils.template_params(dm))
    elif format == 'json':
        dm.params.update(condition=dm.condition, sex=dm.sex, region=dm.region, year=dm.year)
        data_str = json.dumps({'params': dm.params, 'data': [[d.id, d.params] for d in dm.data.all()]},
                              sort_keys=True, indent=2)
        return HttpResponse(data_str, view_utils.MIMETYPE[format])
    else:
        raise Http404


class NewDiseaseModelForm(forms.Form):
    model_json = forms.CharField(required=True, widget=forms.widgets.Textarea(), help_text='See source for details')
    def clean_model_json(self):
        # TODO: make a special form derived from CharField that does this split with the csv package
        #import pdb; pdb.set_trace()
        model_json = self.cleaned_data['model_json']
        model_dict = json.loads(model_json)
        if not model_dict.get('params'):
            raise forms.ValidationError('missing params')
        if not model_dict['params'].get('condition'):
            raise forms.ValidationError('missing params.condition')
        if not model_dict['params'].get('sex'):
            raise forms.ValidationError('missing params.sex')
        if not model_dict['params'].get('region'):
            raise forms.ValidationError('missing params.region')
        if not model_dict['params'].get('year'):
            raise forms.ValidationError('missing params.year')
        return model_dict

@login_required
def disease_model_new(request):
    if request.method == 'GET':     # no form data is associated with page, yet
        form = NewDiseaseModelForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = NewDiseaseModelForm(request.POST) # A form bound to the POST data

        if form.is_valid(): # All validation rules pass, create new data based on the form contents
            model_dict = form.cleaned_data['model_json']
            params = model_dict['params']

            args = {}
            args['condition'] = params['condition']
            args['sex'] = params['sex']
            args['region'] = params['region']
            args['year'] = params['year']

            dm = DiseaseModel.objects.create(**args)
            for d_id, d_data in model_dict.get('data', []):
                dm.data.add(d_id)

            dm.params = params
            dm.cache_params()
            dm.save()
            
            return HttpResponseRedirect(dm.get_absolute_url()) # Redirect after POST

    return render_to_response('disease_model_new.html', {'form': form})
