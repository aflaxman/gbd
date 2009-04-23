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

def clean(str):
    return str.lower().replace(' ', '')

class DataCreationForm(forms.Form):
    tab_separated_values = forms.CharField(required=True, widget=forms.widgets.Textarea(), help_text='See "data input spec.doc" for details')

    def clean_tab_separated_values(self):
        # TODO: make a special form derived from CharField that does this split with the csv package
        tab_separated_values = self.cleaned_data['tab_separated_values']
        lines = tab_separated_values.split('\n')

        col_names = [clean(col) for col in lines.pop(0).split('\t')]

        # check that required fields appear
        for field in ['GBD Cause', 'Region', 'Parameter', 'Sex', 'Country',
                      'Age Start', 'Age End', 'Estimate Year Start', 'Estimate Year End',
                      'Parameter Value', 'Lower Value', 'Upper Value', 'Units', 'Type of Bounds', ]:
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
        #import pdb; pdb.set_trace()
        for r in data_list:
            r['parameter'] = fields.standardize_data_type[r['parameter']]
            r['sex'] = fields.standardize_sex[r['sex']]
            r['agestart'] = int(r['agestart'])
            r['ageend'] = int(r['ageend'] or fields.MISSING)
            r['estimateyearstart'] = int(r['estimateyearstart'])
            r['estimateyearend'] = int(r['estimateyearend'])
            r['parametervalue'] = float(r['parametervalue'])
            r['lowervalue'] = float(r['lowervalue'] or fields.MISSING)
            r['uppervalue'] = float(r['uppervalue'] or fields.MISSING)
            r['units'] = float((re.findall('([\d\.]+)', r['units']) or [fields.MISSING])[0])
            r['typeofbounds'] = float((re.findall('([\d\.]+)', r['typeofbounds']) or [fields.MISSING])[0])
            # TODO: catch ValueError and KeyError, and raise informative error instead, forms.ValidationError('useful msg here')
            # and write tests for this, too
        return data_list


@login_required
def data_new(request):
    if request.method == 'POST': # If the form has been submitted...
        form = DataCreationForm(request.POST) # A form bound to the POST data
        import pdb; pdb.set_trace()
        if form.is_valid(): # All validation rules pass
            data_table = form.cleaned_data['tab_separated_values']

            # make rates from rate_list
            data_list = []
            for d in data_table:
                # add a rate, save it on a list
                args = {}
                args['condition'] = d['gbdcause']
                args['gbd_region'] = d['region']
                args['region'] = d['country']
                args['data_type'] = d['parameter']
                args['sex'] = d['sex']
                args['age_end'] = d['ageend']
                args['age_start'] = d['agestart']
                args['time_start'] = d['estimateyearstart']
                args['time_end'] = d['estimateyearend']

                # TODO: deal with the standard error correctly
                if d['units'] == fields.MISSING:
                    d['units'] = 1.0
                args['value'] = d['parametervalue']
                args['standard_error'] = .01
                
                args['params_json'] = json.dumps(d)

                d = Data(**args)
                d.save()
                data_list.append(d)
                
            # TODO: collect information on new data, use it to display
            # all related data together
            return HttpResponseRedirect(d.get_edit_url()) # Redirect after POST
    else:
        form = DataCreationForm()

    return render_to_response('data_new.html', {'form': form})


@login_required
def data_show(request, id):
    data = get_object_or_404(Data, pk=id)
    data.view_list = [[_('Condition'), data.condition],
                      [_('Data Type'), data.get_data_type_display()],
                      [_('Sex'), data.get_sex_display()],
                      [_('GBD Region'), data.gbd_region],
                      [_('Region'), data.region],
                      [_('Age Range'), data.age_str()],
                      [_('Time Range'), data.time_str()],
                      [_('Value'), data.value_str()],
                      ]
    return render_to_response('data_show.html', view_utils.template_params(data))

