from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

import pymc.gp as gp
import numpy as np
import pylab as pl
import csv
from StringIO import StringIO

import gbd.fields
import gbd.view_utils as view_utils
from gbd.unicode_csv_reader import unicode_csv_reader
import dismod3

from models import *
from gbd.dismod3.utils import clean

from dismod3.settings import DISMOD_TWITTER_NAME

class NewDataForm(forms.Form):
    file  = forms.FileField()
    required_data_fields = ['GBD Cause', 'Region', 'Parameter', 'Sex', 'Country',
                            'Age Start', 'Age End', 'Year Start', 'Year End',
                            'Parameter Value', 'Standard Error', 'Units', ]

    tab_separated_values = \
        forms.CharField(required=False,
                        widget=forms.Textarea(attrs={'rows':20, 'cols':80, 'wrap': 'off'}),
                        help_text=_('See <a href="/public/file_formats.html">file format specification</a> for details.'))
    file  = forms.FileField(required=False)
    
    def clean_tab_separated_values(self):
        tab_separated_values = self.cleaned_data['tab_separated_values']
        if not tab_separated_values:
            if not self.files.has_key('file'):
                raise forms.ValidationError(_('TSV field and file field cannot both be blank'))
            else:
                return tab_separated_values
        lines = unicode_csv_reader(StringIO(tab_separated_values), dialect='excel-tab')
        return self.validate(lines)

    def clean_file(self):
        if self.file:
            file_data = self.file.read()
            lines = unicode_csv_reader(StringIO(file_data), dialect='excel-tab')
            return self.validate(lines)

    def validate(self, lines):
        col_names = [clean(col) for col in lines.next()]

        # check that required fields appear
        for field in NewDataForm.required_data_fields:
            if not clean(field) in col_names:
                raise forms.ValidationError(_('Column "%s" is missing') % field)

        data_list = []
        for ii, cells in enumerate(lines):
            # skip blank lines
            if sum([cell == '' for cell in cells]) == len(cells):
                continue
            
            # ensure that something appears for each column
            if len(cells) != len(col_names):
                raise forms.ValidationError(
                    _('Error loading row %d:  found %d fields (expected %d))')
                    % (ii+2, len(cells), len(col_names)))

            # make an associative array from the row data
            data = {}
            for key, val in zip(col_names, cells):
                data[clean(key)] = val.strip()
            data['_row'] = ii+2

            data_list.append(data)

        # ensure that certain cells are the right format
        error_str = _('Row %d:  could not understand entry for %s')
        for r in data_list:
            try:
                r['parameter'] = gbd.fields.standardize_data_type[r['parameter']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter'))
            try:
                r['sex'] = gbd.fields.standardize_sex[r['sex']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Sex'))
            try:
                r['age_start'] = int(r['age_start'])
                # some people think it is a good idea to use 99 as a missing value
                if r['age_start'] == 99:
                    r['age_start'] = 0
                    
                r['age_end'] = int(r['age_end'] or dismod3.MISSING)
                r['year_start'] = int(r['year_start'])
                r['year_end'] = int(r['year_end'])
            except (ValueError, KeyError):
                raise forms.ValidationError(
                    error_str % (r['_row'],
                                 'at least one of Age Start, Age End, Year Start, Year End'))
            try:
                r['parameter_value'] = float(r['parameter_value'])
            except ValueError:
                r['parameter_value'] = dismod3.MISSING

            try:
                r['standard_error'] = float(r['standard_error'])
            except ValueError:
                r['standard_error'] = dismod3.MISSING
                # raise forms.ValidationError(error_str % (r['_row'], 'Standard Error'))
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Standard Error'))
        return data_list


@login_required
def data_upload(request, id=-1):
    # if id != -1, append the data to DiseaseModel.get(id=id)
    if id == -1:
        dm = None
    else:
        dm = get_object_or_404(DiseaseModel, id=id)
        
    if request.method == 'GET':  # no form data is associated with page, yet
        form = NewDataForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = NewDataForm(request.POST, request.FILES)  # A form bound to the POST data

        form.file = request.FILES.get('file')

        if form.is_valid():
            # All validation rules pass, so create new data based on the
            # form contents
            if request.FILES.get('file'):
                data_table = form.cleaned_data['file']
            else:
                data_table = form.cleaned_data['tab_separated_values']

            # make rates from rate_list
            data_list = []
            for d in data_table:
                # add a data point, save it on the data list
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

                # copy mapped data back into d, so that it appears in
                # params
                d.update(args)
                args['params_json'] = json.dumps(d)

                d = Data.objects.create(**args)

                # calculate age weights only when needed, to speed up interface
                # d.calculate_age_weights()
                d.cache_params()
                d.save()
                data_list.append(d)

            # collect this data together into a new model
            args = {}
            args['condition'] = clean(', '.join(set([d.condition for d in data_list])))
            args['sex'] = 'all' #', '.join(set([d.sex for d in data_list]))
            args['region'] = 'global' #'; '.join(set([d.region for d in data_list]))
            args['year'] = '1990-2005' #max_min_str([d.year_start for d in data_list] + [d.year_end for d in data_list])
            if dm:
                dm_json = dm.to_json()
                dm = create_disease_model(dm_json)
            else:
                dm = DiseaseModel.objects.create(**args)
            for d in data_list:
                dm.data.add(d)
            dm.cache_params()
            dm.save()
            return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_summary', args=[dm.id])) # Redirect after POST

    return render_to_response('data_upload.html', {'form': form, 'dm': dm})


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
def dismod_list(request, format='html'):
    dm_filter = DiseaseModel.objects.all().order_by('-id')
    if format == 'html':
        return render_to_response('dismod_list.html',
                                  {'paginated_models': view_utils.paginated_models(request, dm_filter)})
    else:
        raise Http404

@login_required
def dismod_show(request, id, format='html'):
    if isinstance(id, DiseaseModel):
        dm = id
    else:
        dm = get_object_or_404(DiseaseModel, id=id)

    if format == 'html':
        dm.px_hash = dismod3.sparkplot_boxes(dm.to_json())
        return render_to_response('dismod_show.html', {'dm': dm})
    elif format == 'json':
        return HttpResponse(dm.to_json(), view_utils.MIMETYPE[format])
    elif format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.tile_plot_disease_model(dm.to_json(),
                                        dismod3.utils.gbd_keys(type_list=dismod3.utils.output_data_types))
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    else:
        raise Http404

@login_required
def dismod_show_by_region_year_sex(request, id, region, year, sex, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)

    if format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.tile_plot_disease_model(dm.to_json(),
                                        dismod3.utils.gbd_keys(
                type_list=dismod3.utils.output_data_types,
                region_list=[region],
                year_list=[year],
                sex_list=[sex]))
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    else:
        raise Http404

@login_required
def dismod_show_by_region(request, id, region, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)

    if format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.tile_plot_disease_model(dm.to_json(),
                                        dismod3.utils.gbd_keys(
                type_list=dismod3.utils.output_data_types,
                region_list=[region]))
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    else:
        raise Http404

    
@login_required
def dismod_find_and_show(request, condition, format='html'):
    try:
        dm = DiseaseModel.objects.filter(condition=condition).latest('id')
    except DiseaseModel.DoesNotExist:
        raise Http404
    return dismod_show(request, dm, format)

@login_required
def dismod_sparkplot(request, id, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)
    if format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.sparkplot_disease_model(dm.to_json())
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    else:
        raise Http404

@login_required
def dismod_overlay_plot(request, id, condition, type, region, year, sex, format='png'):
    if not format in ['png', 'svg', 'eps', 'pdf']:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dismod3.overlay_plot_disease_model(dm.to_json(), keys)
    pl.title('%s; %s; %s; %s' % (dismod3.plotting.prettify(condition),
                                 dismod3.plotting.prettify(region), year, sex))
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])


@login_required
def dismod_tile_plot(request, id, condition, type, region, year, sex, format='png'):
    if not format in ['png', 'svg', 'eps', 'pdf']:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
    dismod3.tile_plot_disease_model(dm.to_json(), keys)
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])


@login_required
def dismod_summary(request, id, format='html'):
    if not format in ['html']:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    data = dm.data.all()
    data_counts = []
    for r in dismod3.gbd_regions:
        c = {}

        c['region'] = r
        c['clean_region'] = clean(r)
        
        for type, data_type in [['i', 'incidence data'],
                                ['p', 'prevalence data'],
                                ['r', 'remission data'],
                                ['cf', 'case-fatality data']]:
            c[type] = \
                len([d for d in data if d.relevant_to(data_type, r, year='all', sex='all')])

        # also count relative-risk, mortality, and smr data as case-fatality data
        type = 'cf'
        for data_type in ['relative-risk data', 'smr data', 'mortality data']:
            c[type] += \
                    len([d for d in data if clean(d.data_type) == clean(data_type)
                         and clean(d.gbd_region) == clean(r)])
        
        data_counts.append(c)
    data_counts = sorted(data_counts, reverse=True,
                         key=lambda c: c['i'] + c['p'] + c['r'] + c['cf'])
    total = {}
    for type in ['i', 'p', 'r', 'cf']:
        total[type] = sum([d[type] for d in data_counts])
        
    if format == 'html':
        dm.px_hash = dismod3.sparkplot_boxes(dm.to_json())
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total})
    else:
        raise Http404

    
class NewDiseaseModelForm(forms.Form):
    model_json = \
        forms.CharField(required=True,
                        widget=forms.Textarea(attrs={'rows':20, 'cols':80, 'wrap': 'off'}),
                        help_text=_('See <a href="/public/dismod_data_json.html">dismod json specification</a> for details.'))
    def clean_model_json(self):
        model_json = self.cleaned_data['model_json']
        try:
            model_dict = json.loads(model_json)
        except ValueError:
            raise forms.ValidationError('JSON object could not be decoded')
        if not model_dict.get('params'):
            raise forms.ValidationError('missing params')
        for key in ['condition', 'sex', 'region', 'year']:
            if not model_dict['params'].get(key):
                raise forms.ValidationError('missing params.%s' % key)

        # store the model dict for future use
        self.cleaned_data['model_dict'] = model_dict
        return model_json

@login_required
def dismod_upload(request):
    if request.method == 'GET':  # no form data is associated with page, yet
        form = NewDiseaseModelForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = NewDiseaseModelForm(request.POST)  # A form bound to the POST data

        if form.is_valid():
            # All validation rules pass, so update or create new disease model
            model_dict = form.cleaned_data['model_dict']
            id = model_dict['params'].get('id', -1)
            if id > 0:
                dm = get_object_or_404(DiseaseModel, id=id)
                for key,val in model_dict['params'].items():
                    if type(val) == dict and dm.params.has_key(key):
                        dm.params[key].update(val)
                    else:
                        dm.params[key] = val
                dm.cache_params()
                dm.save()
            else:
                dm = create_disease_model(form.cleaned_data['model_json'])

            return HttpResponseRedirect(dm.get_absolute_url()) # Redirect after POST

    return render_to_response('dismod_upload.html', {'form': form})

@login_required
def job_queue_list(request):
    # accept format specified in url
    format = request.GET.get('format', 'html')

    dm_list = DiseaseModel.objects.filter(needs_to_run=True)
    if format == 'json':
        return HttpResponse(json.dumps([ dm.id for dm in dm_list ]),
                            view_utils.MIMETYPE[format])
    else:
        # more formats shall be added one day
        raise Http404
        
class JobRemovalForm(forms.Form):
    id = forms.IntegerField()
    
@login_required
def job_queue_remove(request):
    if request.method == 'GET':  # no form data is associated with page, yet
        form = JobRemovalForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = JobRemovalForm(request.POST)  # A form bound to the POST data

        if form.is_valid():
            dm = get_object_or_404(DiseaseModel, id=form.cleaned_data['id'])
            if dm.needs_to_run:
                dm.needs_to_run = False
                dm.save()

            return HttpResponseRedirect(
                reverse('gbd.dismod_data_server.views.job_queue_list') + '?format=json')
    return render_to_response('job_queue_remove.html', {'form': form})

@login_required
def job_queue_add(request, id):
    # only react to POST requests
    if request.method != 'POST':
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)
    dm.needs_to_run = True
    if request.POST.has_key('estimate_type'):
        dm.params['estimate_type'] = request.POST['estimate_type']
    dm.cache_params()
    dm.save()

    return HttpResponseRedirect('http://twitter.com/' + DISMOD_TWITTER_NAME)

@login_required
def dismod_run(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    return render_to_response('dismod_run.html', {'dm': dm})

@login_required
def dismod_update_covariates(request, id):
    # this is slow, so it has been spun off into a separate view that
    # can be run only when necessary

    # only react to POST requests, since this may changes database data
    if request.method != 'POST':
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)
    for d in dm.data.all():
        d.age_weights()  # will cache value if it is not already cached
    
    return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_run', args=[dm.id])) # Redirect after POST

@login_required
def dismod_adjust(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)

    if request.method == 'GET':
        return render_to_response('dismod_adjust.html', {'dm': dm, 'sessionid': request.COOKIES['sessionid']})
    elif request.method == 'POST':
        dm.params['global_priors_json'] = request.POST['JSON']
        dm.cache_params()

        dj = dismod3.disease_json.DiseaseJson(dm.to_json())
        dj.extract_params_from_global_priors()
        new_dm = create_disease_model(dj.to_json())

        return HttpResponse(reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))

@login_required
def dismod_preview_priors(request, id, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)
    dm = dismod3.disease_json.DiseaseJson(dm.to_json())

    if request.method == 'POST':
        dm.params['global_priors_json'] = request.POST['JSON']
        dm.extract_params_from_global_priors()

    if format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.plot_prior_preview(dm)
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    else:
        raise Http404

def my_prior_str(dict, smooth_key, conf_key, zero_before_key, zero_after_key):
    s = ''
    if dict.get(smooth_key) and dict[smooth_key] != '(none)':
        s += 'smooth %s, ' % dict[smooth_key]
    if dict.get(conf_key) and dict[conf_key] != '(none)':
        s += 'confidence %s, ' % dict[conf_key]
    if dict.get(zero_before_key):
        s += 'zero 0 %s, ' % dict[zero_before_key]
    if dict.get(zero_after_key):
        s += 'zero %s %d, ' % (dict[zero_after_key], dismod3.utils.MAX_AGE)

    return s
