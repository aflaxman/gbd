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
import time
import os
import socket
from shutil import rmtree
import math

import gbd.fields
import gbd.view_utils as view_utils
from gbd.unicode_csv_reader import unicode_csv_reader
import dismod3

from models import *
from gbd.dismod3.utils import clean
from gbd.dismod3.settings import JOB_LOG_DIR, JOB_WORKING_DIR, SERVER_LOAD_STATUS_HOST, SERVER_LOAD_STATUS_PORT, SERVER_LOAD_STATUS_SIZE
from gbd.dismod3.table import population_by_region_year_sex
from gbd.dismod3.neg_binom_model import countries_for
import fcntl

class NewDataForm(forms.Form):
    file  = forms.FileField()
    required_data_fields = ['GBD Cause', 'Region', 'Parameter', 'Sex', 'Country ISO3 Code',
                            'Age Start', 'Age End', 'Year Start', 'Year End',
                            'Parameter Value', 'Units', ]

    tab_separated_values = \
        forms.CharField(required=False,
                        widget=forms.Textarea(attrs={'rows':20, 'cols':110, 'wrap': 'off'}),
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
        """
Required data fields:
--------------------------------------------------------------------------------
Name                               Type    Limit
--------------------------------------------------------------------------------
GBD Cause                          str     one of the GBD causes
Region                             str     one of the GBD regions
Parameter                          str     standardize_data_type
Sex                                str     standardize_sex
Country ISO3 Code                  str     an ISO3 code in the region
Age Start                          int     [0, 100], <= Age End
Age End                            int     [0, 100], >= Age Start
Year Start                         int     [1950, 2010], <= Year End
Year End                           int     [1950, 2010], >= Year Start
Parameter Value                    float   >= 0
Units                              float   >= 1

Recommended data fields:
--------------------------------------------------------------------------------
Name                               Type    Limit
--------------------------------------------------------------------------------
Study ID                           empty or int     >= 0
Sequela                            empty or str     one of the GBD sequela codes
Case Definition                    empty or str     none
Coverage                           empty or float   [0,1]
Study Size N For This Year & Sex   empty or int     > 0, <= Total Study Size N
Lower CI                           empty or float   > 0 <= Parameter Value
Upper CI                           empty or float   >= Parameter Value
Standard Error                     empty or float   > 0
Total Study Size N                 empty or int     > 0
Design Factor                      empty or float   >= 1
Citation                           empty or str     none
Urbanicity                         empty or float   [0, 1]
Ignore                             empty or int     [0, 1]

Optional data fields:
No checks
        """
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
        gbd_cause = ''

        for r in data_list:
            # check required data fields
            try:
                r['gbd_cause'] = str(r['gbd_cause'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'GBD Cause'))
            if gbd_cause == '':
                gbd_cause = r['gbd_cause']
            else:
                if gbd_cause != r['gbd_cause']:
                    raise forms.ValidationError(error_str % (r['_row'], 'GBD Cause (all GBD Causes must be the same)'))

            try:
                r['region'] = str(r['region'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Region'))
            if not clean(r['region']) in [clean(region) for region in dismod3.gbd_regions]:
                raise forms.ValidationError(error_str % (r['_row'], 'Region'))

            try:
                r['parameter'] = gbd.fields.standardize_data_type[r['parameter']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter'))

            try:
                r['sex'] = gbd.fields.standardize_sex[r['sex']]
            except KeyError:
                raise forms.ValidationError(error_str % (r['_row'], 'Sex'))

            try:
                r['country_iso3_code'] = str(r['country_iso3_code'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Country ISO3 Code'))
            if not r['country_iso3_code'] in countries_for[clean(r['region'])]:
                raise forms.ValidationError(error_str % (r['_row'], 'Country ISO3 Code (%s is not in %s)' % (r['country_iso3_code'], r['region'])))

            try:
                r['age_start'] = int(r['age_start'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start'))
            if r['age_start'] < 0 or r['age_start'] > 100:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start (must be in range [0, 100])'))

            try:
                r['age_end'] = int(r['age_end'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Age End'))
            if r['age_end'] < 0 or r['age_end'] > 100:
                raise forms.ValidationError(error_str % (r['_row'], 'Age End (must be in range [0, 100])'))

            if r['age_start'] > r['age_end']:
                raise forms.ValidationError(error_str % (r['_row'], 'Age Start (must be greater than Age End)'))

            try:
                r['year_start'] = int(r['year_start'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start'))
            if r['year_start'] < 1950 or r['year_start'] > 2010:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start (must be in range [1950, 2010])'))

            try:
                r['year_end'] = int(r['year_end'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Year End'))
            if r['year_end'] < 1950 or r['year_end'] > 2010:
                raise forms.ValidationError(error_str % (r['_row'], 'Year End (must be in range [1950, 2010])'))
   
            if r['year_start'] > r['year_end']:
                raise forms.ValidationError(error_str % (r['_row'], 'Year Start (must be greater than Year End)'))

            try:
                r['parameter_value'] = float(r['parameter_value'])
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value'))
            if r['parameter_value'] < 0:
                raise forms.ValidationError(error_str % (r['_row'], 'Parameter Value (must be greater than 0)'))

            units = 0
            try:
                units = float(r['units'].replace(',', '').replace('per ', ''))
            except ValueError:
                raise forms.ValidationError(error_str % (r['_row'], 'Units'))
            if units < 1:
                raise forms.ValidationError(error_str % (r['_row'], 'Units (must be greater than 1)'))

            # check recommended data fields
            if 'study_id' in col_names and r['study_id'] != '':
                try:
                    r['study_id'] = int(r['study_id'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study ID'))
                if r['study_id'] < 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study ID (must be greater than 0)'))

            if 'sequela' in col_names and r['sequela'] != '':
                try:
                    r['sequela'] = str(r['sequela'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Sequela'))

            if 'case_definition' in col_names and r['case_definition'] != '':
                try:
                    r['case_definition'] = str(r['case_definition'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Case Definition'))

            if 'coverage' in col_names and r['coverage'] != '':
                try:
                    r['coverage'] = float(r['coverage'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Coverage'))
                if r['coverage'] < 0 or r['coverage'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Coverage (must be in range [0, 1])'))

            if 'study_size_n_for_this_year_&_sex' in col_names and r['study_size_n_for_this_year_&_sex'] != '':
                try:
                    r['study_size_n_for_this_year_&_sex'] = int(r['study_size_n_for_this_year_&_sex'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study Size N For This Year and Sex'))
                if r['study_size_n_for_this_year_&_sex'] <= 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study Size N For This Year and Sex (must be greater than 0)'))

            if 'lower_ci' in col_names and r['lower_ci'] != '':
                try:
                    r['lower_ci'] = float(r['lower_ci'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Lower CI'))
                if r['lower_ci'] <= 0 or r['lower_ci'] > r['parameter_value']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Lower CI (must be less than parameter value)'))

            if 'upper_ci' in col_names and r['upper_ci'] != '':
                try:
                    r['upper_ci'] = float(r['upper_ci'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Upper CI'))
                if r['upper_ci'] < r['parameter_value']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Upper CI (must be greater than Parameter Value)'))

            if 'standard_error' in col_names:
                if r['standard_error'] != '':
                    try:
                        r['standard_error'] = float(r['standard_error'])
                    except ValueError:
                        raise forms.ValidationError(error_str % (r['_row'], 'Standard Error'))
                    if r['standard_error'] <= 0 and r['standard_error'] != -99:
                        raise forms.ValidationError(error_str % (r['_row'], 'Standard Error (must be greater than 0 or -99 for missing)'))

            if 'total_study_size_n' in col_names and r['total_study_size_n'] != '':
                try:
                    r['total_study_size_n'] = int(r['total_study_size_n'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Total Study Size N'))
                if r['total_study_size_n'] <= 0:
                    raise forms.ValidationError(error_str % (r['_row'], 'Total Study Size N (must be greater than 0)'))

            if 'total_study_size_n' in col_names and 'study_size_n_for_this_year_&_sex' in col_names and r['study_size_n_for_this_year_&_sex'] != '' and r['total_study_size_n'] != '':
                if r['study_size_n_for_this_year_&_sex'] > r['total_study_size_n']:
                    raise forms.ValidationError(error_str % (r['_row'], 'Study Size N For This Year and Sex (must be at most Total Study Size N)'))

            if 'design_factor' in col_names and r['design_factor'] != '':
                try:
                    r['design_factor'] = float(r['design_factor'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Design Factor'))
                if r['design_factor'] < 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Design Factor (must be greater than 1)'))

            #if 'citation' in col_names and r['citation'] != '':
            #    try:
            #        r['citation'] = str(r['citation'])
            #    except ValueError:
            #        raise forms.ValidationError(error_str % (r['_row'], 'Citation'))

            if 'urbanicity' in col_names and r['urbanicity'] != '':
                try:
                    r['urbanicity'] = float(r['urbanicity'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Urbanicity'))
                if r['urbanicity'] < 0 or r['urbanicity'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Urbanicity (must be in range [0, 1])'))

            if 'ignore' in col_names and r['ignore'] != '':
                try:
                    r['ignore'] = int(r['ignore'])
                except ValueError:
                    raise forms.ValidationError(error_str % (r['_row'], 'Ignore'))
                if r['ignore'] < 0 or r['ignore'] > 1:
                    raise forms.ValidationError(error_str % (r['_row'], 'Ignore (must be 0 or 1)'))

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
                args['region'] = d['country_iso3_code']
                args['data_type'] = d['parameter']
                args['sex'] = d['sex']
                args['age_start'] = d['age_start']
                args['age_end'] = d['age_end']
                args['year_start'] = d['year_start']
                args['year_end'] = d['year_end']

                args['value'] = d['parameter_value']
                try:
                    args['standard_error'] = d['standard_error']
                except KeyError:
                    args['standard_error'] = dismod3.MISSING

                if args['standard_error'] == '':
                    args['standard_error'] = dismod3.MISSING

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
            args['creator'] = request.user
            if dm:
                dm_json = dm.to_json()
                dm = create_disease_model(dm_json, request.user)
            else:
                dm = DiseaseModel.objects.create(**args)
            for d in data_list:
                dm.data.add(d)
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
    if request.GET.get('show') == 'all':
        dm_filter = DiseaseModel.objects.all().order_by('-id')
    else:
        dm_filter = DiseaseModel.objects.filter(creator=request.user).order_by('-id')
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
        dm.px_hash = dismod3.sparkplot_boxes(dm.to_json({'key': 'none'}))
        
        return render_to_response('dismod_show.html',
                                  {'dm': dm,
                                  'paginated_models': view_utils.paginated_models(request, dm.data.all())})
    elif format == 'json':
        return HttpResponse(dm.to_json(), view_utils.MIMETYPE[format])
    elif format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.tile_plot_disease_model(dm.to_json(),
                                        dismod3.utils.gbd_keys(type_list=dismod3.utils.output_data_types))
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    elif format == 'xls':
        group_size = int(request.GET.get('group_size', 1))

        content = dismod3.table(dm.to_json(),
                      dismod3.utils.gbd_keys(
                type_list=dismod3.utils.output_data_types), request.user, group_size)
        return HttpResponse(content, mimetype='application/ms-excel')
    elif format == 'dta':
        import subprocess, csv

        # TODO: pick an appropriate temp file name, so that there are not collisions
        from django.core.files.base import File, ContentFile
        from django.core.files.storage import default_storage
        p_csv, is_new = dm.params.get_or_create(key='csv-plot')
        if not p_csv.file:
            fname = default_storage.save('%s/%s_%s.csv' %(dm.id, dm.id, dm.condition), ContentFile(''))
            p_csv.file = File(default_storage.open(fname))
            p_csv.save()
            dm.params.add(p_csv)

        p_dta, is_new = dm.params.get_or_create(key='dta-plot')
        if not p_dta.file:
            fname = default_storage.save('%s/%s_%s.dta' %(dm.id, dm.id, dm.condition), ContentFile(''))
            p_dta.file = File(default_storage.open(fname))
            p_dta.save()
            dm.params.add(p_dta)

        if is_new:
            X = ['type, region, sex, year, age, prior, posterior, upper, lower'.split(', ')]
            dm = dismod3.disease_json.DiseaseJson(dm.to_json())
            for t in dismod3.utils.output_data_types:
                for r in dismod3.settings.gbd_regions:
                    r = clean(r)
                    for s in ['male', 'female']:
                        for y in [1990, 2005]:
                            k = dismod3.utils.gbd_key_for(t, r, y, s)

                            prior = dm.get_mcmc('emp_prior_mean', k)
                            if len(prior) == 0:
                                prior = -99 * np.ones(100)

                            posterior = dm.get_mcmc('median', k)
                            lower = dm.get_mcmc('lower_ui', k)
                            upper = dm.get_mcmc('upper_ui', k)
                            if len(posterior) == 0:
                                posterior = -99 * np.ones(100)
                                lower = -99 * np.ones(100)
                                upper = -99 * np.ones(100)
                            for a in range(100):
                                X.append([t, r, s, y, a,
                                         prior[a],
                                         posterior[a],
                                         upper[a],
                                         lower[a]
                                         ])

            f = default_storage.open(p_csv.file, 'w')
            csv.writer(f).writerows(X)
            f.close()

            convert_cmd = 'echo \'library(foreign); X=read.csv("%s"); write.dta(X, "%s")\' | %s --no-save' % (p_csv.file.path, p_dta.file.path, dismod3.settings.R_PATH)
            ret = subprocess.call(convert_cmd, shell=True)
            assert ret == 0, 'return code %d' % ret
        
        return HttpResponse(open(p_dta.file.path).read(), mimetype='application/x-stata')
    else:
        raise Http404


@login_required
def dismod_show_by_region(request, id, region, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)
    return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_plot',
                                        args=(id, dm.condition, 'all', region, 'all', 'all', format)))

    
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

    style = 'sparkplot'
    plot_key = '%s-plot-%s' % (style, format)  # make sure plot is in the name of all cached plots (and not other things)
    filter = dm.params.filter(key=plot_key)

    if filter.count() != 0:
        plot = filter[0]

    else:
        if format in ['png', 'svg', 'eps', 'pdf']:
            dismod3.sparkplot_disease_model(dm.to_json())

            # save the results of the plot for faster access
            plot = DiseaseModelParameter(key=plot_key)
            fig_data = view_utils.figure_data(format)
            fname = '%s/%s_%s+%s.%s' % (dm.id, style, dm.id, dm.condition, format)
            dm.save_param_data(plot, fname, fig_data)
        else:
            raise Http404

    # return the plot (which is now cached)
    return HttpResponse(open(plot.file.path).read(),
                        view_utils.MIMETYPE[format])

@login_required
def dismod_plot(request, id, condition, type, region, year, sex, format='png', style='tile'):
    if not format in ['png', 'svg', 'eps', 'pdf']:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)
    plot_key = '%s-plot-%s'%(style,format)

    # if the request includes options for the plot (ymax=.01, fontsize=10, etc...) regenerate the plot
    if request.GET:
        for p in dm.params.filter(key=plot_key, type=type, region=region, year=year, sex=sex):
            p.delete()

    filter = dm.params.filter(key=plot_key, type=type, region=region, year=year, sex=sex)

    if filter.count() != 0:
        plot = filter[0]

    else:
        # generate the plot with matplotlib (using code in dismod3.plotting)
        param_filter = dict(region=region, year=year, sex=sex)
        if type == 'all':
            if region == 'all':
                keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
            elif sex == 'all' and year == 'all':
                # plot tiles for each year and sex
                keys = dismod3.utils.gbd_keys(region_list=[region], year_list=['1990', '2005'], sex_list=['male', 'female'])
                param_filter = dict(region=region)
            else:
                keys = dismod3.utils.gbd_keys(region_list=[region], year_list=[year], sex_list=[sex])
        else:
            keys = dismod3.utils.gbd_keys(type_list=[type], region_list=[region], year_list=[year], sex_list=[sex])

        pl.title('%s; %s; %s; %s' % (dismod3.plotting.prettify(condition),
                                     dismod3.plotting.prettify(region), year, sex))
        if style == 'tile':
            dismod3.tile_plot_disease_model(dm.to_json(param_filter), keys, defaults=request.GET)
        elif style == 'overlay':
            dismod3.overlay_plot_disease_model([dm.to_json(param_filter)], keys)
        elif style == 'bar':
            dismod3.bar_plot_disease_model(dm.to_json(param_filter), keys)
        elif style == 'sparkline':
            dismod3.plotting.sparkline_plot_disease_model(dm.to_json(param_filter), keys)
        else:
            raise Http404

        # save the results of the plot for faster access
        plot = DiseaseModelParameter(key=plot_key, type=type, region=region, year=year, sex=sex)
        fig_data = view_utils.figure_data(format)
        fname = '%s/%s_%s+%s+%s+%s+%s+%s.%s' % (dm.id, style, dm.id, dm.condition, type, region, year, sex, format)
        dm.save_param_data(plot, fname, fig_data)

    # return the plot (which is now cached)
    return HttpResponse(open(plot.file.path).read(),
                        view_utils.MIMETYPE[format])

@login_required
def dismod_summary(request, id, format='html'):
    if not format in ['html']:
        raise Http404
    dm = get_object_or_404(DiseaseModel, id=id)
    data_counts, total = count_data(dm)
        
    if format == 'html':
        dm.px_hash = dismod3.sparkplot_boxes(dm.to_json())
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total})
    else:
        raise Http404

def count_data(dm):
    filter = dm.params.filter(key='data_counts')
    if filter.count() != 0:
        data_counts = json.loads(filter[0].json)
    else:
        data = dm.data.all()
        data_counts = []
        for r in dismod3.gbd_regions:
            c = {}

            c['region'] = r
            c['clean_region'] = clean(r)

            for type, data_type in [['i', 'incidence data'],
                                    ['p', 'prevalence data'],
                                    ['r', 'remission data'],
                                    ['em', 'excess-mortality data']]:
                c[type] = \
                    len([d for d in data if d.relevant_to(data_type, r, year='all', sex='all')])

            # also count relative-risk, mortality, and smr data as excess mortality data
            type = 'em'
            for data_type in ['relative-risk data', 'smr data', 'mortality data']:
                c[type] += \
                        len([d for d in data if d.relevant_to(data_type, r, year='all', sex='all')])

            c['total'] = c['i'] + c['p'] + c['r'] + c['em']


            data_counts.append(c)
        data_counts = sorted(data_counts, reverse=True,
                             key=lambda c: (c['total'], c['region']))
        p = DiseaseModelParameter(key='data_counts', json=json.dumps(data_counts))
        p.save()
        dm.params.add(p)
    total = {}
    for type in ['i', 'p', 'r', 'em']:
        total[type] = sum([d[type] for d in data_counts])
    return data_counts, total

@login_required
def dismod_show_map(request, id):
    year = request.POST.get('year')
    sex = request.POST.get('sex')
    map = request.POST.get('map')
    type = request.POST.get('type')
    age_from = request.POST.get('age_from')
    age_to = request.POST.get('age_to')
    weight = request.POST.get('weight')
    count = request.POST.get('count')
    if sex == 'all' and map != 'data':
        sex = 'total'
    if request.POST.get('data_count') == 'Show Map':
        map = 'data_count'
    
    dm = get_object_or_404(DiseaseModel, id=id)

    age_start = 0
    age_end = 100
    if map != 'data_count':
        try:
            age_start = int(age_from)
            age_end = int(age_to)
        except ValueError:
            error = 'Ages must be integers.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'error': error})

        if age_start < 0 or age_end > 100 or age_start > age_end:
            error = 'Ages must be in [0, 100] and ages must be either equal or in ascending order.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'error': error})

    population_world = []
    if weight == 'world':
        for age in range(age_start, age_end + 1):
            population_age = 0
            for region in dismod3.gbd_regions:
                population_age += population_by_region_year_sex(clean(region), year, sex)[age]
            population_world.append(population_age)

    data = dm.data.all()
    vals = {}
    data_type = 'float'
    for r in dismod3.gbd_regions:
        if map == 'data_count':
            vals[clean(r)] = len([d for d in data if d.relevant_to(type=count, region=r, year='all', sex='all')])
            data_type = 'int'
        elif map == 'data':
            d_list = [d for d in data if d.relevant_to(type=type + ' data', region=r, year=year, sex=sex)]
            age_n = age_end - age_start + 1 
            rate = np.empty([age_n])
            for i in range(age_n):
                rate[i] = 'nan'
            for i in range(age_n):
                n = 0
                for j in range(len(d_list)):
                    if d_list[j].age_start <= age_start + i and d_list[j].age_end >= age_start + i:
                        if str(rate[i]) == 'nan':
                            rate[i] = d_list[j].value / float(d_list[j].params['units'])
                        else:
                            rate[i] += d_list[j].value / float(d_list[j].params['units'])
                        n += 1
                if n > 0:
                    rate[i] /= n

            if weight == 'direct':
                sum = 0
                n = 0
                for i in range(age_n):
                    if str(rate[i]) != 'nan':
                        sum += rate[i]
                        n += 1
                if n > 0:
                    vals[clean(r)] = sum / n
                else:
                    vals[clean(r)] = 'nan'
            elif weight == 'region':
                population_region = population_by_region_year_sex(clean(r), year, sex)[age_start:age_end + 1]
                population_sum = 0
                data_sum = 0
                for i in range(len(rate)):
                    if str(rate[i]) != 'nan':
                        data_sum += rate[i] * population_region[i]
                        population_sum += population_region[i]
                if population_sum > 0:
                    vals[clean(r)] = data_sum / population_sum
                else:
                    vals[clean(r)] = 'nan'
            elif weight == 'world':
                population_sum = 0
                data_sum = 0
                for i in range(len(rate)):
                    if str(rate[i]) != 'nan':
                        data_sum += rate[i] * population_world[i]
                        population_sum += population_world[i]
                if population_sum > 0:
                    vals[clean(r)] = data_sum / population_sum
                else:
                    vals[clean(r)] = 'nan'
 
        elif map == 'emp-prior':
            if dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'})).get_empirical_prior(type) != 'empty':
                priors = dict([[p.key, json.loads(json.loads(p.json))] for p in dm.params.filter(key__contains='empirical_prior')])
                if priors == 'empty':
                    return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})
                try:
                    rate = dismod3.neg_binom_model.predict_region_rate('%s+%s+%s+%s' % (type, clean(r), year, sex),
                                               priors['empirical_prior_' + type]['alpha'],
                                               priors['empirical_prior_' + type]['beta'],
                                               priors['empirical_prior_' + type]['gamma'],
                                               dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'})).get_covariates())[age_start:age_end + 1]
                    set_region_value_dict(vals, r, rate, weight, year, sex, age_start, age_end, population_world)
                except KeyError:
                    return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})
            else:
                return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})

        elif map == 'posterior':
            try:
                t = type
                if type == 'with-condition-mortality':
                    t = 'mortality'
                rate = []
                if sex == 'total':
                    rate_m = dismod3.disease_json.DiseaseJson(dm.to_json()).get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, 'male'))[age_start:age_end + 1]
                    if len(rate_m) > 0:
                        rate_f = dismod3.disease_json.DiseaseJson(dm.to_json()).get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, 'female'))[age_start:age_end + 1]
                        if len(rate_f) > 0:
                            population_m = population_by_region_year_sex(clean(r), year, 'male')[age_start:age_end + 1]
                            population_f = population_by_region_year_sex(clean(r), year, 'female')[age_start:age_end + 1]
                            for i in range(age_end - age_start + 1):
                                rate.append((rate_m[i] * population_m[i] + rate_f[i] * population_f[i]) / (population_m[i] + population_f[i]))
                else:
                    rate = dismod3.disease_json.DiseaseJson(dm.to_json()).get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, sex))[age_start:age_end + 1]
                if len(rate) != 0:
                    set_region_value_dict(vals, r, rate, weight, year, sex, age_start, age_end, population_world)
                else:
                    vals[clean(r)] = 'nan'
            except KeyError:
                return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})
        else:
            raise Http404

    title = ''
    if map == 'data_count':
        title = 'Data Count: model #'+id+' '+dm.condition+' (all types, all years, all sexes, all ages)'
    elif map == 'data':
        title = 'Data Value: model #'+id+' '+dm.condition+' ('+type+', '+sex+', age '+str(age_start)+'-'+str(age_end)+', '+year+')'
    elif map == 'emp-prior':
        title = 'Empirical Prior: model #'+id+' '+dm.condition+' ('+type+', '+sex+', age '+str(age_start)+'-'+str(age_end)+', '+year+')'
    elif map == 'posterior':
        title = 'Posterior: model #'+id+' '+dm.condition+' ('+type+', '+sex+', age '+str(age_start)+'-'+str(age_end)+', '+year+')'
    else:
        raise Http404
    if weight == 'region':
        title += ' weighted by region population'
    elif weight == 'world':
        title += ' weighted by world population'
        
    map_info = dismod3.plotting.choropleth_dict(title, vals, data_type=data_type)
    if map_info == None:
        return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})
    return render_to_response('dismod_map.svg',  map_info, mimetype=view_utils.MIMETYPE['svg'])

def set_region_value_dict(vals, region, rate, weight, year, sex, age_start, age_end, population_world):
    if weight == 'direct':
        vals[clean(region)] = np.mean(rate)
    elif weight == 'region':
        population_region = population_by_region_year_sex(clean(region), year, sex)[age_start:age_end + 1]
        data_sum = 0
        population_sum = 0
        for i in range(len(rate)):
            data_sum += rate[i] * population_region[i]
            population_sum += population_region[i]
        vals[clean(region)] = data_sum / population_sum
    elif weight == 'world':
        population_sum = np.sum(population_world)
        data_sum = 0
        for i in range(len(rate)):
            data_sum += rate[i] * population_world[i]
        vals[clean(region)] = data_sum / population_sum

@login_required
def dismod_show_emp_priors(request, id, format='html', effect='alpha'):
    if not format in ['html', 'json', 'png', 'svg', 'eps', 'pdf', 'csv']:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)
    priors = dict([[p.key, json.loads(json.loads(p.json))] for p in dm.params.filter(key__contains='empirical_prior')])

    if format == 'json':
        return HttpResponse(json.dumps(priors),
                            view_utils.MIMETYPE[format])
    elif format == 'csv':
        X_head = 'type, param, index, value'.split(', ')
        X = []
        for t, p in priors.items():
            for param, vals in p.items():
                if type(vals) == list:
                    for age, val in enumerate(vals):
                        X.append([t, param, age, val])
                else:
                    X.append([t, param, '', vals])
        return HttpResponse(view_utils.csv_str(X_head, X), view_utils.MIMETYPE[format])

    elif format in ['png', 'svg', 'eps', 'pdf']:
        dm = dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'}))
        dismod3.plotting.plot_empirical_prior_effects([dm], effect)
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])

    elif format == 'html':
        return render_to_response('dismod_show_emp_priors.html', {'dm': dm})
    
    else:
        raise Http404

@login_required
def dismod_compare(request):
    id1 = request.GET.get('id1')
    id2 = request.GET.get('id2')
    if not id1:
        filter = DiseaseModel.objects.all().order_by('-id')
        return render_to_response('dismod_compare.html', {'paginated_models': view_utils.paginated_models(request, filter)})
    elif not id2:
        dm1 = get_object_or_404(DiseaseModel, id=id1)
        dm1.notes()
        filter = DiseaseModel.objects.filter(condition=dm1.condition).order_by('-id')
        paginated_models = view_utils.paginated_models(request, filter)
        return render_to_response('dismod_compare.html', {'id1': id1, 'dm': dm1,
                                                          'paginated_models': paginated_models})
    else:
        dm1 = get_object_or_404(DiseaseModel, id=id1)
        dm2 = get_object_or_404(DiseaseModel, id=id2)

        return render_to_response('dismod_comparison.html', {'id1': id1, 'id2': id2, 'dm': dm1, 'regions': [clean(r) for r in dismod3.gbd_regions]})
    

@login_required
def dismod_comparison_plot(request, id1=-1, id2=-1, type='alpha', format='png'):
    if not format in ['png', 'svg', 'eps', 'pdf']:
        raise Http404

    dm1 = get_object_or_404(DiseaseModel, id=id1)
    dm2 = get_object_or_404(DiseaseModel, id=id2)

    if type in ['alpha', 'beta', 'gamma', 'delta']:
        dm_list = [dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'})) for dm in [dm1, dm2]]
        dismod3.plotting.plot_empirical_prior_effects(dm_list, type)
    elif type.startswith('overlay'):
        plot_type, rate_type, region, year, sex = type.split('+')
        dm_list = [dismod3.disease_json.DiseaseJson(dm.to_json({'region': region, 'sex': sex, 'year': year})) for dm in [dm1, dm2]]
        dismod3.overlay_plot_disease_model(dm_list, ['%s+%s+%s+%s' % (rate_type, region, year, sex)], defaults=request.GET)

    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])

    
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
        if not model_dict.has_key('id'):
            raise forms.ValidationError('missing model id' % key)

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
            id = model_dict['id']
            if id > 0:
                dm = get_object_or_404(DiseaseModel, id=id)
                for key,val in model_dict['params'].items():
                    if isinstance(val, dict):
                        for subkey in val:
                            t,r,y,s = dismod3.type_region_year_sex_from_key(subkey)
                            if t != 'unknown':
                                param, flag = dm.params.get_or_create(key=key, type=t, region=r, sex=s, year=y)
                                param.json = json.dumps(val[subkey])
                            else:
                                param, flag = dm.params.get_or_create(key=key, type=t)

                                pd=json.loads(param.json)
                                pd[subkey] = val[subkey]
                                param.json = json.dumps(pd)

                            param.save()

                    else:
                        param, flag = dm.params.get_or_create(key=key)
                        param.json = json.dumps(val)
                        param.save()
            else:
                dm = create_disease_model(form.cleaned_data['model_json'], request.user)

            # TODO: clear cache of images by type, region, year, and sex
            # (and test that it works)
            for p in dm.params.filter(key__contains='plot'):
                p.delete()

            # remember to clear anything else that is cached as a param file here too

            return HttpResponseRedirect(dm.get_absolute_url()) # Redirect after POST

    return render_to_response('dismod_upload.html', {'form': form})

@login_required
def job_queue_list(request):
    # accept format specified in url
    format = request.GET.get('format', 'html')

    to_run_list = DiseaseModelParameter.objects.filter(key='needs_to_run')
    if format == 'json':
        return HttpResponse(json.dumps([ param.id for param in to_run_list ]),
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
            param = get_object_or_404(DiseaseModelParameter, id=form.cleaned_data['id'])
            param_val = json.loads(param.json)

            param.key = 'run_status'
            param_val['run_status'] = '%s started at %s' % (param_val.get('estimate_type', ''), time.strftime('%H:%M on %m/%d/%Y'))
            param.json = json.dumps(param_val)
            param.save()

            return HttpResponse(param.json, view_utils.MIMETYPE['json'])
    return render_to_response('job_queue_remove.html', {'form': form})

@login_required
def job_queue_add(request, id):
    # only react to POST requests
    if request.method != 'POST':
        raise Http404
    dm = get_object_or_404(DiseaseModel, id=id)

    # TODO: add logic here for emp prior vs. posterior runs, and for enqueuing selected region/year/sex 

    param = DiseaseModelParameter(key='needs_to_run')
    param_val = {}
    param_val['dm_id'] = id
    # TODO: add details of region/year/sex to param_val dict
    param_val['estimate_type'] = request.POST.get('estimate_type', '')
    estimate_type = ''
    if param_val['estimate_type'].find('posterior') != -1:
        estimate_type = 'posterior'

        if request.POST.get('requested_by', '') == 'run_page':
            # check if priors estimation is completed
            dir_log = dismod3.settings.JOB_LOG_DIR % int(id)
            filename = '%s/%s/status' % (dir_log, 'empirical_priors')
            if os.path.exists(filename):
                f = open(filename, 'r')
                status = f.read()
                f.close()
                if status.find('prevalence::Completed') == -1 or status.find('incidence::Completed') == -1 or \
                    status.find('remission::Completed') == -1 or status.find('excess-mortality::Completed') == -1:
                    error = 'The empirical priors estimation has not been completed.'
                    return render_to_response('dismod_run.html', {'dm': dm, 'error': error})
            else:
                error = 'The empirical priors have not been estimated.'
                return render_to_response('dismod_run.html', {'dm': dm, 'error': error})

        param_val['regions_to_fit'] = []
        for key in request.POST:
            if key == 'all_regions':
                param_val['regions_to_fit'].append(key)
        if(len(param_val['regions_to_fit']) == 0):
            for key in request.POST:
                if key != 'estimate_type' and key != 'requested_by':
                    param_val['regions_to_fit'].append(key)

        if request.POST.get('requested_by', '') == 'run_page':
            if len(param_val['regions_to_fit']) == 0:
                error = 'Please select at least one GBD region.'
                return render_to_response('dismod_run.html', {'dm': dm, 'error': error})
         
    elif param_val['estimate_type'].find('empirical priors') != -1:
        estimate_type = 'empirical_priors'
    elif param_val['estimate_type'].find('') != -1:
        estimate_type = 'fit each region/year/sex individually'
    else:
        error = 'unrecognized estimate type: %s' % estimate_type
        return render_to_response('dismod_run.html', {'dm': dm, 'error': error})

    param_val['run_status'] = '%s queued at %s' % (param_val['estimate_type'], time.strftime('%H:%M on %m/%d/%Y'))
    param.json = json.dumps(param_val)

    d = '%s/%s' % (dismod3.settings.JOB_LOG_DIR % int(id), estimate_type)
    if os.path.exists(d):
        rmtree(d)

    param.save()
    dm.params.add(param)

    if request.POST.get('requested_by', '') == 'run_page':
        return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_show_status', args=[dm.id]) + '?estimate_type=%s' % estimate_type + '&called_by=auto')
    else:
        return HttpResponse('none')

@login_required
def dismod_run(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    error = ''
    return render_to_response('dismod_run.html', {'dm': dm, 'error': error})

@login_required
def dismod_show_status(request, id):
    dir_log = JOB_LOG_DIR % int(id)
    dir_working = JOB_WORKING_DIR % int(id)
    if request.method == 'GET':
        dm = get_object_or_404(DiseaseModel, id=id)
        estimate_type = request.GET.get('estimate_type', 1)
        called_by = request.GET.get('called_by', 2)
        filename = '%s/%s/status' % (dir_log, estimate_type)
        status = 'unavailable'
        if os.path.exists(filename):
            files = os.listdir('%s/%s/stderr' % (dir_working, estimate_type))
            for x in files:
                p = '%s/%s/stderr/%s' % (dir_working, estimate_type, x)
                if os.path.getsize(p) > 0:
                    f = open(filename, 'a+')
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                    if f.read().find('%s::Failed' % x) == -1:
                        f.write('%s::Failed::%s\n' % (x, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.stat(p).st_atime))))
                    f.close()
            f = open(filename, 'r')
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            status = f.read()
            f.close()
            if status == '':
                status = 'none'
        return render_to_response('dismod_show_status.html', {'dm': dm, 'estimate_type': estimate_type, 'status': status, 'called_by': called_by, 'sessionid': request.COOKIES['sessionid']})
    elif request.method == 'POST':
        estimate_type = request.POST['estimate_type']
        filename = '%s/%s/status' % (dir_log, estimate_type)
        status = 'unavailable'
        if os.path.exists(filename):
            files = os.listdir('%s/%s/stderr' % (dir_working, estimate_type))
            for x in files:
                p = '%s/%s/stderr/%s' % (dir_working, estimate_type, x)
                if os.path.getsize(p) > 0:
                    f = open(filename, 'a+')
                    fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                    if f.read().find('%s::Failed' % x) == -1:
                        f.write('%s::Failed::%s\n' % (x, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.stat(p).st_atime))))
                    f.close()
            f = open(filename, 'r')
            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            status = f.read()
            f.close()
            if status == '':
                status = 'none'
        task = request.POST['TASK']
        if task == '':
            stdout = 'Fitting task not selected'
            stderr = 'Fitting task not selected'
        else:
            filename = '%s/%s/stdout/%s' % (dir_working, estimate_type, task)
            stdout = 'unavailable'
            if os.path.exists(filename):
                f = open(filename, 'r')
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                stdout = f.read()
                f.close()
                if stdout == '':
                    stdout = 'none'
            filename = '%s/%s/stderr/%s' % (dir_working, estimate_type, task)
            stderr = 'unavailable'
            if os.path.exists(filename):
                f = open(filename, 'r')
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                stderr = f.read()
                f.close()
                if stderr == '':
                    stderr = 'none'
        return HttpResponse('%s&&&%s&&&%s' % (status, stdout, stderr))

@login_required
def dismod_export(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    return render_to_response('dismod_export.html', {'dm': dm})

@login_required
def dismod_update_covariates(request, id):
    # this is slow, so it has been spun off into a separate view that
    # can be run only when necessary

    # only react to POST requests, since this may changes database data
    if request.method != 'POST':
        raise Http404
    dm = get_object_or_404(DiseaseModel, id=id)
    cov_f = dm.params.filter(key='covariates')
    if len(cov_f) != 0:
        cov_dict = json.loads(cov_f[0].json)
    else:
        cov_dict = {'Country_level': {}}
    
    for d in dm.data.all():
        d.age_weights()  # will cache value if it is not already cached

        for cov_type in cov_dict['Country_level']:
            if cov_dict['Country_level'][cov_type]['rate']['value'] == 1:
                d.calculate_covariate(cov_type)
    
    return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_run', args=[dm.id])) # Redirect after POST

@login_required
def dismod_set_covariates(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    if request.method == 'GET':
        covariates, is_new = dm.params.get_or_create(key='covariates')
        if is_new:
            # extract covariates from data and save them in covariate json
            covariates.json = json.dumps(
                {'Study_level': dm.study_level_covariates(),
                 'Country_level': dm.country_level_covariates()
                 }
                )
            covariates.save()
        return render_to_response('dismod_set_covariates.html', {'dm': dm, 'sessionid': request.COOKIES['sessionid'], 'covariates': covariates})
    elif request.method == 'POST':
        dj = dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'}))

        #exclude fit specific keys from new model
        for key in dj.params.keys():
            if key.find('empirical_prior_') == 0 or key.find('mcmc_') == 0 or key == 'map' or key == 'initial_value':
                dj.params.pop(key)

        cov = json.loads(request.POST['JSON'].replace('\n', ''))
        dj.set_covariates(cov)
        new_dm = create_disease_model(dj.to_json(), request.user)
        
        return HttpResponse(reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))

@login_required
def dismod_adjust_priors(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    if request.method == 'GET':
        return render_to_response('dismod_adjust_priors.html', {'dm': dm, 'global_priors': dm.params.filter(key='global_priors'), 'sessionid': request.COOKIES['sessionid']})
    elif request.method == 'POST':
        dj = dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'}))

        #exclude fit specific keys from new model
        for key in dj.params.keys():
            if key.find('empirical_prior_') == 0 or key.find('mcmc_') == 0 or key == 'map' or key == 'initial_value':
                dj.params.pop(key)

        new_dm = create_disease_model(dj.to_json(), request.user)

        global_priors, flag = new_dm.params.get_or_create(key='global_priors')
        global_priors.json = json.dumps(json.loads(request.POST['JSON']))
        global_priors.save()
        new_dm.params.add(global_priors)
        
        return HttpResponse(reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))

@login_required
def dismod_preview_priors(request, id, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)
    dm = dismod3.disease_json.DiseaseJson(dm.to_json({'region': 'none'}))
    
    if request.method == 'POST':
        dm.params['global_priors_json'] = request.POST['JSON']
        dm.params['global_priors'] = json.loads(request.POST['JSON'])
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
        s += 'heterogeneity %s, ' % dict[conf_key]
    if dict.get(zero_before_key):
        s += 'zero 0 %s, ' % dict[zero_before_key]
    if dict.get(zero_after_key):
        s += 'zero %s %d, ' % (dict[zero_after_key], dismod3.utils.MAX_AGE)

    return s

def dismod_init_log(request, id, estimate_type, param_id):
    dir_log = dismod3.settings.JOB_LOG_DIR % int(id)
    d = '%s/%s' % (dir_log, estimate_type)
    if not os.path.exists(d):
        os.makedirs(d)
    filename = '%s/status' % d
    if os.path.exists(filename):
        os.remove(filename)
    f = open(filename, 'a')
    fcntl.flock(f.fileno(), fcntl.LOCK_EX)
    if estimate_type == 'posterior':
        param = get_object_or_404(DiseaseModelParameter, id=param_id)
        regions_to_fit = json.loads(param.json)['regions_to_fit']
        if regions_to_fit[0] == 'all_regions':
            regions_to_fit = dismod3.gbd_regions
        f.write('%d\n' % (len(regions_to_fit) * len(dismod3.gbd_sexes) * len(dismod3.gbd_years)))
        for r in regions_to_fit:
            for s in dismod3.gbd_sexes:
                for y in dismod3.gbd_years:
                    f.write('%s+%s+%s::Queued::%s\n' % (clean(r), s, y, time.strftime("%Y-%m-%d %H:%M:%S")))
    elif estimate_type == 'within_each_region':
        f.write('%d\n' % len(dismod3.gbd_regions))
        for r in dismod3.gbd_regions:
            f.write('%s::Queued::%s\n' % (clean(r), time.strftime("%Y-%m-%d %H:%M:%S")))
    elif estimate_type == 'across_all_regions':
        f.write('1\n')
        f.write('all_regions::Queued::%s\n' % (time.strftime("%Y-%m-%d %H:%M:%S")))
    elif estimate_type == 'empirical_priors':
        f.write('4\n')
        for t in ['excess-mortality', 'remission', 'incidence', 'prevalence']:
            f.write('%s::Queued::%s\n' % (t, time.strftime("%Y-%m-%d %H:%M:%S")))
    f.close()
    return HttpResponse('')

def dismod_log_status(request, id, estimate_type, fitting_task, state):
    dir_log = dismod3.settings.JOB_LOG_DIR % int(id)
    f = open('%s/%s/status' % (dir_log, estimate_type), 'a')
    fcntl.flock(f.fileno(), fcntl.LOCK_EX)
    f.write('%s::%s::%s\n' % (fitting_task.replace('--', '+'), state, time.strftime("%Y-%m-%d %H:%M:%S")))
    f.close()
    return HttpResponse('')

def dismod_server_load(request):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect((SERVER_LOAD_STATUS_HOST, SERVER_LOAD_STATUS_PORT))
    data = ''
    while(1):
        income = s.recv(SERVER_LOAD_STATUS_SIZE)
        if income == '':
            break;
        data = '%s%s' % (data, income.strip())
    s.close()
    return HttpResponse(data)

