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

from forms import *


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
                dj = dm.to_djson(region='none')
                
                #exclude fit specific keys from new model
                for key in dj.params.keys():
                    if key.find('empirical_prior_') == 0 or key.find('mcmc_') == 0 or key == 'map' or key == 'initial_value':
                        dj.params.pop(key)

                    # also remove any cached data counts
                    if key == 'data_counts':
                        dj.params.pop(key)

                # TODO:  extract this into a function for DiseaseJson/DiseaseModel, since it is a sometimes helpful task
                # and include it as a command-line and/or GUI task, to deal with times when the system has a hiccup
                # for key in ['empirical_prior_', 'mcmc_', 'map', 'initial_value', 'plot']:
                #     for p in dm.params.filter(key__contains==key):
                #         try: p.delete()
                #         except: pass
                        
                dm = create_disease_model(dj, request.user)
            else:
                dm = DiseaseModel.objects.create(**args)

            for d in data_list:
                dm.data.add(d)

            # extract covariates from covariate_data_server and save them in covariate json
            covariates, is_new = dm.params.get_or_create(key='covariates')
            covariates.json = json.dumps(
                {'Study_level': dm.study_level_covariates(),
                 'Country_level': dm.country_level_covariates()
                 }
                )
            covariates.save()

            # set expert priors to the defaults
            priors, is_new = dm.params.get_or_create(key='global_priors')
            priors.json = json.dumps(dismod3.settings.default_priors)
            priors.save()
                

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
def dismod_list(request, format='html', show='cur_user'):
    if show == 'all':
        dm_filter = DiseaseModel.objects.all().order_by('-id')
    else:
        dm_filter = DiseaseModel.objects.filter(creator=request.user).order_by('-id')
        if dm_filter.count() == 0:
            dm_filter = DiseaseModel.objects.all().order_by('-id')
    if format == 'html':
        return render_to_response('dismod_list.html',
                                  {'paginated_models': view_utils.paginated_models(request, dm_filter)})
    else:
        raise Http404

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
@login_required
def dismod_show(request, id, format='html'):
    if isinstance(id, DiseaseModel):
        dm = id
    else:
        dm = get_object_or_404(DiseaseModel, id=id)

    if format == 'html':
        dm.px_hash = dismod3.sparkplot_boxes(dm.to_djson(region='none'))
                    
        return render_to_response('dismod_show.html',
                                  {'dm': dm,
                                  'paginated_models': view_utils.paginated_models(request, dm.data.all()), 'page_description': 'Full Data from'})
    elif format == 'json':
        return HttpResponse(dm.to_djson().to_json(), view_utils.MIMETYPE[format])
    elif format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.tile_plot_disease_model(dm.to_djson(),
                                        dismod3.utils.gbd_keys(type_list=dismod3.utils.output_data_types))
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    elif format == 'xls':
        group_size = int(request.GET.get('group_size', 1))

        content = dismod3.table(dm.to_djson(),
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
            dm = dm.to_djson()
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

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
@login_required
def dismod_show_selected_regions(request, id, format='png'):
    params = request.GET
    type = request.GET.get('type')
    year = request.GET.get('year')
    sex = request.GET.get('sex')
    xmin = request.GET.get('xmin')
    xmax = request.GET.get('xmax')
    ymin = request.GET.get('ymin')
    ymax = request.GET.get('ymax')
    grid = request.GET.__contains__('grid')
    linewidth = request.GET.get('linewidth')

    if type == None or year == None or sex == None or xmin == None or xmax == None or ymin == None or ymax == None or linewidth == None:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    selected_regions = []
    for key in request.GET:
        if key == 'all_regions':
            selected_regions = dismod3.gbd_regions
    if(len(selected_regions) == 0):
        for key in request.GET:
            if key != 'type' and key != 'year' and key != 'sex' and key != 'xmin' and key != 'xmax' and key != 'ymin' and key != 'ymax' and key != 'grid' and key != 'linewidth':
                selected_regions.append(str(key).replace('+', ' ').replace('%2C', ',').replace('%2F', '/'))
    if len(selected_regions) == 0:
        raise Http404, 'No region is selected.'

    try:
        xmin = int(xmin)
    except ValueError:
        raise Http404, 'X axis lower bound must be an integer.'

    try:
        xmax = int(xmax)
    except ValueError:
        raise Http404, 'X axis upper bound must be an integer.'

    if xmin < 0 or xmax > 100:
        raise Http404, 'X must be in the range [0, 100].'

    if ymin != 'auto':
        try:
            minY = float(ymin)
        except ValueError:
            raise Http404, 'Y axis lower bound must be a floating point number.'

    if ymax != 'auto':
        try:
            maxY = float(ymax)
        except ValueError:
            raise Http404, 'Y axis upper bound must be a floating point number.'

    try:
        linewidth = float(linewidth)
    except ValueError:
        raise Http404, 'Line Width must be a floating point number.'

    dm_json = dm.to_djson()
    t = type
    if t == 'with-condition-mortality':
        t = 'mortality'
    region_value_dict = {}
    for r in selected_regions:
        rate = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, sex))
        if len(rate) == dismod3.MAX_AGE:
            region_value_dict[r] = rate

    if len(region_value_dict) == 0:
        # views which normally return a graphic should not
        # return html when the parameters are incorrect; this breaks
        # webpages.  instead, return a graphic that says the message
        # this change applies to all message/render_to_response patterns
        message = 'Estimation result of the selected posterior is not available.'
        view_utils.clear_plot()
        view_utils.plot_text(message)
        return HttpResponse(view_utils.figure_data(format), view_utils.MIMETYPE[format])
        
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    dismod3.plotting.plot_posterior_selected_regions(region_value_dict, dm_json.get_param_age_mesh(),
                                                     dm.condition, type, year, sex, dm_json.get_estimate_age_mesh(),
                                                     xmin, xmax, ymin, ymax, grid, linewidth)

    # return the plot (which is now cached)
    return HttpResponse(view_utils.figure_data(format), view_utils.MIMETYPE[format])

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
@login_required
def dismod_show_all_years(request, id, format='png'):
    type = request.GET.get('type')
    region = request.GET.get('region')
    sex = request.GET.get('sex')
    xmin = request.GET.get('xmin')
    xmax = request.GET.get('xmax')
    ymin = request.GET.get('ymin')
    ymax = request.GET.get('ymax')
    grid = request.GET.__contains__('grid')
    linewidth = request.GET.get('linewidth')

    if type == None or region == None or sex == None or xmin == None or xmax == None or ymin == None or ymax == None or linewidth == None:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    region = str(region).replace('+', ' ').replace('%2C', ',').replace('%2F', '/')
    if region == None:
        message = 'No region is selected.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        xmin = int(xmin)
    except ValueError:
        message = 'X axis lower bound must be an integer.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        xmax = int(xmax)
    except ValueError:
        message = 'X axis upper bound must be an integer.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if xmin < 0 or xmax > 100:
        message = 'X must be in the range [0, 100].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymin != 'auto':
        try:
            minY = float(ymin)
        except ValueError:
            message = 'Y axis lower bound must be a floating point number.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})
        if minY < 0:
            message = 'Y axis lower bound must be >= 0.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymax != 'auto':
        try:
            maxY = float(ymax)
        except ValueError:
            message = 'Y axis upper bound must be a floating point number.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})
        if maxY < 0:
            message = 'Y axis upper bound must be >= 0.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymin != 'auto' and ymax != 'auto':
        if  float(ymin) >= float(ymax):
            message = 'Y axis bounds are out of order.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        linewidth = float(linewidth)
    except ValueError:
        message = 'Line Width must be a floating point number.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if linewidth < 0.1 or linewidth > 10:
        message = 'linewidth must be in range [0.1, 10].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if linewidth < 0.1 or linewidth > 10:
        message = 'linewidth must be in range [0.1, 10].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    dm_json = dm.to_djson()
    t = type
    if t == 'with-condition-mortality':
        t = 'mortality'
    ages = dm_json.get_estimate_age_mesh()
    year_value_dict = {}
    for year in ['1990', '2005']:
        if sex == 'total':
            rate_m = dm_json.get_mcmc('mean', '%s+%s+%s+male' % (type, clean(region), year))
            if len(rate_m) == dismod3.MAX_AGE:
                rate_f = dm_json.get_mcmc('mean', '%s+%s+%s+female' % (type, clean(region), year))
                if len(rate_f) == dismod3.MAX_AGE:
                    population_m = population_by_region_year_sex(clean(region), year, 'male')
                    population_f = population_by_region_year_sex(clean(region), year, 'female')
                    year_value_dict[year] = [(r_m * p_m + r_f * p_f) / (p_m + p_f)
                                             for r_m, p_m, r_f, p_f in zip(rate_m, population_m, rate_f, population_f)]
        else:
            rate = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(region), year, sex))
            if len(rate) == dismod3.MAX_AGE:
                year_value_dict[year] = rate

    if len(year_value_dict) == 0:
        message = 'Estimation result of the selected posterior is not available.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    dismod3.plotting.plot_posterior_region(year_value_dict, dm.condition, type, region, sex, ages, xmin, xmax, ymin, ymax, grid, linewidth)

    # return the plot (which is now cached)
    return HttpResponse(view_utils.figure_data(format), view_utils.MIMETYPE[format])

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
@login_required
def dismod_show_all_sexes(request, id, format='png'):
    type = request.GET.get('type')
    region = request.GET.get('region')
    year = request.GET.get('year')
    xmin = request.GET.get('xmin')
    xmax = request.GET.get('xmax')
    ymin = request.GET.get('ymin')
    ymax = request.GET.get('ymax')
    grid = request.GET.__contains__('grid')
    linewidth = request.GET.get('linewidth')

    if type == None or year == None or region == None or xmin == None or xmax == None or ymin == None or ymax == None or linewidth == None:
        raise Http404

    dm = get_object_or_404(DiseaseModel, id=id)

    region = str(region).replace('+', ' ').replace('%2C', ',').replace('%2F', '/')
    if region == None:
        message = 'No region is selected.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        xmin = int(xmin)
    except ValueError:
        message = 'X axis lower bound must be an integer.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        xmax = int(xmax)
    except ValueError:
        message = 'X axis upper bound must be an integer.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if xmin < 0 or xmax > 100:
        message = 'X must be in the range [0, 100].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymin != 'auto':
        try:
            minY = float(ymin)
        except ValueError:
            message = 'Y axis lower bound must be a floating point number.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})
        if minY < 0:
            message = 'Y axis lower bound must be >= 0.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymax != 'auto':
        try:
            maxY = float(ymax)
        except ValueError:
            message = 'Y axis upper bound must be a floating point number.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})
        if maxY < 0:
            message = 'Y axis upper bound must be >= 0.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if ymin != 'auto' and ymax != 'auto':
        if  float(ymin) >= float(ymax):
            message = 'Y axis bounds are out of order.'
            data_counts, total = count_data(dm)
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    try:
        linewidth = float(linewidth)
    except ValueError:
        message = 'Line Width must be a floating point number.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if linewidth < 0.1 or linewidth > 10:
        message = 'linewidth must be in range [0.1, 10].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    if linewidth < 0.1 or linewidth > 10:
        message = 'linewidth must be in range [0.1, 10].'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    dm_json = dm.to_djson()
    t = type
    if t == 'with-condition-mortality':
        t = 'mortality'
    ages = dm_json.get_estimate_age_mesh()
    sex_value_dict = {}
    for sex in ['male', 'female', 'total']:
        if sex == 'total':
            rate_m = dm_json.get_mcmc('mean', '%s+%s+%s+male' % (type, clean(region), year))
            if len(rate_m) == dismod3.MAX_AGE:
                rate_f = dm_json.get_mcmc('mean', '%s+%s+%s+female' % (type, clean(region), year))
                if len(rate_f) == dismod3.MAX_AGE:
                    population_m = population_by_region_year_sex(clean(region), year, 'male')
                    population_f = population_by_region_year_sex(clean(region), year, 'female')
                    sex_value_dict[sex] = [(r_m * p_m + r_f * p_f) / (p_m + p_f)
                                          for r_m, p_m, r_f, p_f in zip(rate_m, population_m, rate_f, population_f)]
        else:
            rate = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(region), year, sex))
            if len(rate) == dismod3.MAX_AGE:
                sex_value_dict[sex] = rate

    if len(sex_value_dict) == 0:
        message = 'Estimation result of the selected posterior is not available.'
        data_counts, total = count_data(dm)
        return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'message': message})

    dismod3.plotting.plot_posterior_region(sex_value_dict, dm.condition, type, region, year, ages, xmin, xmax, ymin, ymax, grid, linewidth)

    # return the plot (which is now cached)
    return HttpResponse(view_utils.figure_data(format), view_utils.MIMETYPE[format])
    
@login_required
def dismod_find_and_show(request, condition, format='html'):
    try:
        dm = DiseaseModel.objects.filter(condition=condition).latest('id')
    except DiseaseModel.DoesNotExist:
        raise Http404
    return dismod_show(request, dm, format)

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
@login_required
def dismod_sparkplot(request, id, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)

    style = 'sparkplot'
    plot_key = '%s-plot-%s' % (style, format)  # make sure plot is in the name of all cached plots (and not other things)
    # if the request includes options regenerate the plot
    if request.GET:
        for p in dm.params.filter(key=plot_key):
            p.delete()

    filter = dm.params.filter(key=plot_key)

    if filter.count() != 0:
        plot = filter[0]
        try:
            return HttpResponse(open(plot.file.path).read(),
                                view_utils.MIMETYPE[format])
        except IOError:
            pass

    if format in ['png', 'svg', 'eps', 'pdf']:
        dismod3.sparkplot_disease_model(dm.to_djson())

        # save the results of the plot for faster access
        plot = DiseaseModelParameter(key=plot_key)
        fig_data = view_utils.figure_data(format)
        fname = '%s/%s_%s+%s.%s' % (dm.id, style, dm.id, dm.condition, format)
        dm.save_param_data(plot, fname, fig_data)
    else:
        raise Http404

    # return the plot (which is now cached)
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])

# TODO: change image caching to use cluster computation when possible,
# and to load dm_json from cluster when recomputing is necessary
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

        try:
            # return the plot (which is now cached)
            return HttpResponse(open(plot.file.path).read(),
                                view_utils.MIMETYPE[format])
        except IOError, e:
            pass
            
    ## generate the plot with matplotlib (using code in dismod3.plotting)
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
        dismod3.tile_plot_disease_model(dm.to_djson(region), keys, defaults=request.GET)
    elif style == 'overlay':
        dismod3.overlay_plot_disease_model([dm.to_djson(region)], keys)
    elif style == 'bar':
        dismod3.bar_plot_disease_model(dm.to_djson(region), keys)
    elif style == 'sparkline':
        dismod3.plotting.sparkline_plot_disease_model(dm.to_djson(region), keys)
    else:
        raise Http404

    # save the results of the plot for faster access
    plot = DiseaseModelParameter(key=plot_key, type=type, region=region, year=year, sex=sex)
    fig_data = view_utils.figure_data(format)
    fname = '%s/%s_%s+%s+%s+%s+%s+%s.%s' % (dm.id, style, dm.id, dm.condition, type, region, year, sex, format)
    try:
        dm.save_param_data(plot, fname, fig_data)
    except:
        pass

    # return the plot (which is now cached)
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])

@login_required
def dismod_summary(request, id, rate_type='all', format='html'):
    if not format in ['html']:
        raise Http404
    dm = get_object_or_404(DiseaseModel, id=id)

    if rate_type == 'all':
        data_counts, total = count_data(dm)
        
        if format == 'html':
            if request.GET.get('recalc_table'):
                dm.params.filter(key='data_counts').delete()
                data_counts, total = count_data(dm)
            dm.px_hash = dismod3.sparkplot_boxes(dm.to_djson('none'))
            return render_to_response('dismod_summary.html', {'dm': dm, 'counts': data_counts, 'total': total, 'page_description': 'Summary of'})

    elif format == 'html':
        ymax = request.GET.get('ymax', 'auto')
        return render_to_response('dismod_subsummary.html', {'dm': dm, 'rate_type': rate_type, 'ymax': ymax})
    
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

            # also count relative-risk, mortality, smr, csmr, pf data as excess mortality data
            type = 'em'
            for data_type in ['relative-risk data', 'smr data', 'mortality data', 'prevalence x excess-mortality data', 'cause-specific mortality data']:
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

# TODO: clean up this view
@login_required
def dismod_show_map(request, id):
    year = request.GET.get('year')
    sex = request.GET.get('sex')
    map = request.GET.get('map')
    type = request.GET.get('type')
    age_from = request.GET.get('age_from')
    age_to = request.GET.get('age_to')
    weight = request.GET.get('weight')
    count = request.GET.get('count')
    scheme = request.GET.get('scheme')
    data_count = request.GET.get('data_count')

    if not ((count != None and data_count != None) or (year != None and sex != None and map != None and type != None and age_from != None and age_to != None and weight != None and scheme != None)):
        raise Http404

    if sex != None and sex == 'all' and map != 'data':
        sex = 'total'
    if data_count != None and data_count == 'Show_Map':
        map = 'data_count'
    
    dm = get_object_or_404(DiseaseModel, id=id)

    if map != None:
        if map == 'emp-prior':
            dm_json = dm.to_djson()
        elif map == 'posterior':
            dm_json = dm.to_djson()

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

    population_world = np.zeros(age_end - age_start + 1)
    if weight == 'world':
        for region in dismod3.gbd_regions:
            population_region = population_by_region_year_sex(clean(region), year, sex)[age_start:age_end + 1]
            for age in range(age_end - age_start + 1):
                population_world[age] += population_region[age]

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
                        if np.isnan(rate[i]):
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
                    if not np.isnan(rate[i]):
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
                    if not np.isnan(rate[i]):
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
                    if not np.isnan(rate[i]):
                        data_sum += rate[i] * population_world[i]
                        population_sum += population_world[i]
                if population_sum > 0:
                    vals[clean(r)] = data_sum / population_sum
                else:
                    vals[clean(r)] = 'nan'
        elif map == 'emp-prior':
            if dm_json.get_empirical_prior(type) != 'empty':
                dm_json.vars = dismod3.neg_binom_model.setup(dm_json, type, [])
                priors = dict([[p.key, json.loads(json.loads(p.json))] for p in dm.params.filter(key__contains='empirical_prior')])
                if priors == 'empty':
                    return render_to_response('dismod_message.html', {'type': type, 'year': year, 'sex': sex, 'map': map})
                try:
                    rate = dismod3.neg_binom_model.predict_region_rate('%s+%s+%s+%s' % (type, clean(r), year, sex),
                                               priors['empirical_prior_' + type]['alpha'],
                                               priors['empirical_prior_' + type]['beta'],
                                               priors['empirical_prior_' + type]['gamma'],
                                               dm_json.get_covariates(),
                                               dm_json.vars['bounds_func'],
                                               dm_json.get_estimate_age_mesh())[age_start:age_end + 1]  # TODO: load from fs instead from db                    
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
                    rate_m = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, 'male'))[age_start:age_end + 1]
                    if len(rate_m) > 0:
                        rate_f = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, 'female'))[age_start:age_end + 1]
                        if len(rate_f) > 0:
                            population_m = population_by_region_year_sex(clean(r), year, 'male')[age_start:age_end + 1]
                            population_f = population_by_region_year_sex(clean(r), year, 'female')[age_start:age_end + 1]
                            for i in range(age_end - age_start + 1):
                                rate.append((rate_m[i] * population_m[i] + rate_f[i] * population_f[i]) / (population_m[i] + population_f[i]))
                else:
                    rate = dm_json.get_mcmc('mean', '%s+%s+%s+%s' % (t, clean(r), year, sex))[age_start:age_end + 1]
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

    map_info = dismod3.plotting.choropleth_dict(title, vals, scheme, data_type=data_type)
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
    priors = dict([[p.key, json.loads(json.loads(p.json))] for p in dm.params.filter(key__contains='empirical_prior')])  # TODO: load from fs instead from db

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
        dm = dm.to_djson(region='none')
        dismod3.plotting.plot_empirical_prior_effects([dm], effect)
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])

    elif format == 'html':
        return render_to_response('dismod_show_emp_priors.html', {'dm': dm, 'page_description': 'Empirical Priors for'})
    
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
        dm = get_object_or_404(DiseaseModel, id=id1)
        dm2 = get_object_or_404(DiseaseModel, id=id2)
        pm = dict(object_list = [dm, dm2])

        return render_to_response('dismod_comparison.html', {'id1': id1, 'id2': id2, 'dm': dm, 'paginated_models': pm, 'regions': [clean(r) for r in dismod3.gbd_regions]})
    

@login_required
def dismod_comparison_plot(request, id1=-1, id2=-1, type='alpha', format='png'):
    if not format in ['png', 'svg', 'eps', 'pdf']:
        raise Http404

    dm1 = get_object_or_404(DiseaseModel, id=id1)
    dm2 = get_object_or_404(DiseaseModel, id=id2)

    if type in ['alpha', 'beta', 'gamma', 'delta']:
        dm_list = [dm.to_djson('none') for dm in [dm1, dm2]]  # TODO: load from fs instead from db
        dismod3.plotting.plot_empirical_prior_effects(dm_list, type)
    elif type.startswith('overlay'):
        plot_type, rate_type, region, year, sex = type.split('+')
        dm_list = [dm.to_djson(region) for dm in [dm1, dm2]]  # TODO: load from fs instead from db
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
                from dismod3.disease_json import DiseaseJson
                dj = DiseaseJson(form.cleaned_data['model_json'])
                dm = create_disease_model(dj, request.user)

            # clear cache of images by type, region, year, and sex
            for p in dm.params.filter(key__contains='plot'):
                p.delete()

            # TODO: remember to clear anything else that is cached as a param file here too

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


# TODO: clean up this view
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

    if param_val['estimate_type'] == 'Fit continuous single parameter model':
        param_val['run_status'] = '%s queued at %s' % (param_val['estimate_type'], time.strftime('%H:%M on %m/%d/%Y'))
        param_val['dir'] = '%s/%s' % (dismod3.settings.JOB_LOG_DIR % int(id), 'spm')
        param.json = json.dumps(param_val)

        param.save()
        dm.params.add(param)
                                
        return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_spm_monitor', args=[dm.id]))


    
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
    # TODO: handle 'error' in a more standard django manner
    error = ''
    return render_to_response('dismod_run.html', {'dm': dm, 'error': error})

# TODO: clean up this view
@login_required
def dismod_show_status(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)

    dir = dismod3.settings.JOB_WORKING_DIR % int(id)
    std = {'out': 'FITTING MODEL\n\n', 'err': ''}
    import glob
    try:
        for path in [dir + '/empirical_priors/std%s/*', dir + '/posterior/std%s/*']:
            for out in ['out', 'err']:
                for fname in glob.glob(path % out):
                    f = open(fname)
                    std[out] += '%40s:   ' % fname.split('/')[-1] + ''.join(f.readlines()[-1:]).strip() + '\n'
                    #std[out] += '\n\n****\n\n' + ''.join(f.readlines()[-10:])
                    f.close()
    except IOError:
        std = {'out': 'FITTING MODEL\n\n', 'err': ''}
    return render_to_response('spm_monitor.html', {'dm': dm, 'stdout': std['out'], 'stderr': std['err']})



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
        # FIXME: putting the session id in the html like this is probably insecure
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
    # FIXME: putting the session id in the html like this is probably insecure
    return render_to_response('dismod_export.html', {'dm': dm, 'sessionid': request.COOKIES['sessionid']})

# TODO: make cov_dict an object, and make it easier to get the needed parts out with appropriate methods
@login_required
def dismod_update_covariates(request, id):
    # this is slow, so it has been spun off into a separate view that
    # can be run only when necessary

    # only react to POST requests, since this changes database data
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

    dm.params.filter(key='derived_covariate').delete()
    for cov_type in cov_dict['Country_level']:
        if cov_dict['Country_level'][cov_type]['rate']['value'] == 1:
            dm.cache_derived_covariates(cov_type)
    
    return HttpResponseRedirect(reverse('gbd.dismod_data_server.views.dismod_run', args=[dm.id])) # Redirect after POST

@login_required
def dismod_set_covariates(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    if request.method == 'GET':
        # FIXME: there should not be duplicate covariate params being created in the first place.
        all_covs = dm.params.filter(key='covariates')
        if len(all_covs) > 1:
            for cv in all_covs[1:]:
                cv.delete()
                
        covariates, is_new = dm.params.get_or_create(key='covariates')
        if is_new:
            # extract covariates from data and save them in covariate json
            covariates.json = json.dumps(
                {'Study_level': dm.study_level_covariates(),
                 'Country_level': dm.country_level_covariates()
                 }
                )
            covariates.save()
        # FIXME: putting the session id in the html like this is probably insecure
        return render_to_response('dismod_set_covariates.html', {'dm': dm, 'sessionid': request.COOKIES['sessionid'], 'covariates': covariates})
    elif request.method == 'POST':
        dj = dm.to_djson(region='none')

        # TODO: this exclude block may be unnecessary when results are stored in filesystem
        #exclude fit specific keys from new model
        for key in dj.params.keys():
            if key.find('empirical_prior_') == 0 or key.find('mcmc_') == 0 or key == 'map' or key == 'initial_value':
                dj.params.pop(key)

        cov = json.loads(request.POST['JSON'].replace('\n', ''))
        dj.set_covariates(cov)
        new_dm = create_disease_model(dj, request.user)
        
        return HttpResponse(reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))

@login_required
def dismod_adjust_priors(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    if request.method == 'GET':
        # FIXME: putting the session id in the html like this is probably insecure
        return render_to_response('dismod_adjust_priors.html', {'dm': dm, 'global_priors': dm.params.filter(key='global_priors'), 'sessionid': request.COOKIES['sessionid']})
    elif request.method == 'POST':
        dj = dm.to_djson('none')

        # TODO: this exclude block may be unnecessary when results are stored in filesystem
        #exclude fit specific keys from new model
        for key in dj.params.keys():
            if key.find('empirical_prior_') == 0 or key.find('mcmc_') == 0 or key == 'map' or key == 'initial_value':
                dj.params.pop(key)

        new_dm = create_disease_model(dj, request.user)

        global_priors, flag = new_dm.params.get_or_create(key='global_priors')
        global_priors.json = json.dumps(json.loads(request.POST['JSON']))
        global_priors.save()
        new_dm.params.add(global_priors)
        
        return HttpResponse(reverse('gbd.dismod_data_server.views.dismod_run', args=[new_dm.id]))

@login_required
def dismod_preview_priors(request, id, format='png'):
    dm = get_object_or_404(DiseaseModel, id=id)
    dm = dm.to_djson(region='none')
    
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

# TODO: clean this up and comment it well
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

@login_required
def dismod_experimental(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    return render_to_response('dismod_experimental.html', {'dm': dm})

@login_required
def dismod_spm_monitor(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)

    dir = dismod3.settings.JOB_WORKING_DIR % int(id)

    try:
        f = file('%s/continuous_spm.stdout' % dir)
        stdout = f.read()
        f.close()
    except IOError, e:
        stdout = 'Warning: could not load output\n%s' % e

    try:
        f = file('%s/continuous_spm.stderr' % dir)
        stderr = f.read()
        f.close()
    except IOError, e:
        stderr = 'Warning: could not load stderr\n%s' % e
        
    return render_to_response('spm_monitor.html', {'dm': dm, 'stdout': stdout, 'stderr': stderr})

@login_required
def dismod_spm_view_results(request, id):
    dm = get_object_or_404(DiseaseModel, id=id)
    return render_to_response('spm_view_results.html', {'dm': dm, 'regions': [dismod3.utils.clean(r) for r in dismod3.settings.gbd_regions]})

