from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse

import pymc.gp as gp
import numpy as np
import pylab as pl

from dismod3.models import Population
import view_utils

@login_required
def population_show(request, id):
    pop = get_object_or_404(Population, pk=id)
    return render_to_response('population/show.html',
                              view_utils.template_params(pop))

@login_required
def population_redirect(request, id, action):
    pop = get_object_or_404(Population, pk=id)
    if action == 'edit':
        url = '/admin/dismod3/population/%d' % pop.id
    elif action in view_utils.command_list['move']:
        url = reverse('dismod3.views.population_show', args=(pop.id+view_utils.id_delta[action],))
    elif action in view_utils.command_list['sex']:
        url = Population.objects.get(country=pop.country, sex=action, year=pop.year).get_absolute_url()
    elif action in view_utils.command_list['format']:
        url = '%s.%s' % (pop.get_absolute_url(), action)
    else:
        raise Http404
    
    return HttpResponseRedirect(url)

@login_required
def population_plot(request, id, format):
    pop = get_object_or_404(Population, pk=id)

    M,C = pop.gaussian_process()
    x = np.arange(0.,100.,1.)

    # handle json & csv formats, which are not technically a plot
    if format in ['json', 'csv']:
        if format == 'json':
            data_str = pop.data_json
        elif format == 'csv':
            headings = ['Age (years)', 'Population (thousands)']
            rows = [[a, p] for a,p in zip(pop.data['mesh'], pop.data['vals'])]
            data_str = view_utils.csv_str(headings, rows)
        return HttpResponse(data_str, view_utils.MIMETYPE[format])
        

    view_utils.clear_plot()
    pl.plot(x,np.maximum(0.0,M(x)))
    view_utils.label_plot("%s, %d, %s" % (pop.country, pop.year, pop.sex))
    pl.ylabel('Population (thousands)')
    
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])
