from django.http import *
from django.shortcuts import render_to_response, get_object_or_404

import pymc.gp as gp
import numpy as np
import pylab as pl
import simplejson as json

from gbd.population_data_server.models import Population
import gbd.view_utils as view_utils

def population_show(request, id, format='png'):
    """ Serve a representation of the selected population curve

    Parameters::

      id : int
        the id of the population curve to display
      format : str, optional
        the format to return the results in, may be one of the following:
        json, csv, png, pdf

    Results:

    An HTTP response with population data in the requested form

    Notes:

    It would be cool if this were selected by region and year, and
    showed a population pyramid; and even cooler if it would give a
    visual representation of how the population curve changes over
    time.
    """
    
    pop = get_object_or_404(Population, pk=id)

    M,C = pop.gaussian_process()
    
    x = np.arange(0.,100.,1.)
    p = np.maximum(0., M(x))

    if format == 'json':
        response = json.dumps({'age': x, 'population': p})
    elif format == 'csv':
        headings = ['Age (years)', 'Population (thousands)']
        rows = [[age, val] for age, val in zip(x, p)]
        response = view_utils.csv_str(headings, rows)
    else:
        view_utils.clear_plot()

        # plot bars
        params = {}
        params['left'] = pop.params['left']
        params['width'] = pop.params['width']

        params['height'] = pop.params['vals']
        
        color = '#5cbe5c' # light green
        #color = '#be5c5c' # light red
        params['color'] = color
        params['edgecolor'] = color
        pl.bar(**params)

        # plot interpolated curve
        pl.plot(x, p, linewidth=4, alpha=.5, color='black')

        view_utils.label_plot("%s, %d, %s" % (pop.region, pop.year, pop.sex))
        pl.ylabel('Population (thousands)')
        response = view_utils.figure_data(format)
    
    return HttpResponse(response, view_utils.MIMETYPE[format])

