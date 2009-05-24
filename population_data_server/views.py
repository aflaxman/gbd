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

    It would be cool if this were selected by region and year; and
    even cooler if it would give a visual representation of how the
    population curve changes over time.
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

        if pop.sex == 'total':
            # plot with age on x axis

            # vertical bars
            params = {}
            params['left'] = pop.params['interval_start']
            params['width'] = pop.params['interval_length']

            params['height'] = pop.params['vals']
        
            color = '#5cbe5c' # light green
            params['color'] = color
            params['edgecolor'] = color
            pl.bar(**params)

            # interpolated curve
            pl.plot(x, p, linewidth=4, alpha=.75, color='#126612')

            view_utils.label_plot("%s, %d, %s" % (pop.region, pop.year, pop.sex))
            pl.ylabel('Population (thousands)')
        else:
            # plot as population pyramid (with age on y axis)
            if pop.sex == 'male':
                male_pop = pop
                mp = p
                female_pop = Population.objects.get(region=pop.region,
                                                    year=pop.year,
                                                    sex='female')
                M,C = female_pop.gaussian_process()
                fp = np.maximum(0., M(x))
            else:
                female_pop = pop
                fp = p
                male_pop = Population.objects.get(region=pop.region,
                                                    year=pop.year,
                                                    sex='male')
                M,C = male_pop.gaussian_process()
                mp = np.maximum(0., M(x))

            # horizontal bars
            params = {}
            params['bottom'] = male_pop.params['interval_start']
            params['height'] = male_pop.params['interval_length']

            params['width'] = -np.array(male_pop.params['vals'])
        
            color = '#5c5cbe' # light blue
            params['color'] = color
            params['edgecolor'] = color
            pl.barh(**params)

            params = {}
            params['bottom'] = female_pop.params['interval_start']
            params['height'] = female_pop.params['interval_length']

            params['width'] = np.array(female_pop.params['vals'])
        
            color = '#be5c5c' # light red
            params['color'] = color
            params['edgecolor'] = color
            pl.barh(**params)
            
            # interpolated curves
            pl.plot(-mp, x, linewidth=4, alpha=.5, color='#121266')
            pl.plot(fp, x, linewidth=4, alpha=.5, color='#661212')

            view_utils.label_plot("%s, %d" % (pop.region, pop.year))
            pl.xlabel('Population (thousands)')
            loc, labels = pl.xticks()
            pl.xticks(loc, np.abs(loc))
            
        response = view_utils.figure_data(format)
    
    return HttpResponse(response, view_utils.MIMETYPE[format])

