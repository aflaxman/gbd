from django.http import *
from django.shortcuts import render_to_response, get_object_or_404

import pylab as pl

from gbd.covariate_data_server.models import *
from gbd import view_utils

def covariate_show(request, id, format='png'):
    """ Serve a representation of the selected covariate

    Parameters::

      id : int
        the id of the covariate to display
      format : str, optional
        the format to return the results in, may be one of the following:
        json, csv, png, pdf
    """
    cov = get_object_or_404(Covariate, id=id)
    X = pl.array(
        sorted(
        [[c.year, c.value] for c in
         Covariate.objects.filter(type=c.type, iso3=c.iso3, sex=c.sex)
         ]))

    fig_width = 6.
    fig_height = 4.5
    fig = pl.figure(figsize=(fig_width, fig_height), dpi=100)
    pl.plot(X[:,0], X[:,1], '.-')
    pl.ylabel(c.type)
    pl.xlabel('Time (Years)')
    pl.title('%s in %s' % (c.type, c.iso3))
            
    response = view_utils.figure_data(format)
    
    return HttpResponse(response, view_utils.MIMETYPE[format])

