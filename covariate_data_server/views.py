from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

import pylab as pl

from gbd.covariate_data_server.models import *
from gbd import view_utils

def covariate_type_show(request, id):
    """ Show an index page for the selected covariate type

    Parameters::

      id : int
        the id of the covariate to display
    """
    ct = get_object_or_404(CovariateType, id=id)
    
    return render_to_response('covariate_type_show.html',
                              {'ct': ct,
                               'paginated_models': view_utils.paginated_models(request, ct.covariate_set.all().distinct().values('iso3', 'sex'))})

def covariate_show(request, type, iso3, sex, format='png'):
    """ Serve a representation of the selected covariate

    Parameters::

      type : str
        the covariate_type
      iso3 : str
        the country code
      sex : str
        the sex to display ('male', 'female', 'all', or '')
      format : str, optional
        the format to return the results in, may be one of the following:
        json, csv, png, pdf
    """
    ct = get_object_or_404(CovariateType, slug=type)
    X = pl.array(
        sorted(
        [[c.year, c.value] for c in
         ct.covariate_set.filter(iso3=iso3, sex=sex)
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


class NewDataForm(forms.Form):
    type = forms.CharField(max_length=50)
    file  = forms.FileField()
    
    def clean_file(self):
        import csv
        data = [d for d in csv.DictReader(self.file)]
        return data
    
    def clean(self):
        data = self.cleaned_data.get('file')
        type_slug = self.cleaned_data.get('type')
        
        if data and type_slug:
            for i, d in enumerate(data):
                try:
                    d['value'] = float(d[type_slug])
                except KeyError:
                    raise forms.ValidationError('Could not find column %s (is it spelled correctly?)' % type_slug)
                except ValueError:
                    raise forms.ValidationError('Could not interpret value for %s in line %d' % (type_slug, i+2))

                if d.has_key('year'):
                    try:
                        d['year'] = int(d['year'])
                    except ValueError:
                        raise forms.ValidationError('Could not interpret year in line %d' % (ii+2))
                else:
                    d['year'] = gbd.fields.ALL_YEARS
                        

                d['sex'] = d.get('sex', '')
                if not d['sex'] in ['male', 'female', 'total', '']:
                    raise forms.ValidationError('Could not interpret sex in line %d' % (ii+2))

        # Always return the full collection of cleaned data.
        return self.cleaned_data
    
@login_required
def covariate_upload(request):
    if request.method == 'GET':  # no form data is associated with page, yet
        form = NewDataForm()
    elif request.method == 'POST':  # If the form has been submitted...
        form = NewDataForm(request.POST, request.FILES)  # A form bound to the POST data
        form.file = request.FILES.get('file')
        if form.is_valid():
            # All validation rules pass, so create new data based on the
            # form contents
            type_slug = form.cleaned_data['type']
            cov_type, is_new = CovariateType.objects.get_or_create(slug=type_slug)
            cov_data = form.cleaned_data['file']

            # make rates from rate_list
            vals = [d['value'] for d in cov_data]
            mu = pl.mean(vals)
            std = pl.std(vals)
            
            for d in cov_data:
                # if sex == '' add a covariate for male, female, and total
                if d['sex'] == '':
                    sex_list = ['male', 'female', 'total']
                else:
                    sex_list = [d['sex']]
                for sex in sex_list:
                    # add a data point, save it on the data list
                    cov, is_new = Covariate.objects.get_or_create(type=cov_type,
                                                                  iso3=d['iso3'],
                                                                  year=d['year'],
                                                                  sex=sex,
                                                                  defaults={'value': 0.})
                    cov.value = (d['value'] - mu) / std
                    cov.save()

            return HttpResponseRedirect(reverse('gbd.covariate_data_server.views.covariate_type_show', args=[cov_type.id])) # Redirect after POST

    return render_to_response('covariate_upload.html', {'form': form})

