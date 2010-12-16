from django.shortcuts import render_to_response, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms
from datetime import datetime
import pylab as pl
import csv
from gbd.covariate_data_server.models import *
from gbd import view_utils
from gbd.dismod3.settings import CSV_PATH
from gbd.dismod3.settings import gbd_regions
from gbd.dismod3.utils import clean

@login_required
def covariate_type_list_show(request):
    """ Show the list of all uploaded covariate types """
    ct = CovariateType.objects.all()
    for c in ct:
        if c.region_only:
            c.region_only = 'yes'
        else:
            c.region_only = 'no'

    return render_to_response('covariate_type_list_show.html', {'ct': ct})

def covariate_notes_show(request, id):
    """ Show description for the selected covariate type

    Parameters:
    -----------
      id : int
        the id of the covariate type to display
    """
    ct = get_object_or_404(CovariateType, id=id)
    return render_to_response('covariate_notes_show.html', {'ct': ct})

def covariate_data_count_show(request, id):
    """ Show amount of data for each country of the selected covariate type

    Parameters:
    -----------
      id : int
        the id of the covariate type to display
    """
    ct = get_object_or_404(CovariateType, id=id)
    
    if ct.region_only:
        pm = ct.covariate_set.all().distinct().values('region')

        for c in pm:
            c['clean_region'] = clean(c['region'])
            c['count'] = ct.covariate_set.filter(region=c['region']).count()
            if c['count'] < (ct.year_end - ct.year_start + 1) * 3:
                c['color'] = 'class=highlight'
            else:
                c['color'] = ''
        
        if len(pm) != 21:
            error = 'Total number of regions are wrong.  Found ' + str(len(pm)) + '.  Should be 21.'
        else:
            error = ''

        return render_to_response('covariate_data_count_show.html',
                                  {'ct': ct, 'level': 'region', 'error': error,
                                   'paginated_models': view_utils.paginated_models(request, pm)})
    else:
        pm = ct.covariate_set.all().distinct().values('iso3')

        for c in pm:
            c['count'] = ct.covariate_set.filter(iso3=c['iso3']).count()
            if c['count'] < (ct.year_end - ct.year_start + 1) * 3:
                c['color'] = 'class=highlight'
            else:
                c['color'] = ''

        return render_to_response('covariate_data_count_show.html',
                                  {'ct': ct, 'level': 'country',
                                   'paginated_models': view_utils.paginated_models(request, pm)})

def covariate_type_show(request, id):
    """ Show an index page for the selected covariate type

    Parameters:
    -----------
      id : int
        the id of the covariate type to display
    """
    ct = get_object_or_404(CovariateType, id=id)

    if ct.region_only:
        pm = ct.covariate_set.all().distinct().values('region', 'sex')

        for c in pm:
            c['clean_region'] = clean(c['region'])
            c['count'] = ct.covariate_set.filter(region=c['region'], sex=c['sex']).count()
            if c['count'] < ct.year_end - ct.year_start + 1:
                c['color'] = 'class=highlight'
            else:
                c['color'] = ''

        return render_to_response('covariate_type_show.html',
                                  {'ct': ct, 'level': 'region',
                                   'paginated_models': view_utils.paginated_models(request, pm)})
    else:
        pm = ct.covariate_set.all().distinct().values('iso3', 'sex')

        for c in pm:
            c['count'] = ct.covariate_set.filter(iso3=c['iso3'], sex=c['sex']).count()
            if c['count'] < ct.year_end - ct.year_start + 1:
                c['color'] = 'class=highlight'
            else:
                c['color'] = ''

        return render_to_response('covariate_type_show.html',
                                  {'ct': ct, 'level': 'country',
                                   'paginated_models': view_utils.paginated_models(request, pm)})

def covariate_data_value_show(request, type, area, format='png'):
    """ Serve a representation of the covariate for the specified type and country

    Parameters:
    -----------
      type : str
        the covariate_type
      area : str
        either the country code or the gbd region
      format : str, optional
        the format to return the results in, may be one of the following:
        json, csv, png, pdf
    """
    ct = get_object_or_404(CovariateType, slug=type)
    fig_width = 18.
    fig_height = 4.5
    sexes = ['male', 'female', 'total']
    pl.figure(figsize=(fig_width, fig_height), dpi=100)
    if len(area) == 3:
        for i, s in enumerate(sexes):
            pl.subplot(1, 3, i + 1)
            X = pl.array(
                sorted([[c.year, c.value] for c in ct.covariate_set.filter(iso3=area, sex=s)]))
            if len(X) > 0:
                pl.plot(X[:,0], X[:,1], '.-')
                pl.ylabel(c.type)
                pl.xlabel('Time (Years)')
                pl.title('%s for %s in %s' % (c.type, s, c.iso3))
    else:
        region_dict = {}
        for r in gbd_regions:
            region_dict[clean(r)] = r

        for i, s in enumerate(sexes):
            pl.subplot(1, 3, i + 1)
            X = pl.array(
                sorted([[c.year, c.value] for c in ct.covariate_set.filter(region=region_dict[area], sex=s)]))
            if len(X) > 0:
                pl.plot(X[:,0], X[:,1], '.-')
                pl.ylabel(c.type)
                pl.xlabel('Time (Years)')
                pl.title('%s for %s in %s' % (c.type, s, c.region))

    response = view_utils.figure_data(format)
    
    return HttpResponse(response, view_utils.MIMETYPE[format])


def covariate_show(request, type, area, sex, format='png'):
    """ Serve a representation of the selected covariate

    Parameters:
    -----------
      type : str
        the covariate_type
      area : str
        either the country code or the gbd region
      sex : str
        the sex to display ('male', 'female', 'all', or '')
      format : str, optional
        the format to return the results in, may be one of the following:
        json, csv, png, pdf
    """
    ct = get_object_or_404(CovariateType, slug=type)
    fig_width = 6.
    fig_height = 4.5
    pl.figure(figsize=(fig_width, fig_height), dpi=100)

    if len(area) == 3:
        X = pl.array(sorted([[c.year, c.value] for c in ct.covariate_set.filter(iso3=area, sex=sex)]))
        if len(X) > 0:
            pl.plot(X[:,0], X[:,1], '.-')
            pl.ylabel(c.type)
            pl.xlabel('Time (Years)')
            pl.title('%s in %s' % (c.type, c.iso3))
    else:
        region_dict = {}
        for r in gbd_regions:
            region_dict[clean(r)] = r
        X = pl.array(sorted([[c.year, c.value] for c in ct.covariate_set.filter(region=region_dict[area], sex=sex)]))
        if len(X) > 0:
            pl.plot(X[:,0], X[:,1], '.-')
            pl.ylabel(c.type)
            pl.xlabel('Time (Years)')
            pl.title('%s in %s' % (c.type, c.region))

    response = view_utils.figure_data(format)
    
    return HttpResponse(response, view_utils.MIMETYPE[format])


class NewDataForm(forms.Form):
    type = forms.CharField(max_length=50)
    file  = forms.FileField()
    rescale = forms.BooleanField(initial=True, required=False, help_text='Transform data to have mean o and variance 1?')
    source = forms.CharField(required=False, widget=forms.Textarea(attrs={'rows':1, 'cols':90, 'wrap': 'off'}))
    yearStart = forms.IntegerField(required=False)
    yearEnd = forms.IntegerField(required=False)
    regionOnly = forms.BooleanField(initial=False, required=False, help_text='Does the covariate data file have region column and no iso3 column?')
        
    def clean_file(self):
        import csv
        data = [d for d in csv.DictReader(self.file)]
        return data
    
    def clean(self):
        data = self.cleaned_data.get('file')
        type_slug = self.cleaned_data.get('type')
        source = self.cleaned_data.get('source')
        year_start = self.cleaned_data.get('yearStart')
        year_end = self.cleaned_data.get('yearEnd')
        region_only = self.cleaned_data.get('regionOnly')

        if data and type_slug:
            if year_start != None:
                try:
                    year = int(year_start)
                except ValueError:
                    raise forms.ValidationError('Could not interpret YearStart')
        
            if year_end != None:     
                try:
                    year = int(year_end)
                except ValueError:
                    raise forms.ValidationError('Could not interpret YearEnd')
            try:
                cov_type = CovariateType.objects.get(slug=type_slug)
            except CovariateType.DoesNotExist:
                if source == '':
                    raise forms.ValidationError('Source is missing')
                if year_start == None:
                    raise forms.ValidationError('YearStart is missing')
                if year_end == None:
                    raise forms.ValidationError('YearEnd is missing')

            if data == None:
                raise forms.ValidationError('Data is missing')       

            # make an iso3 list
            iso3_data = [x[1:] for x in csv.reader(open(CSV_PATH + 'country_region.csv'))]
            iso3_list = []
            for r in iso3_data:
                iso3_list += r

            # make a region list
            region_list = [x[0] for x in csv.reader(open(CSV_PATH + 'country_region.csv'))]

            # make a region_country_dict
            region_country_dict = {}
            for x in csv.reader(open(CSV_PATH + 'country_region.csv')):
                region_country_dict[x[0]] = x[1:]
        
            for ii, d in enumerate(data):
                try:
                    d['value'] = float(d[type_slug])
                except KeyError:
                    raise forms.ValidationError('Could not find column %s (is it spelled correctly?)' % type_slug)
                except ValueError:
                    raise forms.ValidationError('Could not interpret value for %s in line %d' % (type_slug, ii+2))

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

                if region_only and not d.has_key('region'):
                    raise forms.ValidationError('Could not find column region (is it spelled correctly?)')

                if not d.has_key('iso3') and not d.has_key('region'):
                    raise forms.ValidationError('Could not find either column iso3 or column region (is it spelled correctly?)')

                if d.has_key('iso3') and not d['iso3'] in iso3_list:
                    raise forms.ValidationError('Could not interpret iso3 in line %d' % (ii+2))

                if d.has_key('region') and not d['region'] in region_list:
                    raise forms.ValidationError('Could not interpret region in line %d' % (ii+2))

                if d.has_key('iso3') and d.has_key('region') and d['iso3'] not in region_country_dict[d['region']]:
                    raise forms.ValidationError('The iso3 and the region are inconsistent in line %d' % (ii+2))
                    
                if d.has_key('age'):
                    try:
                        int(d['age'])
                    except ValueError:
                        raise forms.ValidationError('Could not interpret age in line %d' % (ii+2))

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
            cov_data = form.cleaned_data['file']
            cov_type, is_new = CovariateType.objects.get_or_create(slug=type_slug, defaults={'year_start': 0, 'year_end': 0})

            if is_new and request.POST.get('notes') == '':
                 return render_to_response('covariate_upload.html', {'form': form, 'error': 'Notes are missing.'})

            cov_type.uploader = str(request.user)

            if form.cleaned_data['source'] != '':
                cov_type.source = form.cleaned_data['source']                
                
            cov_type.last_modified_time = datetime.now()
                
            if request.POST.get('notes') != '':
                cov_type.description = request.POST.get('notes')

            if form.cleaned_data['yearStart'] != None:
                cov_type.year_start = form.cleaned_data['yearStart']
        
            if form.cleaned_data['yearEnd'] != None:
                cov_type.year_end = form.cleaned_data['yearEnd']

            cov_type.region_only = form.cleaned_data['regionOnly']
        
            cov_type.save()

            # make rates from rate_list
            vals = [d['value'] for d in cov_data]

            if form.cleaned_data['rescale']:
                shift = pl.mean(vals)
                scale = pl.std(vals)
            else:
                shift = 0.
                scale = 1.
            
            for d in cov_data:
                # if sex == '' add a covariate for male, female, and total
                if d['sex'] == '':
                    sex_list = ['male', 'female', 'total']
                else:
                    sex_list = [d['sex']]
                for sex in sex_list:
                    # add a data point, save it on the data list
                    if d.has_key('iso3'):
                        cov, is_new = Covariate.objects.get_or_create(type=cov_type,
                                                                      iso3=d['iso3'],
                                                                      year=d['year'],
                                                                      sex=sex,
                                                                      defaults={'value': 0.})
                        if d.has_key('region'):
                            cov.region = d['region']
                    else:
                        cov, is_new = Covariate.objects.get_or_create(type=cov_type,
                                                                      region=d['region'],
                                                                      year=d['year'],
                                                                      sex=sex,
                                                                      defaults={'value': 0.})
                    cov.value = (d['value'] - shift) / scale
                    
                    if d.has_key('age'):
                        cov.age = d['age']
                    
                    cov.save()

            return HttpResponseRedirect(reverse('gbd.covariate_data_server.views.covariate_type_list_show')) # Redirect after POST

    return render_to_response('covariate_upload.html', {'form': form, 'error': ''})

