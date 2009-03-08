from django.shortcuts import render_to_response, get_object_or_404
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms
from django.core.paginator import Paginator, InvalidPage, EmptyPage

from dismod3.models import *
import view_utils

class DiseaseModelForm(forms.Form):
    notes = forms.CharField(required=False)

def disease_model_create(request):
    # http customs dictate using POSTs for any interaction which will
    # change the database
    if request.method == 'POST':
        form = DiseaseModelForm(request.POST)
        if form.is_valid():
            # TODO: create new disease model
            return HttpResponseRedirect(new_disease_model.get_absolute_url())
    else:
        form = DiseaseModelForm()

    return render_to_response('disease_model/create.html', {'form': form})
    
def disease_model_show(request, id, format='html'):
    dm = get_object_or_404(DiseaseModel, id=id)

    if format == 'html':
        return render_to_response('disease_model/show.html',
                                  view_utils.template_params(dm, id_str=dm.get_asrf_id_str(), asrfs=dm.rates.all()))

    # TODO: handle json & csv formats
    if format in ['json', 'csv']:
        if format == 'json':
            data_str = json.dumps([[rf.id, rf.fit] for rf in dm.rates.all()])
        elif format == 'csv':
            headings = {}
            rows = {}
            data_str = ''

            for rf in dm.rates.all():
                headings[rf] = ['Age (years)', 'MAP Rate (per 1.0)']
                rows[rf] = [[a, p] for a,p in zip(rf.fit['out_age_mesh'], rf.fit['map'])]
                data_str += view_utils.csv_str(headings[rf], rows[rf])
        return HttpResponse(data_str, view_utils.MIMETYPE[format])

    else:
        raise Http404

class DiseaseModelCreationForm(forms.Form):
    disease = forms.ModelChoiceField(Disease.objects.all())
    region = forms.ModelChoiceField(Region.objects.all(), required=False)
    rate_type = forms.ChoiceField(fields.ALL_OPTION + fields.RATE_TYPE_CHOICES, required=False)
    sex = forms.ChoiceField(fields.ALL_OPTION + fields.SEX_CHOICES)
    notes = forms.CharField(required=False, widget=forms.widgets.Textarea())

def disease_model_index(request):
    if request.method == 'POST': # If the form has been submitted...
        form = DiseaseModelCreationForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            dm = DiseaseModel(**form.cleaned_data)
            return HttpResponseRedirect(dm.get_absolute_url()) # Redirect after POST
    else:
        form = DiseaseModelCreationForm()

    dm_list = DiseaseModel.objects.all().order_by('-id')
    paginator = Paginator(dm_list, per_page=10)
    
    # Make sure page request is an int. If not, deliver first page.
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1
                            
    try:
        dms = paginator.page(page)
    except (EmptyPage, InvalidPage):
        dms = paginator.page(paginator.num_pages)

    return render_to_response('disease_model/index.html', {'form': form, 'disease_models': dms})
