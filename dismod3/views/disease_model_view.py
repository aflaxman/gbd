from django.shortcuts import render_to_response, get_object_or_404
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

from dismod3.models import *
import view_utils
from age_specific_rate_function_view import plot_map_fit, plot_mcmc_fit, plot_truth, plot_prior

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

    dms = DiseaseModel.objects.all().order_by('-id')
    paginated_models = view_utils.paginated_models(request, dms)

    return render_to_response('disease_model/index.html', {'form': form, 'paginated_models': paginated_models})

def disease_model_sparkplot(request, id, format):
    dm = get_object_or_404(DiseaseModel, id=id)

    width = 1
    height = .5
    
    fig = view_utils.clear_plot(width,height)

    ax = None
    for ii, rf in enumerate(dm.rates.all()):
        ax = pl.axes([ii/4., 0., 1., (ii+1)/4.], frameon=False)
        pl.subplot(4,1,ii+1)
        plot_map_fit(rf)
        #plot_mcmc_fit(rf)
        plot_truth(rf)
        #plot_prior(rf)
        pl.xticks([])
        pl.yticks([])
        #pl.delaxes()
    pl.subplots_adjust(left=0, bottom=0, right=1, top=1,
                    wspace=0, hspace=0)

    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])
    
