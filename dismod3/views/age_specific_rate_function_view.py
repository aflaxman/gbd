from django.shortcuts import render_to_response, get_object_or_404
from django.http import *
from django.core.urlresolvers import reverse
from django.utils.translation import ugettext as _
from django import forms

import pymc.gp as gp
import numpy as np
import pylab as pl

from dismod3.models import *
import view_utils

def age_specific_rate_function_show(request, id_str):
    asrfs = view_utils.id_str_to_objects(id_str, AgeSpecificRateFunction)
    return render_to_response('age_specific_rate_function/show.html',
                              view_utils.template_params(asrfs[0], asrfs=asrfs, id_str=id_str))

def age_specific_rate_function_redirect(request, id_str, action):
    asrfs = view_utils.id_str_to_objects(id_str, AgeSpecificRateFunction)

    if action == 'edit':
        url = '/admin/dismod3/agespecificratefunction/%d' % asrfs[0].id
    elif action in view_utils.command_list['move']:
        url = reverse('dismod3.views.age_specific_rate_function_show', args=(asrfs[0].id+view_utils.id_delta[action],))
    elif action in view_utils.command_list['sex']:
        url = AgeSpecificRateFunction.objects.filter(disease=asrfs[0].disease, region=asrfs[0].region, rate_type=asrfs[0].rate_type, sex=action)[0].get_absolute_url()
    elif action in view_utils.command_list['format']:
        url = '%s.%s' % (asrfs[0].get_absolute_url(), action)
    else:
        raise Http404
    
    return HttpResponseRedirect(url)


class NotesForm(forms.Form):
    notes = forms.CharField(required=False)

def age_specific_rate_function_clone(request, id):
    asrf = get_object_or_404(AgeSpecificRateFunction, id=id)

    # http customs dictate using POSTs for any interaction which will
    # change the database
    if request.method == 'POST':
        form = NotesForm(request.POST)
        if form.is_valid():
            new_asrf = asrf.clone(**form.cleaned_data)
            return HttpResponseRedirect(new_asrf.get_absolute_url())
    else:
        form = NotesForm()

    return render_to_response('age_specific_rate_function/clone.html', {'rf': asrf, 'form': form})

def url_params_to_dict(string):
    return {}


def asrf_posterior_predictive_check_scatter(request, id, format):
    return asrf_posterior_predictive_check(request, id, format, style='scatter')

def asrf_posterior_predictive_check_intervals(request, id, format):
    return asrf_posterior_predictive_check(request, id, format, style='intervals')

def same_side_of_xy(x, y):
    x0, x1 = x
    y0, y1 = y
    vals = np.array([x0-y0, x0-y1, x1-y0, x1-y1])
    return (vals >= 0.).all() or (vals <= 0.).all()

def asrf_posterior_predictive_check(request, id, format, style):
    """
    generate a posterior predictive check of the model's
    'goodness-of-fit', which is to say a scatterplot of the observed
    rates versus the rates predicted by the model parameters
    """
    rf = get_object_or_404(AgeSpecificRateFunction, id=id)


    params = {'fontsize': 7}
    #params.update(url_params_to_dict(param_str))

    fig= view_utils.clear_plot(width=4,height=3)
    ax = fig.add_subplot(111)
    
    
    if rf.rates.count() == 0:
        pl.figtext(.5, .5, 'no rates found')
    
        return HttpResponse(view_utils.figure_data(format),
                            view_utils.MIMETYPE[format])
    
    observed_rates = np.array([ float(rate.numerator)/float(rate.denominator) for rate in rf.rates.all() ])
    observed_cls = np.transpose([ rate.ci() for rate in rf.rates.all() ])

    predicted_rates = np.array([ probabilistic_utils.predict_rate_from_asrf(rf, rate) for rate in rf.rates.all() ])
    predicted_cls = np.array([ [ probabilistic_utils.predict_rate_from_asrf(rf, rate, 'mcmc_lower_cl') for rate in rf.rates.all() ],
                                [ probabilistic_utils.predict_rate_from_asrf(rf, rate, 'mcmc_upper_cl') for rate in rf.rates.all() ] ])

    max_x = np.max(observed_cls[1,:])
    max_y = np.max(predicted_cls[1,:])
    max_t = max(max_x, max_y, probabilistic_utils.NEARLY_ZERO)

    observed_errs = np.abs(observed_cls - observed_rates)
    predicted_errs = np.abs(predicted_cls - predicted_rates)

    if style == 'scatter':
        pl.plot([probabilistic_utils.NEARLY_ZERO,1.], [probabilistic_utils.NEARLY_ZERO,1.], linestyle='dashed', linewidth=2, color='black', alpha=.75)
        pl.errorbar(x=observed_rates + np.random.rand(len(observed_rates)) * max_t / 250., xerr=observed_errs,
                    y=predicted_rates + np.random.rand(len(predicted_rates)) * max_t / 250., yerr=predicted_errs, fmt='bo')

        from matplotlib.patches import Ellipse
        
        for ii in range(len(observed_rates)):
            e = Ellipse(xy=[np.mean(observed_cls[:,ii]),np.mean(predicted_cls[:,ii])],
                        width=(observed_cls[1,ii] - observed_cls[0,ii]),
                        height=(predicted_cls[1,ii] - predicted_cls[0,ii]),
                        alpha=.5)
            if same_side_of_xy(observed_cls[:,ii], predicted_cls[:,ii]):
                e.set_facecolor((1.,0.,0.))
            ax.add_artist(e)
            
            
        
        pl.axis([0., max_x + 0.001, 0., max_y + 0.001])

        tick_list, tick_objs = pl.xticks()
        pl.xticks([0,tick_list[-1]], **params)
        pl.xlabel('Observed Rates', **params)

        tick_list, tick_objs = pl.yticks()
        pl.yticks([0,tick_list[-1]], **params)
        pl.ylabel('Predicted Rates', **params)
    elif style == 'intervals':
        rate_list = rf.rates.all()

        pl.subplot(2,1,1)
        plot_intervals(rf, rate_list, color='green', **params)
        pl.axis([0., 100., 0., max_t])
        pl.xticks([])
        pl.yticks([])
        pl.xlabel('')
        pl.ylabel('')

        for x,r in zip(predicted_rates, rate_list):
            r.numerator = x * r.denominator

        pl.subplot(2,1,2)
        plot_intervals(rf, rate_list, color='red', **params)
        pl.axis([0., 100., 0., max_t])
        pl.xticks([])
        pl.yticks([])
        pl.ylabel('')
        pl.title('')
        
        pl.subplot(2,1,1)
    else:
        raise Exception('ERROR: plot style "%s" unrecognized' % style)
    pl.title('Predictive Check for %s' % rf, **params)
    
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])
    
def age_specific_rate_function_plot(request, id_str, format):
    asrfs = view_utils.id_str_to_objects(id_str, AgeSpecificRateFunction)

    # handle json & csv formats, which are not technically plots
    if format in ['json', 'csv']:
        if format == 'json':
            data_str = json.dumps([[rf.id, rf.fit] for rf in asrfs])
        elif format == 'csv':
            headings = {}
            rows = {}
            data_str = ''

            for rf in asrfs:
                headings[rf] = ['Age (years)', 'MAP Rate (per 1.0)']
                rows[rf] = [[a, p] for a,p in zip(rf.fit['out_age_mesh'], rf.fit['map'])]
                data_str += view_utils.csv_str(headings[rf], rows[rf])
        return HttpResponse(data_str, view_utils.MIMETYPE[format])

    cnt = asrfs.count()
    cols = (cnt / 10) + 1
    rows = (cnt / cols)
#    rows = 1

    subplot_width = 6
    subplot_height = 4
    
    view_utils.clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    for ii, rf in enumerate(asrfs):
        pl.subplot(rows,cols,ii+1)
        plot_intervals(rf, rf.rates.all(), fontsize=12)
        #plot_map_fit(rf)
        plot_mcmc_fit(rf)
        plot_prior(rf)

        max_rate = np.max([.0001] + [r.rate for r in rf.rates.all()])
        pl.axis([0, 100, 0, 1.25*max_rate])
        if ii % cols != 0:
            pl.ylabel('')
        #pl.yticks([])
        if (ii + cols) < cnt:
            pl.xlabel('')
        #pl.xticks([])
    
    
    return HttpResponse(view_utils.figure_data(format),
                        view_utils.MIMETYPE[format])


class ASRFCreationForm(forms.Form):
    disease = forms.ModelChoiceField(Disease.objects.all())
    region = forms.ModelChoiceField(Region.objects.all(), required=False)
    rate_type = forms.ChoiceField(fields.ALL_OPTION + fields.RATE_TYPE_CHOICES, required=False)
    sex = forms.ChoiceField(fields.ALL_OPTION + fields.SEX_CHOICES)
    notes = forms.CharField(required=False, widget=forms.widgets.Textarea())

def age_specific_rate_function_index(request):
    if request.method == 'POST': # If the form has been submitted...
        form = ASRFCreationForm(request.POST) # A form bound to the POST data
        if form.is_valid(): # All validation rules pass
            asrfs = age_specific_rate_function.create_multiple(**form.cleaned_data)
            return HttpResponseRedirect('/age_specific_rate_function/%s' % view_utils.objects_to_id_str(asrfs)) # Redirect after POST
    else:
        form = ASRFCreationForm()

    return render_to_response('age_specific_rate_function/index.html', {'form': form})

def plot_intervals(rf, rate_list, alpha=.75, color=(.0,.5,.0), text_color=(.0,.3,.0), fontsize=12):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of RateData that will be used in fitting this
    age-specific rate function
    """
    for r in rate_list:
        rate_val = float(r.numerator)/float(r.denominator)
        x_jitter = 0.*np.random.rand()
        y_jitter = 0.
        pl.plot(.5*np.array([(r.age_start+r.age_end+1)]*2)+x_jitter,
                r.ci(),
                color=color, alpha=alpha, linewidth=1)
        pl.plot(np.array([r.age_start, r.age_end+1.]),
                np.array([rate_val,rate_val]),
                color=color, alpha=alpha, linewidth=5,
                )
        # pl.text(r.age_end+x_jitter, rate_val+y_jitter,
        #         "n=%d" % r.rate_denominator,
        #         color=text_color, alpha=alpha, fontsize=8)
    view_utils.label_plot(rf, fontsize=fontsize)
    
def plot_fit(rf, fit_name, color='blue', linestyle='solid'):
    try:
        rate = rf.fit[fit_name]
        pl.plot(rf.fit['out_age_mesh'], rate, color=color, linewidth=2, linestyle=linestyle,
                alpha=.75, label=fit_name)
    except (KeyError, ValueError):
        pl.figtext(0.4,0.2, 'No %s data Found' % fit_name)

def plot_normal_approx(rf):
    plot_fit(rf, 'normal_approx', 'blue')

def plot_map_fit(rf):
    plot_fit(rf, 'map', 'black', 'dashed')

def plot_mcmc_fit(rf, detailed_legend=False, color='black'):
    try:
        x = np.concatenate((rf.fit['out_age_mesh'], rf.fit['out_age_mesh'][::-1]))
        y = np.concatenate((rf.fit['mcmc_lower_cl'], rf.fit['mcmc_upper_cl'][::-1]))

        pl.fill(x, y, facecolor='.2', edgecolor=color, alpha=.5)

        mcmc_mean = rf.fit['mcmc_mean']

        if detailed_legend:
            label = str(rf.region)
            color = np.random.rand(3)
        else:
            label = 'MCMC Fit'
            color = color

        pl.plot(rf.fit['out_age_mesh'], mcmc_mean, color=color, linewidth=3, alpha=.75, label=label)
    except (KeyError, ValueError):
        pl.figtext(0.4,0.4, 'No MCMC Fit Found')

def plot_prior(rf):
    pl.plot([0,5],[0,0], color='red', linewidth=15, alpha=.75)
