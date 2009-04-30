import pylab as pl
import numpy as np
import pymc as mc
from pymc import gp

from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

import generic_disease_model as default_prob_model

MISSING = -99

from disease_json import *
from model_utils import *

def fit(disease_model, probabilistic_model=default_prob_model):
    """
    fit a single data_type from the model
    """
    dm = get_disease_model(disease_model)

    # filter out all data with type != data_type
    # dm.data = dm.filter_data(data_type=data_type)

    # store the probabilistic model code for future reference
    dm.append_model_source(probabilistic_model)
    dm.set_param_age_mesh([0.0, 0.5, 3.0, 10.0, 20.0, 30.0, 40.0,
                           50.0, 60.0, 70.0, 80.0, 90.0, 100.0])

    dm.set_estimate_age_mesh(range(MAX_AGE))

    probabilistic_model.initialize(dm)
    vars = probabilistic_model.setup(dm)
    map = probabilistic_model.map_fit(dm, vars)
    mcmc = probabilistic_model.mcmc_fit(dm, vars)

    url = post_disease_model(dm)
    print 'url for fit:\n\t%s' % url

    return {'vars': vars,
            'map': map,
            'mcmc': mcmc,
            'disease_model': dm}






def plot_disease_model(dm_json):
    dm = DiseaseJson(dm_json)
    # divide up disease_model data by data_type
    data_by_type = {}
    for d in dm.data:
        data_by_type[d['data_type']] = data_by_type.get(d['data_type'], []) + [d]

    types = set(data_by_type.keys()) | set(dm.params.get('map', {}).keys())
    cnt = max(1, len(types))
    cols = min(2, cnt)
    rows = int(np.ceil(float(cnt) / float(cols)))

    subplot_width = 6
    subplot_height = 4
    
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    for ii, type in enumerate(types):
        data = data_by_type.get(type, [])
        
        pl.subplot(rows, cols, ii + 1)
        plot_intervals(dm, data, fontsize=12, alpha=.5)
        plot_normal_approx(dm, type)
        plot_map_fit(dm, type)
        plot_mcmc_fit(dm, type)
        plot_truth(dm, type)
        plot_prior(dm, type)
        label_plot(dm, type, fontsize=10)
        
        # max_rate = np.max([.0001] + [d['value']*extract_units(d) for d in data] + list(dm.get_map(type)) + list(dm.get_mcmc('upper_ui', type)))
        if len(data) > -1:
            max_rate = np.max([.01] + [dm.value_per_1(d) for d in data] + list(dm.get_map(type)))
            ages = dm.get_estimate_age_mesh()
            xmin = ages[0]
            xmax = ages[-1]
            ymin = 0.
            ymax = 1.25*max_rate
            pl.axis([xmin, xmax, ymin, ymax])

def plot_intervals(dm, data, alpha=.75, color=(.0,.5,.0), text_color=(.0,.3,.0), fontsize=12):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of Data
    """
    for d in data:
        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE

        val = dm.value_per_1(d)
        se = dm.se_per_1(d)
        
        if se > 0.:  # don't draw error bars if standard error is zero or MISSING
            lower_ci = max(0., val - 1.98 * se)
            upper_ci = min(1., val + 1.98 * se)
            pl.plot([.5 * (d['age_start']+d['age_end']+1)]*2,
                    [lower_ci, upper_ci],
                    color=color, alpha=alpha, linewidth=1)
        pl.plot(np.array([d['age_start'], d['age_end']+1.]),
                np.array([val, val]),
                color=color, alpha=alpha, linewidth=5)
    
def plot_fit(dm, fit_name, data_type, **params):
    fit = dm.params.get(fit_name, {}).get(data_type)
    age = dm.get_estimate_age_mesh()
    if fit and age:
        pl.plot(age, fit, **params)

def plot_normal_approx(dm, type):
    plot_fit(dm, 'normal_approx', type, color='blue', alpha=.5)

def plot_truth(dm, type):
    plot_fit(dm, 'truth', type, linestyle='dotted', color='green', alpha=.95, linewidth=2)

def plot_map_fit(dm, type, **params):
    default_params = {'color': 'blue',
                      #'linestyle': 'dashed',
                      'linewidth': 2,
                      'alpha': .9,
                      }
    default_params.update(**params)
    plot_fit(dm, 'map', type, **default_params)

def plot_mcmc_fit(dm, type, color=(.2,.2,.2)):
    age = dm.get_estimate_age_mesh()
    lb = dm.get_mcmc('lower_ui', type)
    ub = dm.get_mcmc('upper_ui', type)

    if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
        x = np.concatenate((age, age[::-1]))
        y = np.concatenate((lb, ub[::-1]))
        pl.fill(x, y, facecolor='.2', edgecolor=color, alpha=.5)

    val = dm.get_mcmc('mean', type)

    if len(age) > 0 and len(age) == len(val):
        pl.plot(age, val, ':', color=color, linewidth=1, alpha=.75)

def plot_prior(dm, type):
    # show 'zero' priors
    for prior_str in dm.get_priors(type).split('\n'):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])

            pl.plot([age_start, age_end], [0, 0], color='red', linewidth=15, alpha=.75)

    # write out details of priors in a friendly font as well
    a0 = dm.get_estimate_age_mesh()[0]
    v0 = 0.
    pl.text(a0, v0, dm.get_priors(type).replace('\r\n', '\n'), color='black', family='monospace', fontsize=8, alpha=.75)
    

def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(dm, type, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('%s %s' % (type, dm.get_units(type)), **params)
    pl.title('%s; %s; %s; %s' % \
                 (dm.params['condition'],
                  dm.params['sex'], dm.params['region'],
                  dm.params['year']), **params)

