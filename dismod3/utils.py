import inspect

import pylab as pl
import numpy as np
import pymc as mc
from pymc import gp

import twill.commands as twc
import simplejson as json

from dismod3.settings import *
from model_utils import *
from bayesian_models import probabilistic_utils
from bayesian_models.probabilistic_utils import trim, uninformative_prior_gp, NEARLY_ZERO, MAX_AGE

import beta_binomial_model as probabilistic_model

MISSING = -99


def get_disease_model(disease_model_id):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_LOGIN_URL)
    twc.fv('1', 'username', DISMOD_USERNAME)
    twc.fv('1', 'password', DISMOD_PASSWORD)
    twc.submit()
    twc.url('accounts/profile')

    twc.go(DISMOD_DOWNLOAD_URL % disease_model_id)

    return json.loads(twc.show())

def post_disease_model(disease_model):
    """
    fetch specificed disease model data from
    dismod server given in settings.py
    """
    
    twc.go(DISMOD_UPLOAD_URL)
    twc.fv('1', 'model_json', json.dumps(disease_model))
    twc.submit()

    return twc.browser.get_url()


def fit(disease_model, data_type='prevalence data'):
    """
    fit a single data_type from the model
    """
    dm = get_disease_model(disease_model)

    # filter out all data with type != data_type
    dm['data'] = [d for d in dm['data'] if d['data_type'] == data_type]

    # store the probabilistic model code for future reference
    dm['params'].update(
        model_source=inspect.getsource(probabilistic_model),
        age_mesh=[0.0, 0.5, 3.0, 10.0, 20.0, 30.0, 40.0,
                  50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
        out_age_mesh=range(MAX_AGE),
        map={},
        units={})

    dm['params']['units'][data_type] = '(per 1)'
    # do normal approximation first, to generate a good starting point
    fit_normal_approx(dm, data_type)
    
    # define the model variables
    model = probabilistic_model.setup(dm, data_type)

    map = mc.MAP(model)
    map.fit(method='fmin_powell', iterlim=500, tol=.001, verbose=1)
    dm['params']['map'][data_type] = list(model.rate.value)
    
    mcmc = mc.MCMC(model)
    mcmc.sample(iter=4000, burn=1000, thin=2, verbose=1)
    store_mcmc_fit(dm, mcmc, data_type)

    url = post_disease_model(dm)
    print 'url for fit: %s' % url

    model.dm = dm
    return model

def store_mcmc_fit(dm, mcmc, type):
    p = dm['params']
    p['mcmc_lower_ui'] = {}
    p['mcmc_median'] = {}
    p['mcmc_upper_ui'] = {}
    p['mcmc_mean'] = {}
    
    rate = mcmc.rate.trace()
    trace_len = len(rate)
    age_len = len(p['out_age_mesh'])
    
    sr = []
    for ii in xrange(age_len):
        sr.append(sorted(rate[:,ii]))
    p['mcmc_lower_ui'][type] = [sr[ii][int(.025*trace_len)] for ii in xrange(age_len)]
    p['mcmc_median'][type] = [sr[ii][int(.5*trace_len)] for ii in xrange(age_len)]
    p['mcmc_upper_ui'][type] = [sr[ii][int(.975*trace_len)] for ii in xrange(age_len)]
    p['mcmc_mean'][type] = list(np.mean(rate, 0))







def plot_disease_model(dm):
    
    # divide up disease_model data by data_type
    #import pdb; pdb.set_trace()
    data_by_type = {}
    for d in dm['data']:
        data_by_type[d['data_type']] = data_by_type.get(d['data_type'], []) + [d]

    types = set(data_by_type.keys()) | set(dm['params'].get('map', {}).keys())
    cnt = max(1, len(types))
    cols = min(2, cnt)
    rows = int(np.ceil(float(cnt) / float(cols)))

    subplot_width = 6
    subplot_height = 4
    
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    for ii, type in enumerate(types):
        data = data_by_type.get(type, [])
        
        pl.subplot(rows, cols, ii + 1)
        plot_intervals(data, fontsize=12, alpha=.5)
        plot_normal_approx(dm, type)
        plot_map_fit(dm, type)
        plot_mcmc_fit(dm, type)
        plot_truth(dm, type)
        #plot_prior(dm, type)
        label_plot(dm, type, fontsize=10)
        
        max_rate = np.max([.0001] + [d['value']*extract_units(d) for d in data])
        xmin = 0.
        xmax = 85.
        ymin = 0.
        ymax = 1.25*max_rate
        pl.axis([xmin, xmax, ymin, ymax])

        if ii % cols != 0:
            pl.ylabel('')
        if (ii + cols) < cnt:
            pl.xlabel('')

def plot_intervals(data, alpha=.75, color=(.0,.5,.0), text_color=(.0,.3,.0), fontsize=12):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of Data
    """
    for d in data:
        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE

        scale = extract_units(d)
        val = d['value'] * scale
        lower_ci = max(0., val - 1.98 * d['standard_error'] * scale)
        upper_ci = min(1., val + 1.98 * d['standard_error'] * scale)
        pl.plot([.5 * (d['age_start']+d['age_end']+1)]*2,
                [lower_ci, upper_ci],
                color=color, alpha=alpha, linewidth=1)
        pl.plot(np.array([d['age_start'], d['age_end']+1.]),
                np.array([val, val]),
                color=color, alpha=alpha, linewidth=5,
                )
    
def plot_fit(dm, fit_name, data_type, **params):
    fit = dm['params'].get(fit_name, {}).get(data_type)
    age = dm['params'].get('out_age_mesh')
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
    age = dm['params'].get('out_age_mesh')
    lb = dm['params'].get('mcmc_lower_ui', {}).get(type)
    ub = dm['params'].get('mcmc_upper_ui', {}).get(type)

    if age and lb and ub:
        x = np.concatenate((age, age[::-1]))
        y = np.concatenate((lb, ub[::-1]))
        pl.fill(x, y, facecolor='.2', edgecolor=color, alpha=.5)

    val = dm['params'].get('mcmc_mean', {}).get(type)

    if age and val:
        pl.plot(age, val, ':', color=color, linewidth=1, alpha=.75)


def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(dm, type, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('%s %s' % (type, dm['params'].get('units', {}).get(type, '')), **params)
    pl.title('%s; %s; %s; %s' % \
                 (dm['params']['condition'],
                  dm['params']['sex'], dm['params']['region'],
                  dm['params']['year']), **params)

