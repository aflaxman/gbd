""" All matplotlib-related methods for visualizing Disease Models

Useful High-level Methods::

    plot_disease_model(dm, [max_intervals])
    sparkplot_disease_model(dm, [max_intervals])

Useful Low-level Methods::

    plot_intervals(dm, data, [fontsize=12], [alpha=.5])
    plot_normal_approx(dm, type)
    plot_map_fit(dm, type)
    plot_mcmc_fit(dm, type)
    plot_truth(dm, type)
    plot_prior(dm, type)
    label_plot(dm, type, [fontsize=10])
    clear_plot()
"""

import pylab as pl
import numpy as np
import random

import dismod3
from disease_json import *


def plot_disease_model(dm_json, max_intervals=50):
    """Make a graphic representation of the disease model data and
    estimates provided

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that is to be plotted
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return
        
    # divide up disease_model data by data_type
    data_by_type = {}
    
    for d in dm.data:
        data_by_type[d['data_type']] = data_by_type.get(d['data_type'], []) + [d]

    types = set(data_by_type.keys()) | set(dm.params.get('initial_value', {}).keys())
    cnt = max(1, len(types))
    cols = min(2, cnt)
    rows = int(np.ceil(float(cnt) / float(cols)))

    subplot_width = 6
    subplot_height = 4
    
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    for ii, type in enumerate(types):
        data = data_by_type.get(type, [])
        if len(data) > max_intervals:
            data = random.sample(data, max_intervals)
        
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
            max_rate = np.max([.01] + [dm.value_per_1(d) for d in data] + list(dm.get_map(type))+ list(dm.get_mcmc('mean', type)))
            ages = dm.get_estimate_age_mesh()
            xmin = ages[0]
            xmax = ages[-1]
            ymin = 0.
            ymax = 1.25*max_rate
            pl.axis([xmin, xmax, ymin, ymax])

color_for = {
    'incidence data': 'cyan',
    'incidence': 'cyan',
    'prevalence data': 'blue',
    'prevalence': 'blue',
    'remission data': 'green',
    'remission': 'green',
    'case-fatality data': 'red',
    'case-fatality': 'red',
    'all-cause mortality data': 'black',
    'all-cause mortality': 'black',
    'duration data': 'orange',
    'duration': 'orange',
    }

def sparkplot_boxes(dm_json):
    """ Find pixels for all boxes in the sparkplot lattice below."""
    rows = len(dismod3.gbd_regions) + 1
    col_list = [[1995, 'male'],
                [2005, 'male'],
                [1995, 'female'],
                [2005, 'female'],
                ]
    cols = len(col_list)

    subplot_width = 1.
    subplot_height = 2./3.
    fig_width = subplot_width*cols
    fig_height = subplot_height*rows

    subplot_px = {}
    for ii, region in enumerate(dismod3.gbd_regions):
        for jj, [year, sex] in enumerate(col_list):
            subplot_px[dismod3.gbd_key_for('all', region, year, sex)] = \
                ', '.join([str(int(100 * (jj) * subplot_width)),
                           str(int(100 * (rows - ii - 1) * subplot_height)),
                           str(int(100 * (jj + 1) * subplot_width)),
                           str(int(100 * (rows - ii) * subplot_height))])

    ii += 1
    subplot_px[dismod3.gbd_key_for('all', 'world', 'total', 'total')] = \
                                          ', '.join(['0',
                                                     str(int(100 * (rows - ii - 1) * subplot_height)),
                                                     str(int(100 * (jj + 1) * subplot_width)),
                                                     str(int(100 * (rows - ii) * subplot_height))])

    return subplot_px

def sparkplot_disease_model(dm_json, max_intervals=50):
    """ Make a lattice of sparkplots for the disease_model, with rows
    corresponding to regions, and columns corresponding to (year,sex)
    pairs.

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that is to be plotted
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return

    class GBDDataHash:
        """ Store and serve data grouped by type, region, year, and sex
        """
        def __init__(self, data):
            self.data = data

        def get(self, type, region, year, sex):
            """ Provide a way to get desired data

            Parameters
            ----------
            type : str, one of the following types
                     'incidence data', 'prevalence data', 'remission data',
                     'case-fatality data', 'all-cause mortality data', 'duration data'
            region : str, one of the 21 gbd regions or 'World'
            year : int, one of 1995, 2005
            sex : str, one of 'male', 'female', 'total'

            Notes
            -----
            TODO:  speed this up by dividing up data once and caching that
            """
            d_list = []
            for d in dm.data:
                if type == '*' or d['data_type'] == type:
                    if region == '*' or d['gbd_region'] == region:
                        if year == '*' or (year == 1995 and d['year_start'] < 2000) \
                                or (year == 2005 and d['year_end'] >= 2000):
                            if sex == '*' or d['sex'] == sex:
                                d_list.append(d)
            return d_list
        
        
    # divide up disease_model data by data_type, region, year, sex
    data_hash = GBDDataHash(dm.data)

    # divide up canvas to fit all the sparklines
    rows = len(dismod3.gbd_regions)+1
    col_list = [[1995, 'male'],
                [2005, 'male'],
                [1995, 'female'],
                [2005, 'female'],
                ]
    cols = len(col_list)

    subplot_width = 1. * .5
    subplot_height = 2./3. * .5
    fig_width = subplot_width*cols
    fig_height = subplot_height*rows

    fig = pl.figure(figsize=(fig_width,fig_height), dpi=100)
    pl.clf()

    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    ymax = 1. #dm.get_ymax()
    
    for ii, region in enumerate(dismod3.gbd_regions):
        for jj, [year, sex] in enumerate(col_list):
            fig.add_axes([jj*subplot_width / fig_width,
                          ii*subplot_height / fig_height,
                          subplot_width / fig_width,
                          subplot_height / fig_height],
                         frameon=False)
            # plot intervals and map_fit for each data type in a different color
            for type in dismod3.data_types:
                data = data_hash.get(type, region, year, sex) + data_hash.get(type, region, year, 'total')
                if len(data) > max_intervals:
                    data = random.sample(data, max_intervals)
        
                plot_intervals(dm, data, alpha=.5, color=color_for[type])
                type = type.replace(' data', '')
                plot_map_fit(dm, dismod3.gbd_key_for(type, region, year, sex),
                             linestyle='-', alpha=.5, color=color_for[type])
            pl.xticks([])
            pl.yticks([])
            pl.axis([xmin, xmax, ymin, ymax])

    ii += 1
    fig.add_axes([0,
                  ii*subplot_height / fig_height,
                  fig_width,
                  subplot_height / fig_height],
                 frameon=False)
    for type in dismod3.data_types:
        type = type.replace(' data', '')
        plot_map_fit(dm, dismod3.gbd_key_for(type, 'world', 'total', 'total'),
                     linestyle='-', alpha=.5, color=color_for[type])
        pl.xticks([])
        pl.yticks([])
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
        if val == MISSING:
            continue
        
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
    
def plot_fit(dm, fit_name, key, **params):
    fit = dm.params.get(fit_name, {}).get(key)
    age = dm.get_estimate_age_mesh()
    if fit and age:
        pl.plot(age, fit, **params)

def plot_normal_approx(dm, type):
    plot_fit(dm, 'normal_approx', type, color='blue', alpha=.5)

def plot_truth(dm, type):
    plot_fit(dm, 'truth', type, linestyle=':', color='green', alpha=.95, linewidth=2)

def plot_map_fit(dm, type, **params):
    default_params = {'color': 'blue',
                      'linestyle': ':',
                      'linewidth': 2,
                      'alpha': .9,
                      'label': 'Max-liklihood',
                      }
    default_params.update(**params)
    plot_fit(dm, 'map', type, **default_params)

def plot_mcmc_fit(dm, type, color=(.2,.2,.2)):
    age = dm.get_estimate_age_mesh()
    param_mesh = dm.get_param_age_mesh()
    
    lb = dm.get_mcmc('lower_ui', type)
    ub = dm.get_mcmc('upper_ui', type)

    if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
        lb = lb[param_mesh]
        ub = ub[param_mesh]

        x = np.concatenate((param_mesh, param_mesh[::-1]))
        y = np.concatenate((lb, ub[::-1]))
        pl.fill(x, y, facecolor='.2', edgecolor=color, alpha=.5, label='MCMC 95% UI')

    val = dm.get_mcmc('median', type)

    if len(age) > 0 and len(age) == len(val):
        val = val[param_mesh]
        pl.plot(param_mesh, val, color=color, linewidth=4, alpha=.75, label='MCMC Median')

def plot_prior(dm, type):
    # show 'zero' priors
    for prior_str in dm.get_priors(type).split('\n'):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2])

            pl.plot([age_start, age_end], [0, 0], color='red', linewidth=15, alpha=.75)

    # write out details of priors in a friendly font as well
    if len(dm.get_estimate_age_mesh()) > 0:
        a0 = dm.get_estimate_age_mesh()[0]
        v0 = 0.
        pl.text(a0, v0, ' Priors:\n' + dm.get_priors(type).replace('\r\n', '\n'), color='black', family='monospace', fontsize=8, alpha=.75)
    

def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(dm, type, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('%s %s' % (type, dm.get_units(type)), **params)
    pl.title('%d: %s; %s; %s; %s' % \
                 (dm.params['id'], dm.params['condition'],
                  dm.params['sex'], dm.params['region'],
                  dm.params['year']), **params)
    pl.legend()
