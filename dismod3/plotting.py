""" All matplotlib-related methods for visualizing Disease Models

Useful High-level Methods::

    plot_disease_model(dm, [max_intervals])
    sparkplot_disease_model(dm, [max_intervals])
    map_plot_int(region_value_dict)
    map_plot_float(region_value_dict)

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

import math
import copy
#import random
import pylab as pl
import numpy as np
from time import strftime
from operator import itemgetter

import dismod3
from dismod3.utils import clean, rate_for_range
from disease_json import *
from dismod3 import settings

color_for = {
    'incidence data': 'magenta',
    'incidence': 'magenta',
    'prevalence data': 'blue',
    'prevalence': 'blue',
    'remission data': 'green',
    'remission': 'green',
    'excess-mortality data': '#a04040',
    'excess-mortality': '#a00000',
    'all-cause mortality data': 'black',
    'all-cause mortality': 'black',
    'duration data': 'orange',
    'duration': 'orange',
    'relative-risk data': '#ff00ff',
    'relative-risk': '#990099',
    'yld': 'black',
    }

def darken(c):
    return (c[0]*.5, c[1]*.5, c[2]*.5)

default_max_for = {
    'incidence': .0005,
    'prevalence': .0001,
    'remission': .1,
    'excess-mortality': .01,
    'mortality': .1,
    'duration': 5,
    'relative-risk': 1,
    'incidence_x_duration': .0001,
    }

def prettify(str):
    """ Turn underscores into spaces"""
    return str.replace('_', ' ')

def overlay_plot_disease_model(dm_json_list, keys, max_intervals=100, defaults={}):
    """ Make a graphic representation of the disease model estimates

    Parameters
    ----------
    dm_json_list : list of strs or DiseaseJson objects
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list
      the keys to include
    """
    dm_list = []
    
    for dm_json in dm_json_list:
        if isinstance(dm_json, DiseaseJson):
            dm_list.append(dm_json)
        else:
            try:
                dm_list.append(DiseaseJson(dm_json))
            except ValueError:
                print 'ERROR: dm_json is not a DiseaseJson object or json string'
                return

    clear_plot(width=6, height=4)
    for ii, dm in enumerate(dm_list):
        data_hash = GBDDataHash(dm.data)

        keys = [k for k in keys if k.split('+')[0] in ['prevalence', 'incidence', 'remission', 'excess-mortality',
                                                       'mortality', 'duration', 'relative-risk', 'incidence_x_duration']]

        for k in sorted(keys, key=lambda k: np.max(list(dm.get_map(k)) + [0]), reverse=True):
            type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)
            data_type = type + ' data'

            # plot the data rectangles for these keys (only on first dm of list)
            if ii == 0:
                data = data_hash.get(data_type, region, year, sex)
                #if len(data) > max_intervals:
                #    data = random.sample(data, max_intervals)
                plot_intervals(dm, data, color=color_for.get(data_type, 'black'))

            # plot the map fit
            #plot_map_fit(dm, k, linestyle='-', linewidth=3, zorder=-5,
            #             color=pl.cm.spectral(ii / float(len(dm_list))),
            #             label=('%d: ' % dm.id) + k.split('+')[0])
            color = pl.cm.spectral(ii / float(len(dm_list)))
            plot_mcmc_fit(dm, k,
                          color=color,
                          show_data_ui=False)

            # plot the empirical prior, if there is any
            plot_empirical_prior(dm, k, linestyle='--', linewidth=3, zorder=-5,
                                 color=color,
                                 label=('%d: ' % dm.id) + k.split('+')[0])
        info_str = '%d: ' % dm.id
        for gof in ['aic', 'bic', 'dic']:
            val = dm.get_mcmc(gof, k)
            if val:
                info_str += '$%s = %.0f$;\t' % (gof.upper(), val)
        pl.text(0., 0., info_str + '\n'*(len(dm_list) - (ii+1)), fontsize=12, va='bottom', ha='left', color=darken(color))

    l,r,b,t, = pl.axis()
    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    try:
        ymax = float(defaults.get('ymax'))
    except (TypeError, ValueError):
        ymax = t
    pl.axis([xmin, xmax, ymin, ymax])

    # if this is a plot of all-cause mortality, make the y-axis log scale
    if dm.params['condition'] == 'all-cause_mortality':
        type='all-cause mortality data'
        data = data_hash.get(type, region, year, sex) + data_hash.get(type, region, year, 'total')

        plot_intervals(dm, data, color=color_for.get(data_type, 'black'))
        pl.plot(dm.get_estimate_age_mesh(), dm.mortality(dismod3.gbd_key_for(type, region, year, sex), data),
                alpha=.5, linestyle='-', color=color_for.get(type, 'black'), label='all-cause mortality')

        pl.semilogy([0.], [0.])
        pl.axis([xmin, xmax, 1.e-4, 1.])

    dm.params['id'] = '/'.join(str(dm.id) for dm in dm_list)
    dm.params['sex'] = sex
    dm.params['region'] = region
    dm.params['year'] = int(year)
    
    label_plot(dm, type, fontsize=10)

    pl.xticks(size='small')
    pl.yticks(size='small')
    pl.subplots_adjust(.15, .1, .95, .95)
            

    leg = pl.legend()
    try:
        # the matplotlib.patches.Rectangle instance surrounding the legend
        frame  = leg.get_frame()  
        frame.set_alpha(.2)    # set the frame face color to light gray
        frame.set_edgecolor('white')    # set the frame face color to light gray
            
        # matplotlib.text.Text instances
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
    except:
        pass

def bar_plot_disease_model(dm_json, keys, max_intervals=50):
    """ Make a barplot representation of the disease model data and
    estimates provided

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list
      the keys to include
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return
    data_hash = GBDDataHash(dm.data)

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] in ['incidence', 'prevalence', 'remission', 'excess-mortality']]

    rows = len(keys)

    subplot_width = 5
    subplot_height = 2
    clear_plot(width=subplot_width,height=subplot_height*rows)
    pl.subplots_adjust(left=.1)



    for i, k in enumerate(keys):
        pl.subplot(rows, 1, i+1)
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)

        if i+1 == 1:
            pl.title('%s: %s; %s; %s; %s' % \
                     (dm.params['id'], prettify(dm.params['condition']),
                      sex, region, year))
        

        data_type = clean(type) + ' data'
        data = data_hash.get(data_type, region, year, sex)
        #if len(data) > max_intervals:
        #    data = random.sample(data, max_intervals)
        plot_intervals(dm, data, color=color_for.get(data_type, 'black'), alpha=.2)
        
        plot_truth(dm, k, color=color_for.get(type, 'black'))
        plot_empirical_prior(dm, k, color=color_for.get(type, 'black'))  # TODO: make this less prominent
        # TODO: plot mcmc_mean over bars

        # plot level bars
        params = {}
        params['left'] = settings.gbd_ages[:-1]
        params['width'] = np.diff(settings.gbd_ages)

        age_intervals = zip(settings.gbd_ages[:-1], settings.gbd_ages[1:])
        mean = dm.get_mcmc('mean', k)
        if len(mean) >= settings.gbd_ages[-1]:
            params['height'] = [rate_for_range(mean, range(a0,a1), np.ones(a1-a0)/(a1-a0)) for a0,a1 in age_intervals]
            color = color_for.get(type, 'black')
            params['color'] =  color
            params['edgecolor'] = color
            params['alpha'] = .25
            pl.bar(**params)

            # plot error bars
            params = {}
            params['x'] = np.mean(age_intervals, 1)
            params['y'] =  [rate_for_range(mean, range(a0,a1), np.ones(a1-a0)/(a1-a0)) for a0,a1 in age_intervals]

            err_below = params['y'] - np.array([rate_for_range(dm.get_mcmc('lower_ui', k), range(a0,a1), np.ones(a1-a0)/(a1-a0)) for a0,a1 in age_intervals])
            err_above = np.array([rate_for_range(dm.get_mcmc('upper_ui', k), range(a0,a1), np.ones(a1-a0)/(a1-a0)) for a0,a1 in age_intervals]) - params['y']
            params['yerr'] = [err_below, err_above]

            params['fmt'] = None

            color = color_for.get(type, 'black')
            params['ecolor'] = color
            params['color'] = color

            pl.errorbar(**params)

        # set axis and decorate figure
        xmin, xmax, ymin, ymax = pl.axis()
        xmin = 0
        xmax = 100
        ymin = 0.
        pl.yticks([ymin, ymax], fontsize=8)
        if type == 'relative-risk':
            ymin = 1.
            ymax = 5.
        elif type == 'duration':
            ymax = 100
        else:
            ymax = max([default_max_for[type]] + [dm.value_per_1(d) for d in data_hash.get(type)])

        pl.axis([xmin, xmax, ymin, ymax])

        pl.text(xmin, ymin, ' %s\n\n'%type, fontsize=8)

        pl.xticks([])
    pl.xticks(range(10,100,10), fontsize=8)
        
            
def tile_plot_disease_model(dm_json, keys, defaults={}):
    """Make a graphic representation of the disease model data and
    estimates provided

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list
      the keys to include
    defaults : dict, optional
      additional parameters for the plot.  allowed (key, value) pairs::

        fontsize, integer > 0
        ticksize, integer > 0
        data_alpha, float in (0,1)
        region_labels, boolean
        ymax, float > 0
        label, string
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return
    data_hash = GBDDataHash(dm.data)

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] != 'bins']

    cnt = len(keys)
    if cnt == 1:
        rows = 1
        cols = 1

        subplot_width = 12
        subplot_height = 8
    else:
        cols = 4
        rows = int(np.ceil(float(cnt) / float(cols)))
        
        subplot_width = 6
        subplot_height = 4

    clear_plot(width=subplot_width*cols,height=subplot_height*rows)

    subplot_by_type = {}
    
    for ii, k in enumerate(keys):
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)

        cur_subplot = pl.subplot(rows, cols, ii + 1, sharey=subplot_by_type.get(type))
        subplot_by_type[type] = cur_subplot

        data_type = clean(type) + ' data'

        # special case for single parameter models
        if type == 'continuous_single_parameter':
            data_type = 'continuous single parameter'
        
        data = data_hash.get(data_type, region, year, sex)
        plot_intervals(dm, data, color=color_for.get(data_type, 'black'), print_sample_size=True, alpha=defaults.get('data_alpha', .8))
        data = data_hash.get(data_type, region, year, 'total')
        plot_intervals(dm, data, color='gray', linewidth=3, alpha=defaults.get('data_alpha', .8))

        # if data_type is excess mortality, also include plot of cause-specific morltaity as a lowerbound
        if clean(type) == 'excess-mortality':
            data_type = 'cause-specific mortality data'
            data = data_hash.get(data_type, region, year, sex)
            plot_intervals(dm, data, color=color_for.get(data_type, 'black'), print_sample_size=True, alpha=.8)
                    
        plot_truth(dm, k, color=color_for.get(type, 'black'))
        plot_empirical_prior(dm, k, color=color_for.get(type, 'black'))
        if region == 'all':
            if dm.params.has_key('empirical_prior_%s' % type):
                emp_prior_effects = json.loads(dm.params['empirical_prior_%s' % type])
                alpha = np.array(emp_prior_effects['alpha'])
                beta = np.array(emp_prior_effects['beta'])
                gamma = np.array(emp_prior_effects['gamma'])
                pl.plot(np.exp(np.mean(alpha)+gamma), color=color_for.get(type, 'black'), alpha=.8, linewidth=2, linestyle='dashed')

                for ii, alpha_a in enumerate(alpha[:-2]):
                    color = pl.cm.Spectral(ii/21.)
                    for alpha_t in [alpha[-2], -alpha[-2]]:
                        for alpha_s in [alpha[-1], -alpha[-1]]:
                            x = np.arange(MAX_AGE)
                            y = np.exp(alpha_a + .5*alpha_s + .1*7*alpha_t + gamma)
                            pl.plot(x, y, color=color, alpha=.75, linewidth=2, linestyle='solid')
                            if defaults.get('region_labels'):
                                pl.text(x[-1], y[-1], dismod3.gbd_regions[ii], va='top', ha='right', alpha=.7, color=np.array(color)/2)
                
        #if not dm.has_mcmc(k):
        #    plot_map_fit(dm, k, color=color_for.get(type, 'black'))
        plot_mcmc_fit(dm, k, color=color_for.get(type, 'black'))

        rate_list = [default_max_for.get(type, .0001)] + [dm.value_per_1(d) for d in dm.data if dismod3.relevant_to(d, type, 'all', 'all', 'all')]
        max_rate = np.max(rate_list)
        ages = dm.get_estimate_age_mesh()

        xmin, xmax, ymin, ymax = pl.axis()
        xmin = ages[0]
        xmax = ages[-1]
        ymin = max(ymin, 0.)
        ymax = float(defaults.get('ymax', max(ymax, 1.25*max_rate)))
        pl.axis([xmin, xmax, ymin, ymax])

        plot_prior(dm, k)
        if type == 'mortality':
            type = 'with-condition mortality'
        type = defaults.get('label', type)
        fontsize = defaults.get('fontsize', len(keys) == 1 and 20 or 10)
        ticksize = defaults.get('ticksize', len(keys) == 1 and 20 or 10)
        
        label_plot(dm, type, fontsize=fontsize)
        pl.title('%s %s; %s, %s, %s' % (prettify(dm.params['condition']), type, prettify(region), sex, year), fontsize=fontsize)
        pl.axis([xmin, xmax, ymin, ymax])

        pl.xticks(fontsize=ticksize)
        t, n = pl.yticks()
        pl.yticks([t[0], t[len(t)/2], t[-1]], fontsize=ticksize)

        if len(keys) == 1:
            pl.subplots_adjust(.1, .1, .95, .95)
            
def sparkline_plot_disease_model(dm_json, keys, max_intervals=50, defaults={}):
    """Make a small graphic representation of the disease model data and
    estimates provided

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that
      is to be plotted
    keys : list
      the keys to include
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return
    data_hash = GBDDataHash(dm.data)

    keys = [k for k in keys if k.split(KEY_DELIM_CHAR)[0] != 'bins']

    rows = 1
    cols = 1

    subplot_width = 1
    subplot_height = .6
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    
    for ii, k in enumerate(keys):
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)

        if dm.has_mcmc(k):
            plot_mcmc_fit(dm, k, color=color_for.get(type, 'black'))
        else:
            plot_empirical_prior(dm, k, color=color_for.get(type, 'black'))
        pl.axis('off')


def sparkplot_boxes(dm_json):
    """ Find pixels for all boxes in the sparkplot lattice below."""
    return sparkplot_disease_model(dm_json, boxes_only=True)

def sparkplot_disease_model(dm_json, max_intervals=50, boxes_only=False):
    """ Make a lattice of sparkplots for the disease_model, with rows
    corresponding to regions, and columns corresponding to (year,sex)
    pairs.

    Parameters
    ----------
    dm_json : str or DiseaseJson object
      the json string or a thin python wrapper around this data that is to be plotted
    max_intervals : int, optional
    boxes_only : bool, optional
      if boxes_only == True, then don't plot the sparks, just
      calculate the dimensions that they would have (for an imagemap)
    """
    if isinstance(dm_json, DiseaseJson):
        dm = dm_json
    else:
        try:
            dm = DiseaseJson(dm_json)
        except ValueError:
            print 'ERROR: dm_json is not a DiseaseJson object or json string'
            return
        
        
    # divide up disease_model data by data_type, region, year, sex
    data_hash = GBDDataHash(dm.data)

    # divide up canvas to fit all the sparklines
    rows = len(dismod3.gbd_regions)+1
    col_list = [[1990, 'male'],
                [2005, 'male'],
                [1990, 'female'],
                [2005, 'female'],
                ]
    cols = len(col_list)

    subplot_width = 1. * .5
    subplot_height = .5 * .5
    fig_width = subplot_width*cols
    fig_height = subplot_height*rows

    if not boxes_only:
        fig = pl.figure(figsize=(fig_width,fig_height), dpi=100)
        pl.clf()
    subplot_px = {}

    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    rate_list = [.0001] + [dm.value_per_1(d) for d in dm.data if dismod3.relevant_to(d, 'prevalence', 'all', 'all', 'all')]
    ymax = np.max(rate_list)
    
    sorted_regions = sorted(dismod3.gbd_regions, reverse=False,
                            key=lambda r: len(data_hash.get(region=r)))
    for ii, region in enumerate(sorted_regions + ['all']):
        for jj, [year, sex] in enumerate(col_list):
            subplot_px[dismod3.gbd_key_for('all', region, year, sex)] = \
                ', '.join([str(int(100 * (jj) * subplot_width)),
                           str(int(100 * (rows - ii - 1) * subplot_height)),
                           str(int(100 * (jj + 1) * subplot_width)),
                           str(int(100 * (rows - ii) * subplot_height))])
            if boxes_only:
                continue

            fig.add_axes([jj*subplot_width / fig_width,
                          ii*subplot_height / fig_height,
                          subplot_width / fig_width,
                          subplot_height / fig_height],
                         frameon=False)
            # plot data and fit for each data type in a different color
            for type in ['prevalence', 'incidence', 'all-cause mortality']:
                plot_fit(dm, 'mcmc_mean', dismod3.gbd_key_for(type, region, year, sex),
                             linestyle='-', color=color_for.get(type, 'black'), linewidth=1, alpha=.8)
                type = ' '.join([type, 'data'])
                data = data_hash.get(type, region, year, sex) + data_hash.get(type, region, year, 'total')
                #if len(data) > max_intervals:
                #    data = random.sample(data, max_intervals)
                plot_intervals(dm, data, color=color_for.get(type, 'black'), linewidth=1, alpha=.25)
            pl.xticks([])
            pl.yticks([])
            pl.axis([xmin, xmax, ymin, ymax])

    if boxes_only:
        return subplot_px

def plot_prior_preview(dm):
    """ Generate a preview of what a rate function with this prior looks like"""

    fig = pl.figure(figsize=(3.4, 3.4), dpi=100)
    pl.clf()

    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    ymax = dm.get_ymax()

    for ii, type in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality']):
        pl.subplot(2,2,ii+1)
        prior_str = dm.get_global_priors(type)

        vars = dismod3.utils.prior_vals(dm, type)
        key = 'prevalence+asia_southeast+1990+male'  # just pick something, what doesn't matter
        dismod3.utils.generate_prior_potentials(vars, dm.get_priors(key), dm.get_estimate_age_mesh())
        
        dm.vars = vars
        ages = dm.get_estimate_age_mesh()
        color = color_for.get(type, 'black')
        #mu = vars['rate_stoch'].stats()['mean']
        prior_vals = dict(
            alpha=list(dm.vars['region_coeffs']),
            beta=list(dm.vars['study_coeffs']),
            gamma=list(dm.vars['age_coeffs'].value),
            delta=float(dm.vars['dispersion'].value))
        from neg_binom_model import predict_rate, regional_covariates
        mu = predict_rate(regional_covariates(key, dm.get_covariates()),
                          alpha=prior_vals['alpha'],
                          beta=prior_vals['beta'],
                          gamma=prior_vals['gamma'],
                          bounds_func=dm.vars['bounds_func'],
                          ages=dm.get_param_age_mesh())
        dispersion = vars['dispersion'].stats()['mean']

        pl.plot(ages, mu, color=color, linestyle='-', linewidth=2)
        lb = mu*(1 - 1/np.sqrt(dispersion))
        ub = mu*(1 + 1/np.sqrt(dispersion))
        plot_uncertainty(ages, lb, ub, edgecolor=color, alpha=.75)

        plot_prior(dm, type)
        plot_intervals(dm, vars['data'], color=color)

        x0, x1, ymin, ymax = pl.axis()
        ymax = max(ymax, .0001)
        pl.text(xmin, ymax, type, color=color,
                verticalalignment='top', horizontalalignment='left')
        
        pl.axis([xmin, xmax, 0., ymax])
        pl.xticks([])
        pl.yticks(fontsize=8)

def simplify_ticks():
    ticks, labels = pl.xticks()
    pl.xticks(ticks)
    ticks, labels = pl.yticks()
    pl.yticks([ticks[0], 0, ticks[-1]])

def plot_empirical_prior_effects(dm_list, effect, **params):
    types = ['prevalence', 'incidence', 'remission', 'excess-mortality']
    text_size = 8
    label_size=10
    msg_params = dict(fontsize=text_size, alpha=.6, color=(.7,.2,.2))
    
    if effect in ['alpha', 'beta']:
        # setup canvas
        if effect == 'alpha':
            covariate_list = dismod3.gbd_regions + ['Year', 'Sex']
        else:
            from neg_binom_model import covariate_names
            covariate_list=covariate_names(dm_list[0]) # TODO: make this work for comparing models with different covariates

        clear_plot(width=8, height=3. + len(covariate_list)/6.)

        if effect == 'alpha':
            pl.figtext(0, 1, 'Country/Sex/Time Effects', fontsize=15, va='top')
        else:
            pl.figtext(0, 1, 'Additional Covariate Effects', fontsize=15, va='top')

        for jj, t in enumerate(types):
            k = 'empirical_prior_%s' % t
            color = color_for.get(t, 'black')
            msg_params['color'] = (.7,.2,.2)
            pl.text(jj+1, len(covariate_list) - .5, t + '\n', va='center', ha='center',
                    size=label_size, color=color)
            pl.vlines([jj+1], -1, len(covariate_list)+1, color=color, linewidth=2, linestyle='dashed')

            dm_len = float(len(dm_list))
            for ii, dm in enumerate(dm_list):
                # TODO: mix color with something to indicate different models
                #color = mix(color, pl.cm.spectral(ii/dm_len))
                
                emp_p = json.loads(dm.params.get(k, '{}'))
                val = emp_p.get(effect)
                se = emp_p.get('sigma_%s'%effect)

                if not val or not se:
                    pl.text(1+jj+ii/5., 0., (' no empirical prior in model %d' % dm.id) + '\n'*ii, size=label_size, color=color, alpha=.5, rotation=90)
                else:
                    val = np.atleast_1d(val)
                    se = np.atleast_1d(se)

                    err = [np.exp(val) - np.exp(val - 1.96*se), np.exp(val + 1.96*se) - np.exp(val)]
                    pl.errorbar(jj + np.exp(val),
                                np.arange(len(val)) + len(covariate_list)*jj/200. + ii/5.,
                                xerr=err,
                                color=color,
                                fmt='o', alpha=.8)

        # decorate figure
        pl.hlines(range(23), -5, 22.5, linestyle='dotted', alpha=.5, zorder=-1)
        pl.xticks([1,2,3,4], ['1', '1', '1', '1'])
        pl.yticks(range(len(covariate_list)), covariate_list, size=label_size)
        pl.axis([0, jj+2, -.5, len(covariate_list)-.5])
        pl.subplots_adjust(.3, .08, .95, .92)

    elif effect == 'gamma':
        # setup canvas
        clear_plot(height=8, width=8)
        pl.figtext(0, 1, 'Empirical Prior Age Patterns', fontsize=15, va='top')
        
        for jj, t in enumerate(types):
            k = 'empirical_prior_%s' % t
            color = color_for.get(t, 'black')
            msg_params['color'] = (.7,.2,.2)
            pl.text(MAX_AGE-1, jj+1, '\n'+t, va='top', ha='right',
                    size=label_size, color=color)
            pl.hlines([jj+1], 0, MAX_AGE, color=color, linewidth=2, linestyle='dashed', zorder=-1)

            dm_len = float(len(dm_list))
            for ii, dm in enumerate(dm_list):
                # TODO: mix color with something to indicate different models
                #color = mix(color, pl.cm.spectral(ii/dm_len))
                
                emp_p = json.loads(dm.params.get(k, '{}'))
                val = emp_p.get(effect)
                se = emp_p.get('sigma_%s'%effect)

                if not val or not se:
                    pl.text(0., 1+jj, (' no empirical prior in model %d' % dm.id) + '\n'*ii, size=label_size, color=color, alpha=.5)
                else:
                    val = np.atleast_1d(val)
                    se = np.atleast_1d(se)

                    mu = np.mean(np.exp(val))
                    sigma = np.std(np.exp(val)) * 4

                    lb = (np.exp(val - 1.96*se) - mu) / sigma
                    ub = (np.exp(val + 1.96*se) - mu) / sigma
                    val = (np.exp(val) - mu) / sigma

                    pl.plot(np.arange(len(val)),
                            1 + jj + val,
                            color=color,
                            alpha=.8)
                    plot_uncertainty(np.arange(len(val)),
                                     1 + jj + lb,
                                     1 + jj + ub, edgecolor='none')
        # decorate figure
        pl.yticks([])
        pl.axis([0, MAX_AGE-1, 0, jj+2])
        #pl.subplots_adjust(.3, .08, .95, .92)

    elif effect == 'delta':
        fig_width = 12
        fig_height = 4
        fig_rows = 4
        fig_cols = 1

        fig = clear_plot(width=fig_width, height=fig_height)
        ax = pl.subplot(fig_rows, fig_cols, 1)

        pl.figtext(0, 1, 'Over-dispersion Prior and Goodness of Fit Statistics', fontsize=15, va='top')
            
        for i, t in enumerate(['prevalence', 'incidence', 'remission', 'excess-mortality']):
            pl.xticks(fontsize=text_size)
            pl.yticks(fontsize=text_size)
            k = 'empirical_prior_%s' % t

            dm_len = float(len(dm_list))
            for ii, dm in enumerate(dm_list):
                color = pl.cm.spectral(ii/dm_len)

                if dm.params.has_key(k):
                    emp_p = json.loads(dm.params[k])
                    val = emp_p.get(effect)
                    se = emp_p.get('sigma_%s'%effect)
                else:
                    val = None

                if not val or not se:
                    pl.text(0., 0., (' no empirical prior for %s in model %d' % (t, dm.id)) + '\n'*ii, **msg_params)
                    simplify_ticks()
                    if effect == 'delta':
                        pl.axis('off')
                    continue

                val = np.atleast_1d(val)
                se = np.atleast_1d(se)

                if effect == 'alpha':
                    err = [np.exp(val) - np.exp(val - 1.96*se), np.exp(val + 1.96*se) - np.exp(val)]
                    pl.errorbar(np.exp(val), np.arange(len(val)), xerr=err, fmt='o', )

                elif effect == 'beta':
                    pl.errorbar(np.arange(len(val))+.5/dm_len+ii/dm_len, val, 1.96*se, fmt=None, color=color)
                    pl.bar(np.arange(len(val))+ii/dm_len, val, 1/dm_len-.1, alpha=.5, color=color)
                    pl.xticks([], [])
                    simplify_ticks()

                elif effect == 'gamma':
                    pl.errorbar(range(len(val)), val, 1.96*se,
                                fmt=None, alpha=.15, color=color)
                    pl.plot(range(len(val)), val, color=color)
                    l,r,b,t = pl.axis()
                    pl.axis([0, 100, b, t])
                    simplify_ticks()
                elif effect == 'delta':
                    pl.axis('off')
                    info_str = '$\delta_%s^{%s} = %3.2f \pm %.2f$;' % (t[0], dm.id, val, 1.96*se)
                    for gof in ['aic', 'bic', 'dic']:
                        if gof in emp_p:
                            info_str += '\t $%s_%s^{%s} = %.0f$;' % (gof.upper(), t[0], dm.id, emp_p[gof])
                    pl.text(0., 0., info_str + '\n'*ii, fontsize=16, va='bottom', ha='left', color=color)

            if i == max(fig_rows, fig_cols)-1:
                continue
            if fig_rows == 1:
                ax = pl.subplot(fig_rows, fig_cols, i+2, sharex=ax)
                pl.xlabel(t, fontsize=text_size, color=color_for.get(t, 'black'))
                msg_params['rotation'] = 90
            else:
                ax = pl.subplot(fig_rows, fig_cols, i+2, sharey=ax)
                pl.ylabel(t, fontsize=text_size, color=color_for.get(t, 'black'))
                msg_params['rotation'] = 0

    else:
        raise AttributeError('Unknown effect type %s'%effect)

                
def plot_intervals(dm, data, print_sample_size=False, **params):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of Data
    """
    default_params = dict(alpha=.35, color=(.0,.5,.0), linewidth=2)
    default_params.update(**params)
    
    errorbar_params = copy.copy(default_params)
    errorbar_params['linewidth'] = 1

    # include at most 200 data bars
    #if len(data) > 200:
    #    import random
    #    data = random.sample(data, 200)

    for d in data:
        if d['age_end'] == MISSING:
            d['age_end'] = MAX_AGE

        val = dm.value_per_1(d)
        if val == MISSING:
            continue
        
        lb, ub = dm.bounds_per_1(d)
        if lb != ub:  # don't draw error bars if interval is zero or MISSING
            pl.plot([.5 * (d['age_start']+d['age_end'])]*2,
                    [lb, ub],
                    **errorbar_params)
        
        pl.plot(np.array([d['age_start']-.5, d['age_end']+.5]),
                np.array([val, val]),
                **default_params)

        #if clean(d.get('self_reported', '')) == 'true':
        #    pl.text(.5*(d['age_start']+d['age_end']), val, 'self-reported', fontsize=6, horizontalalignment='center', verticalalignment='center')
        if print_sample_size:
            pl.text(.5*(d['age_start']+d['age_end']), val, d.get('effective_sample_size', ''), fontsize=6, horizontalalignment='center', verticalalignment='center')

def plot_fit(dm, fit_name, key, **params):
    fit = dm.params.get(fit_name, {}).get(key)
    age = dm.get_estimate_age_mesh()
    if fit and age:
        pl.plot(age, fit, **params)

def plot_initial_estimate(dm, type, **params):
    default_params = {'color': 'blue', 'alpha': .5}
    default_params.update(**params)
    plot_fit(dm, 'initial_estimate', type, **default_params)

def plot_truth(dm, type, **params):
    default_params = {'color': 'blue',
                      'linestyle': '--',
                      'linewidth': 1,
                      'alpha': .9,
                      'label': 'Ground Truth',
                      }
    default_params.update(**params)
    default_params.update(color='black')
    plot_fit(dm, 'truth', type, zorder=11, **default_params)
    plot_fit(dm, 'truth', type, color='yellow', alpha=.5, linewidth=5, zorder=10)

def plot_map_fit(dm, type, **params):
    default_params = {'color': 'blue',
                      'linestyle': 'solid',
                      'linewidth': .5,
                      'alpha': 1.,
                      'label': 'Max-liklihood',
                      }
    default_params.update(**params)
    plot_fit(dm, 'map', type, **default_params)

def plot_mcmc_fit(dm, type, color=(.2,.2,.2), show_data_ui=True):
    age = dm.get_estimate_age_mesh()
    param_mesh = dm.get_param_age_mesh()

    # adjust ages to make proper stair-step curve
    for ii, a in enumerate(age):
        if a in param_mesh and ii > 0:
            age[ii-1] += .5
            age[ii] -= .5
    
    lb = dm.get_mcmc('lower_ui', type)
    ub = dm.get_mcmc('upper_ui', type)

    if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
        plot_uncertainty(age, lb, ub, facecolor='none', edgecolor=color, linewidth=.5, alpha=1, zorder=2.)

    val = dm.get_mcmc('median', type)

    if len(age) > 0 and len(age) == len(val):
        pl.plot(age, val, color=color, linewidth=1, alpha=1., zorder=100)
    #left = param_mesh[:-1]
    #height = val[left]
    #width = pl.diff(param_mesh)
    #yerr = pl.array([val[left] - lb[left], ub[left] - val[left]])
    #pl.errorbar(left+.5*width, height, yerr=yerr, linestyle='none',
    #            elinewidth=3, color=color, capsize=pl.mean(width), zorder=101)
    #pl.bar(left, height, width,
    #       facecolor='none', linewidth=2,
    #       edgecolor=color, alpha=1.)
def plot_empirical_prior(dm, type, **params):
    default_params = dict(color=(.2,.2,.2), linewidth=3, alpha=1., linestyle='dashed')
    default_params.update(**params)
    
    age = dm.get_estimate_age_mesh()
    #lb = dm.get_mcmc('emp_prior_lower_ui', type)
    #ub = dm.get_mcmc('emp_prior_upper_ui', type)
    #if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
    #    plot_uncertainty(age, lb, ub, linestyle='dotted', fill=True, edgecolor=color, alpha=.5, zorder=0.)

    val = dm.get_mcmc('emp_prior_mean', type)
    if len(age) > 0 and len(age) == len(val):
        pl.plot(age, val, **default_params)

def plot_uncertainty(ages, lower_bound, upper_bound, **params):
    default_params = {'facecolor': '.8'}
    default_params.update(**params)

    x = np.concatenate((ages, ages[::-1]))
    y = np.concatenate((lower_bound, upper_bound[::-1]))
    pl.fill(x, y, **default_params)

def plot_mcmc_diagnostics(rate_stoch):
    e = rate_stoch.trace() - rate_stoch.stats()['mean']
    cols = 3
    for i in range(10):
        pl.subplot(10, cols, cols*i+1)
        pl.acorr(e[:, i*10], normed=True)
        pl.xticks([])
        pl.yticks([0, 1])
        
        pl.subplot(10, cols, cols*i+2)
        pl.plot(e[:, i*10])
        pl.xticks([])
        pl.yticks([])
        
        pl.subplot(10, cols, cols*i+3)
        pl.hist(e[:, i*10])
        pl.xticks([])
        pl.yticks([])

def plot_posterior_predicted_checks(dm, key):
    vars = dm.vars[key]

    if not vars.get('observed_counts'):
        pl.figure()
        pl.text(.5, .5, 'no data')
        return
    
    n = len(vars['observed_counts'].value)
    k = len(vars['predicted_rates'].trace())


    pl.figure(figsize=(max(6, .75*n), 8))

    observed_rates = pl.array(vars['observed_counts'].value)/vars['effective_sample_size']
    observed_std = pl.sqrt(observed_rates * (1 - observed_rates) / vars['effective_sample_size'])

    hpd = vars['predicted_rates'].stats()['95% HPD interval']
    sorted_indices = pl.argsort(
        pl.where((observed_rates < hpd[:,0]) | (observed_rates > hpd[:,1]), 1000., 1.)
        + observed_rates
        )   # sort by observed_rate in predicted_rate 95% HPD, and then by value

    pl.plot([-1], [-1], 'go', mew=0, ms=10, label='Data Predicted Rate')
    pl.plot((pl.outer(pl.ones(k), range(n)) + pl.randn(k, n)*.1).flatten(),
            vars['predicted_rates'].trace()[:, sorted_indices].flatten(),
            'g.', alpha=.5)

    pl.plot([-1], [-1], 'bo', mew=0, ms=10, label='Expected Rate')
    pl.plot((pl.outer(pl.ones(k), range(n)) + pl.randn(k, n)*.05).flatten(),
            vars['expected_rates'].trace()[:, sorted_indices].flatten(),
            'b.', alpha=.5)

    pl.plot([-1], [-1], 'ro', mew=0, ms=10, label='Observed Rate')
    pl.errorbar(range(n), observed_rates[sorted_indices], yerr=observed_std[sorted_indices]*1.96,
                fmt='ro', mew=0, ms=10, alpha=.75,
                barsabove=True, zorder=10)

    l,r,b,t = pl.axis()
    t = max(observed_rates)*1.5

    # visually identify panels where observed_rate falls outside of predicted_rate 95% HPD interval
    coverage = sum((observed_rates >= hpd[:,0]) & (observed_rates <= hpd[:,1]))
    if coverage == len(observed_rates):
        cov_str = '\n No data outside of\n95% HPD interval'
    else:
        cov_str = '\n Data outside of\n95% HPD interval'
    #pl.vlines([coverage - .5], 0, t, color='black')
    pl.text(coverage-.5, t, cov_str, fontsize=8, va='top', ha='left')

    # set axis limits for main panel
    pl.axis([-.5, n-.5, 0, t])
    pl.xticks([])
    pl.legend(loc='upper left', numpoints=1)
    pl.title('%s\nData Posterior Predictions' % key)
    pl.subplots_adjust(left=.1, right=.9, bottom=.3)

    dx = (.9 - .1) / n
    Xa = vars['expected_rates'].parents['Xa']
    Xb = vars['expected_rates'].parents['Xb']
    N = vars['effective_sample_size']
    
    for ii, jj in enumerate(sorted_indices):
        pl.axes([.1 + ii*dx, .1, dx, .2])
        x = vars['expected_rates'].trace()[:, jj]
        pl.acorr(x, normed=True, detrend=pl.mlab.detrend_mean, usevlines=True, maxlags=20,)
        pl.xticks([])
        pl.yticks([])
        pl.axis([-20,20,-.2, 1])
        pl.text(-20, 1.,
                '$X_t=%.1f$\n$X_s=%.1f$\n$X_r=%d$\n$X_b=%s$\n$N=%.1f$' % (Xa[jj][21], Xa[jj][22], pl.argmax(Xa[jj][:21]), Xb[jj], N[jj]),
                fontsize=10, va='top', ha='left', color='red')
        # TODO: refactor following code into a function, and use it to decorate age-pattern plots upon request
        d = vars['data'][jj]
        label_str = ''
        label_str += 'Row %d' % d['_row']
        label_str += '\n' + d['country_iso3_code']
        label_str += '\n' + d['sex']
        label_str += '\n%s-%s' % (d['age_start'], d['age_end'])
        label_str += '\n%s-%s' % (str(d['year_start'])[-2:], str(d['year_end'])[-2:])
        label_str += '\n%s' % (d['parameter'][:4])
        pl.xlabel(label_str, fontsize=8)
        if ii == 0:
            pl.ylabel('mcmc acorr', fontsize=8)

    info_str = '%d: ' % dm.id
    for s, v in [['DIC', dm.mcmc.dic()]]:
        info_str += '$%s = %.0f$;\t' % (s, v)
    pl.figtext(0., 1., info_str, fontsize=12, va='top', ha='left')
        

def plot_prior(dm, type):
    # write out details of priors in a friendly font as well
    if len(dm.get_estimate_age_mesh()) > 0:
        l, r, b, t = pl.axis()
        a0 = dm.get_estimate_age_mesh()[0]
        v0 = b
        pl.text(a0, v0, ' Priors:\n' + dm.get_priors(type).replace(dismod3.PRIOR_SEP_STR, '\n'), color='black', family='monospace', fontsize=8, alpha=.75)

    # show level value priors
    for prior_str in dm.get_priors(type).split(dismod3.PRIOR_SEP_STR):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'level_value':
            level_val = float(prior[1])
            age_start = int(prior[2])
            age_end = int(prior[3])

            pl.plot([age_start, age_end+1], [level_val, level_val], color='red', linewidth=15, alpha=.75)
        
def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    pl.subplots_adjust(.025, .025, .99, .965)
    return fig

def label_plot(dm, type, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('%s %s' % (type, dm.get_units(type) or ''), **params)
    pl.title('%s: %s; %s; %s; %s' % \
                 (dm.params['id'], prettify(dm.params['condition']),
                  dm.params['sex'], prettify(dm.params['region']),
                  dm.params['year']), **params)
    #pl.legend()

def choropleth_dict(title, region_value_dict, scheme, data_type='int'):
    """ make a map plot for integer values

    Parameters 
    ----------
    title : string
      Title of the map
    region_value_dict : dictionary
      GBD region name versus integer value of the region
    scheme : string
    data_type : string
      'int' or 'float', which are treated differently

    Returns
    region_color_dict : dictionary
      GBD region name versus color code of the region
    region_value_dict : dictionary
      GBD region name versus value of the region
    legend : list
      list of legend colors
    bin_name_list : list
      list of bin names
    title : string
       Title of the map
    note : string
       Note of the map
    """
    value_list = [v for v in region_value_dict.values() if not np.isnan(v)]
    if len(value_list) == 0:
        return None

    legend = ['00ffff', '00ff00', 'aad400', 'ffcc00', 'ff7f2a', 'ff0000', 'ffffff']
    if data_type == 'int':
        max_v = max(value_list)
        if max_v < 12:
            max_v = 12
        bin_size = int(math.ceil(float(max_v) / 6))

        region_color_dict = {}
        for key in region_value_dict:
            region_color_dict[key] = legend[int(math.ceil(float(region_value_dict[key]) / bin_size)) - 1]

        bin_name_list = []
        for i in range(6):
            bin_name_list.append('%d - %d' % ((i * bin_size) + 1, (i + 1) * bin_size))
    elif data_type == 'float':
        if scheme == 'uniform':
            bin_size = 0.00001
            max_v = max(value_list)
            min_v = min(value_list)
            if max_v != 0:
                bin_size = float(max_v - min_v) / 6
            region_color_dict = {}
            for key in region_value_dict:
                if np.isnan(region_value_dict[key]):
                    region_color_dict[key] = legend[6]
                else:
                    if bin_size != 0:
                        color_index = int(math.floor(float(region_value_dict[key] - min_v) / bin_size))
                    else:
                        color_index = 0
                    if color_index == 6:
                        color_index = 5
                    region_color_dict[key] = legend[color_index]
            bin_name_list = []
            for i in range(6):
                bin_name_list.append('%g - %g' % (format(min_v + bin_size * i), format(min_v + bin_size * (i + 1))))

        elif scheme == 'colorful':
            # sort the region values
            items = region_value_dict.items()
            items.sort(key = itemgetter(1))

            region_color_dict = {}
            n = len(region_value_dict)
            for i in items:
                if np.isnan(i[1]):
                    n -= 1
                    region_color_dict[i[0]] = legend[6]

            # calculate numbers of regions in each bin
            m = n % 6
            l = n / 6
            n_list = [0, 0, 0, 0, 0, 0]
            for i in range(6):
                if i < m:
                    n_list[i] = l + 1
                else:
                    n_list[i] = l

            # set region color and bin color name
            i = 0
            bin_name_list = ['', '', '', '', '', '']
            for j in range(6):
                for k in range(n_list[j]):
                    if n_list[j] != 0:
                        region_color_dict[items[i][0]] = legend[j]
                        i += 1
                if n_list[j] == 1:
                    bin_name_list[j] = '%g' % format(items[i - n_list[j]][1])
                else:
                    bin_name_list[j] = '%g - %g' % (format(items[ i - n_list[j]][1]), format(items[i - 1][1]))
            for i in range(6):
                if n_list[i] == 0:
                    legend[i] = legend[6]
                    bin_name_list[i] = ''

    # remove dashes from key names, since django templates can't handle them
    for r in region_color_dict.keys():
        rr = r.replace('-', '_')
        region_color_dict[rr] = region_color_dict[r]
        region_value_dict[rr] = region_value_dict[r]
        if rr != r:
            del region_value_dict[r]

    # format float numbers
    if data_type == 'float':
        for r in region_value_dict.keys():
            if str(region_value_dict[r]) != 'nan':
                region_value_dict[r] = '%g' % format(region_value_dict[r])
            else:
                region_value_dict[r] = ''

    # make a note
    note = ''
    for r in region_value_dict.keys():
        if (data_type == 'int' and region_value_dict[r] == 0) or (data_type == 'float' and region_value_dict[r] == ''):
            data = title.split(':')[0]
            if data == 'Data Count':
                data = 'data'
            note = 'Regions in white color have no ' + data.lower()
            break

    #return dict(color=region_color_dict, value=region_value_dict, label=bin_name_list, title=title)
    return dict(color=region_color_dict, value=region_value_dict, legend=legend, label=bin_name_list, title=title, note=note)

def format(v):
    s = float(v)
    if s != 0:
        p = math.pow(10, math.floor(math.log10(s)) - 3)
        return math.ceil(s / p) * p
    else:
        return s

def plot_posterior_selected_regions(region_value_dict, condition, type, year, sex, ages, xmin, xmax, ymin, ymax, grid, linewidth):
    """Make a graphic representation of the posterior of disease model for selected regions

    Parameters
    ----------
    region_value_dict : dictionary
      year or sex versus rate value
    condition : str
      GBD cause
    type : str
      one of the eight parameter types
    year : str
      1990 or 2005
    sex :  str
      male or female
    ages : list
      estimated age mesh
    xmin : int
      X axis lower bound
    xmax : int
      X axis upper bound
    ymin : str
      Y axis lower bound
    ymax : str
      Y axis upper bound
    grid : str
      True or False
    linewidth : float
      line width [.1, 10]
    """
    params = {
        'Asia Pacific, High Income' : dict(color=(1,0,0), linestyle='-', marker='x'),
        'Latin America, Southern' : dict(color=(1,0,0), linestyle='-', marker='^'),
        'Australasia' : dict(color=(1,0,0), linestyle='-', marker='s'),
        'North America, High Income' : dict(color=(1,0,0), linestyle='-', marker='o'),
        'Asia, Central' : dict(color='b', linestyle='-', marker='x'),
        'Asia, East' : dict(color='b', linestyle='-', marker='o'),
        'Asia, South' : dict(color='b', linestyle='-', marker='s'),
        'Asia, Southeast' : dict(color='b', linestyle='-', marker='^'),
        'North Africa/Middle East' : dict(color='y', linestyle='-', marker='o'),
        'Oceania' : dict(color='y', linestyle='-', marker='s'),
        'Caribbean' : dict(color='y', linestyle='-', marker='^'),
        'Europe, Central' : dict(color='g', linestyle='-', marker='^'),
        'Europe, Eastern' : dict(color='g', linestyle='-', marker='s'),
        'Europe, Western' : dict(color='g', linestyle='-', marker='o'),
        'Latin America, Andean' : dict(color='c', linestyle='-', marker='^'),
        'Latin America, Central' : dict(color='c', linestyle='-', marker='s'),
        'Latin America, Tropical' : dict(color='c', linestyle='-', marker='o'),
        'Sub-Saharan Africa, Central' : dict(color='m', linestyle='-', marker='x'),
        'Sub-Saharan Africa, East' : dict(color='m', linestyle='-', marker='o'),
        'Sub-Saharan Africa, Southern' : dict(color='m', linestyle='-', marker='s'),
        'Sub-Saharan Africa, West' : dict(color='m', linestyle='-', marker='^'),
        }

    
    pl.figure(figsize=(11.5, 8))
    ax1 = pl.axes([0.1, 0.1, 0.6, 0.8])
    p = []
    style = ''
    t = type
    if t == 'with-condition-mortality':
        t = 'mortality'
    regions = []
    for i, region in enumerate(sorted(region_value_dict.keys())):
        rate = region_value_dict[region]
        if len(rate) == dismod3.MAX_AGE:
            if i > 6:
                style = 'o-'
            if i > 13:
                style = '^-'
            p.append(ax1.plot(ages, rate, alpha=.7, ms=4, linewidth=linewidth, **params[region]))
            regions.append(region)

    pl.xlabel('Age')
    pl.ylabel(type)
    pl.title('Posterior %s %s %s %s' % (prettify(condition), type, sex, year))
    pl.grid(grid)

    minX, xmaX, minY, maxY = pl.axis()
    if not ymin == 'auto':
        minY = float(ymin)
    if not ymax == 'auto':
        maxY = float(ymax)
    pl.axis([xmin, xmax, minY, maxY])

    ax2 = pl.axes([0.72, 0.1, 0.24, 0.8], frameon=False)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)

    l = ax2.legend(p, regions, mode="expand", ncol=1, borderaxespad=0.) 

    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    llines = leg.get_lines()
    frame  = leg.get_frame()

    frame.set_facecolor('0.90')
    pl.setp(ltext, fontsize='small')
    pl.setp(llines, linewidth=linewidth)

def plot_posterior_region(key_value_dict, condition, type, region, key, ages, xmin, xmax, ymin, ymax, grid, linewidth):
    """Make a graphic representation of the posterior of disease model for a region to compare year or sex effects

    Parameters
    ----------
    key_value_dict : dictionary
      year or sex versus rate value
    condition : str
      GBD cause
    type : str
      one of the eight parameter types
    region : str
      one of the GBD regions
    key :  str
      if year: 1990 or 2005, if sex: male, female or total
    ages : list
      estimated age mesh
    xmin : int
      X axis lower bound
    xmax : int
      X axis upper bound
    ymin : str
      Y axis lower bound
    ymax : str
      Y axis upper bound
    grid : str
      True or False
    linewidth : float
      line width [.1, 10]
    """
    pl.figure(figsize=(11.5, 8))
    ax1 = pl.axes([0.1, 0.1, 0.6, 0.8])
    p = []
    style = ''
    t = type
    if t == 'with-condition-mortality':
        t = 'mortality'
    keys = []
    for i, k in enumerate(key_value_dict.keys()):
        rate = key_value_dict[k]
        if len(rate) == dismod3.MAX_AGE:
            if i > 6:
                style = '--'
            if i > 13:
                style = '-.'
            p.append(ax1.plot(ages, rate, style, linewidth=linewidth))
            keys.append(k)

    pl.xlabel('Age')
    pl.ylabel(type)
    pl.title('Posterior %s %s %s %s' % (prettify(condition), type, region, key))
    pl.grid(grid)

    minX, maxX, minY, maxY = pl.axis()
    if not ymin == 'auto':
        minY = float(ymin)
    if not ymax == 'auto':
        maxY = float(ymax)
    pl.axis([xmin, xmax, minY, maxY])

    ax2 = pl.axes([0.72, 0.1, 0.11, 0.8], frameon=False)
    ax2.xaxis.set_visible(False)
    ax2.yaxis.set_visible(False)

    l = ax2.legend(p, keys, mode="expand", ncol=1, borderaxespad=0.) 

    leg = pl.gca().get_legend()
    ltext  = leg.get_texts()
    llines = leg.get_lines()
    frame  = leg.get_frame()

    frame.set_facecolor('0.90')
    pl.setp(ltext, fontsize='small')
    pl.setp(llines, linewidth=linewidth)

class GBDDataHash:
    """ Store and serve data grouped by type, region, year, and sex
    """
    def __init__(self, data):
        self.data = data
        self.d_hash = {}

    def get(self, type='all', region='all', year='all', sex='all'):
        """ Provide a way to get desired data
        
        Parameters
        ----------
        type : str, one of the following types
          'incidence data', 'prevalence data', 'remission data',
          'excess-mortality data', 'all-cause mortality data', 'duration data', 'cause-specific mortality data', or 'all'
        region : str, one of the 21 gbd regions or 'World' or 'all'
        year : int, one of 1990, 2005, 'all'
        sex : str, one of 'male', 'female', 'total', 'all'
        """
        if not self.d_hash.has_key((type, region, year, sex)):
            self.d_hash[(type, region, year, sex)] = [d for d in self.data if dismod3.relevant_to(d, type, region, year, sex)]
        return self.d_hash[(type, region, year, sex)]
