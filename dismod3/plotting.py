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

import copy
import random
import pylab as pl
import numpy as np

from xlwt import *

from time import strftime

import dismod3
from dismod3.utils import clean
from disease_json import *

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
    'relative-risk data': '#ff00ff',
    'relative-risk': '#990099',
    'yld': 'black',
    }

def prettify(str):
    """ Turn underscores into spaces"""
    return str.replace('_', ' ')

def overlay_plot_disease_model(dm_json, keys, max_intervals=100):
    """ Make a graphic representation of the disease model estimates

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

    keys = [k for k in keys if k.split('+')[0] in ['prevalence', 'incidence', 'remission', 'case-fatality']]

    clear_plot(width=6, height=4)
    for k in sorted(keys, key=lambda k: np.max(list(dm.get_map(k)) + [0]), reverse=True):
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)
        data_type = type + ' data'

        # plot the data rectangles for these keys
        data = data_hash.get(data_type, region, year, sex) \
            + data_hash.get(data_type, region, year, 'total')
        if len(data) > max_intervals:
            data = random.sample(data, max_intervals)
        plot_intervals(dm, data, color=color_for.get(data_type, 'black'))

        # plot the map fit
        plot_map_fit(dm, k, linestyle='-',
                     color=color_for.get(type, 'black'),
                     label=k.split('+')[0])

        plot_mcmc_fit(dm, k,
                      color=color_for.get(type, 'black'), show_data_ui=False)

        # plot the empirical prior, if there is any
        plot_empirical_prior(dm, k,
                             color=color_for.get(type, 'black'))

    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    ymax = dm.get_ymax()
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

    label_plot(dm, k, fontsize=10)
    pl.ylabel('')
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
            
def tile_plot_disease_model(dm_json, keys, max_intervals=50):
    """Make a graphic representation of the disease model data and
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

    cnt = len(keys)
    cols = 4 # int(np.sqrt(cnt) + .99)
    rows = int(np.ceil(float(cnt) / float(cols)))

    subplot_width = 6
    subplot_height = 4
    
    clear_plot(width=subplot_width*cols,height=subplot_height*rows)
    
    for ii, k in enumerate(keys):
        pl.subplot(rows, cols, ii + 1)

        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)

        data_type = clean(type) + ' data'
        data = data_hash.get(data_type, region, year, sex) \
               + data_hash.get(data_type, region, year, 'total')
        if len(data) > max_intervals:
            data = random.sample(data, max_intervals)
        plot_intervals(dm, data, color=color_for.get(data_type, 'black'), alpha=.2)
        
        plot_truth(dm, k, color=color_for.get(type, 'black'))
        plot_map_fit(dm, k, color=color_for.get(type, 'black'))
        plot_mcmc_fit(dm, k, color=color_for.get(type, 'black'))
        plot_empirical_prior(dm, k, color=color_for.get(type, 'black'))
        plot_prior(dm, k)
        label_plot(dm, type, fontsize=10)
        pl.title('%s %s; %s, %s, %s' % (prettify(dm.params['condition']), type, prettify(region), sex, year), fontsize=10)

        max_rate = np.max([.001] + [dm.value_per_1(d) for d in dm.data if dismod3.relevant_to(d, type, region, year, sex)]
                          + list(dm.get_map(k))+ list(dm.get_mcmc('mean', k)) + list(dm.get_mcmc('emp_prior_mean', k)))
        ages = dm.get_estimate_age_mesh()
        xmin = ages[0]
        xmax = ages[-1]
        ymin = 0.
        ymax = 1.25*max_rate
        pl.axis([xmin, xmax, ymin, ymax])

        if type == 'yld':
            est_yld = sum(dm.get_mcmc('median', k))
            est_yld_lower_ui = sum(dm.get_mcmc('lower_ui', k))
            est_yld_upper_ui = sum(dm.get_mcmc('upper_ui', k))
            yld_str = 'Total YLD (preliminary):\n %.2f, (%.2f, %.2f)' % (est_yld, est_yld_lower_ui, est_yld_upper_ui)
            pl.text(.5 * (xmin + xmax), .15 * (ymin + ymax), yld_str)

        if rows == 2 and cols == 4:
            pl.subplot(rows, cols, rows*cols)
            emp_prior = dm.get_empirical_prior(type)
            if emp_prior.has_key('beta'):
                pl.title('emp prior coefficients', fontsize=8)
                y = np.array(emp_prior['beta'])
                pl.plot(y, '.-', alpha=.5, label=type, color=color_for.get(type, 'black'), linewidth=2)

def table_disease_model(dm_json, keys):
    """Make a tabular representation of the disease model data and
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
    type, region, year, sex = keys[0].split(dismod3.utils.KEY_DELIM_CHAR)
    cnt = len(keys)
        
    wb = Workbook()
    ws0 = wb.add_sheet('0')
    
    date = strftime("%Y/%m/%d")
    time = strftime("%H:%M:%S")
    user = "username"
    ws0.write(0, 0, "Dismod III output, date: %s, time: %s, user: %s" % (date, time, user))
    #ws0.write(2, 0, "Disease: %s" % (dm.get_condition())
    ws0.write(2, 0, "Condition: %s" % (dm.params['condition']))
    ws0.write(3, 0, "Region: %s" % (region))
    ws0.write(4, 1, "%s %s" % (sex.capitalize(), year))
    for i in range(1, 5):
        ws0.write(5, i, "Data")
    for i in range(5, 9):
        ws0.write(5, i, "Prior")
    for i in range(9, 31):
        ws0.write(5, i, "Posterior")
    ws0.write(6, 0, "Age")
    ws0.write(6, 1, "Prevalence")
    ws0.write(6, 2, "Incidence")
    ws0.write(6, 3, "Remission")
    ws0.write(6, 4, "Excess Mortality")
    ws0.write(6, 5, "Prevalence")
    ws0.write(6, 6, "Incidence")
    ws0.write(6, 7, "Remission")
    ws0.write(6, 8, "Excess Mortality")
    ws0.write(6, 9, "Prevalence")
    ws0.write(6, 10, "Prevalence")
    ws0.write(6, 11, "Prevalence")
    ws0.write(6, 12, "Incidence")
    ws0.write(6, 13, "Incidence")
    ws0.write(6, 14, "Incidence")
    ws0.write(6, 15, "Remission")
    ws0.write(6, 16, "Remission")
    ws0.write(6, 17, "Remission")
    ws0.write(6, 18, "Excess Mortality")
    ws0.write(6, 19, "Excess Mortality")
    ws0.write(6, 20, "Excess Mortality")
    ws0.write(6, 21, "Duration")
    ws0.write(6, 22, "Duration")
    ws0.write(6, 23, "Duration")
    ws0.write(6, 24, "With-condition")
    ws0.write(6, 25, "With-condition")
    ws0.write(6, 26, "With-condition")
    ws0.write(6, 27, "RR mortality")
    ws0.write(6, 28, "RR mortality")
    ws0.write(6, 29, "RR mortality")
    ws0.write(6, 30, "Age of onset")
    for i in range(1, 10):
        ws0.write(7, i, "(rate)")
    ws0.write(7, 10, "lower ui")
    ws0.write(7, 11, "upper ui")
    ws0.write(7, 12, "(rate)")
    ws0.write(7, 13, "lower ui")
    ws0.write(7, 14, "upper ui")
    ws0.write(7, 15, "(rate)")
    ws0.write(7, 16, "lower ui")
    ws0.write(7, 17, "upper ui")
    ws0.write(7, 18, "(rate)")
    ws0.write(7, 19, "lower ui")
    ws0.write(7, 20, "upper ui")
    ws0.write(7, 21, "(years)")
    ws0.write(7, 22, "lower ui")
    ws0.write(7, 23, "upper ui")
    ws0.write(7, 24, "(rate)")
    ws0.write(7, 25, "lower ui")
    ws0.write(7, 26, "upper ui")
    ws0.write(7, 27, "(rate)")
    ws0.write(7, 28, "lower ui")
    ws0.write(7, 29, "upper ui")
    ws0.write(7, 30, "(years)")
    for j in range(101):
        ws0.write(8 + j, 0, j)
        ws0.write(8 + j, 30, j + .5)

    for ii, k in enumerate(keys):
        type, region, year, sex = k.split(dismod3.utils.KEY_DELIM_CHAR)
        data_type = clean(type) + ' data'
        data = data_hash.get(data_type, region, year, sex) \
               + data_hash.get(data_type, region, year, 'total')
            # ws0.write(8 + jj, 0, "%s-%s" % (start, end))
        column = 0
        if type == 'prevalence':
            column = 1
        if type == 'incidence':
            column = 2
        if type == 'remission':
            column = 3
        if type == 'case-fatality':
            column = 4
        for jj in range(len(data)):
            start = data[jj]['age_start']
            end = data[jj]['age_end'] + 1
            for j in range (start, end):
                ws0.write(8 + j, column, data[jj]['value'])
        if type == 'prevalence':
            column = 5
        if type == 'incidence':
            column = 6
        if type == 'remission':
            column = 7
        if type == 'case-fatality':
            column = 8
        if len(dm.get_mcmc('emp_prior_mean', k)) == 101:
            for j in range(0, 101):
                ws0.write(8 + j, column, dm.get_mcmc('emp_prior_mean', k)[j])
        if type == 'prevalence':
            column = 9
        if type == 'incidence':
            column = 12
        if type == 'remission':
            column = 15
        if type == 'case-fatality':
            column = 18
        if type == 'duration':
            column = 21
        if type == 'relative-risk':
            column = 27

        for j in range(0, 101):
            ws0.write(8 + j, column, dm.get_mcmc('mean', k)[j])
        if type == 'prevalence':
            column = 10
        if type == 'incidence':
            column = 13
        if type == 'remission':
            column = 16
        if type == 'case-fatality':
            column = 19
        if type == 'duration':
            column = 22
        if type == 'relative-risk':
            column = 28
        for j in range(0, 101):
            ws0.write(8 + j, column, dm.get_mcmc('lower_ui', k)[j])
        if type == 'prevalence':
            column = 11
        if type == 'incidence':
            column = 14
        if type == 'remission':
            column = 17
        if type == 'case-fatality':
            column = 20
        if type == 'duration':
            column = 23
        if type == 'relative-risk':
            column = 29
        for j in range(0, 101):
            ws0.write(8 + j, column, dm.get_mcmc('upper_ui', k)[j])
    wb.save('output.xls')

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
    ymax = dm.get_ymax()
    
    sorted_regions = sorted(dismod3.gbd_regions, reverse=False,
                            key=lambda r: len(data_hash.get(region=r)))
    for ii, region in enumerate(sorted_regions):
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
                #plot_empirical_prior(dm, dismod3.gbd_key_for(type, region, year, sex),
                #                     color=color_for.get(type, 'black'))
                type = ' '.join([type, 'data'])
                data = data_hash.get(type, region, year, sex) + data_hash.get(type, region, year, 'total')
                if len(data) > max_intervals:
                    data = random.sample(data, max_intervals)
                plot_intervals(dm, data, color=color_for.get(type, 'black'), linewidth=1, alpha=.25)
            pl.xticks([])
            pl.yticks([])
            pl.axis([xmin, xmax, ymin, ymax])

    ii += 1
    subplot_px[dismod3.gbd_key_for('all', 'world', 1997, 'total')] = \
                                          ', '.join(['0',
                                                     str(int(100 * (rows - ii - 1) * subplot_height)),
                                                     str(int(100 * (jj + 1) * subplot_width)),
                                                     str(int(100 * (rows - ii) * subplot_height))])

    if boxes_only:
        return subplot_px
    
    fig.add_axes([0,
                  ii*subplot_height / fig_height,
                  fig_width,
                  subplot_height / fig_height],
                 frameon=False)
    for type in dismod3.data_types:
        type = type.replace(' data', '')
        plot_map_fit(dm, dismod3.gbd_key_for(type, 'world', 1997, 'total'),
                     linestyle='-', color=color_for.get(type, 'black'))
        pl.xticks([])
        pl.yticks([])
        pl.axis([xmin, xmax, ymin, ymax])

def plot_prior_preview(dm):
    """ Generate a preview of what a rate function with this prior looks like"""

    fig = pl.figure(figsize=(3.4, 3.4), dpi=100)
    pl.clf()

    ages = dm.get_estimate_age_mesh()
    xmin = ages[0]
    xmax = ages[-1]
    ymin = 0.
    ymax = dm.get_ymax()

    for ii, type in enumerate(['prevalence', 'incidence', 'remission', 'case-fatality']):
        pl.subplot(2,2,ii+1)
        prior_str = dm.get_global_priors(type)

        vars = dismod3.utils.prior_vals(dm, type)

        ages = dm.get_estimate_age_mesh()
        color = color_for.get(type, 'black')
        mu = vars['rate_stoch'].value
        dispersion = vars['dispersion'].value

        pl.plot(ages, mu, color=color, linestyle='-', linewidth=2)
        lb = mc.invlogit(mc.logit(mu) - 1.96*dispersion)
        ub = mc.invlogit(mc.logit(mu) + 1.96*dispersion)
        plot_uncertainty(ages, lb, ub, edgecolor=color, alpha=.75)

        pl.text(.9 * xmin + .1 * xmax, .9 * ymax + .1 * ymin, type, color=color)
        plot_prior(dm, type)
        plot_intervals(dm, vars['data'], color=color)
        #pl.semilogy([xmin], [ymax])
        
        pl.xticks([])
        pl.yticks(fontsize=8)
        pl.axis([xmin, xmax, ymin, ymax])

def plot_intervals(dm, data, **params):
    """
    use matplotlib plotting functions to render transparent
    rectangles on the current figure representing each
    piece of Data
    """
    default_params = dict(alpha=.35, color=(.0,.5,.0), linewidth=5)
    default_params.update(**params)
    
    errorbar_params = copy.copy(default_params)
    errorbar_params['linewidth'] = 1

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
                    **errorbar_params)

        pl.plot(np.array([d['age_start'], d['age_end']+1.]),
                np.array([val, val]),
                **default_params)

        if clean(d.get('self_reported', '')) == 'true':
            pl.text(.5*(d['age_start']+d['age_end']), val, 'self-reported', fontsize=6, horizontalalignment='center', verticalalignment='center')

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
    
    lb = dm.get_mcmc('lower_ui', type)
    ub = dm.get_mcmc('upper_ui', type)

    if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
        #lb = lb[param_mesh]
        #ub = ub[param_mesh]
        #x = np.concatenate((param_mesh, param_mesh[::-1]))
        plot_uncertainty(age, lb, ub, edgecolor=color, alpha=1., zorder=2.)

    val = dm.get_mcmc('mendian', type)

    if len(age) > 0 and len(age) == len(val):
        pl.plot(age, val, color=color, linewidth=2, alpha=.75)

    c = dm.get_mcmc('dispersion', type)
    if len(c) == 5:
        #pl.text(age[3*len(age)/5], 0, 'dispersion:\n%.3f (%.3f,%.3f)' % (c[2], c[0], c[4]), fontsize=8)

        # plot dispersion + uncertainty for rough est of 95% uncertainty interval for data
        lb = mc.invlogit(mc.logit(lb) - 1.96*c[2])
        ub = mc.invlogit(mc.logit(ub) + 1.96*c[2])

        if show_data_ui and len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
            plot_uncertainty(age, lb, ub, linestyle='dashed', edgecolor=color, facecolor=(.95,.95,.95), label='Data 95% UI', alpha=.95, zorder=1.)

def plot_empirical_prior(dm, type, color=(.2,.2,.2)):
    age = dm.get_estimate_age_mesh()
    lb = dm.get_mcmc('emp_prior_lower_ui', type)
    ub = dm.get_mcmc('emp_prior_upper_ui', type)
    if len(age) > 0 and len(age) == len(lb) and len(age) == len(ub):
        plot_uncertainty(age, lb, ub, linestyle='dotted', fill=True, edgecolor=color, alpha=.5, zorder=0.)
        #pl.plot(age, lb, linestyle='dotted', color=color, linewidth=1, alpha=.5, zorder=1.5)
        #pl.plot(age, ub, linestyle='dotted', color=color, linewidth=1, alpha=.5, zorder=1.5)

    val = dm.get_mcmc('emp_prior_mean', type)
    if len(age) > 0 and len(age) == len(val):
        pl.plot(age, val, color=color, linewidth=1, alpha=.5, linestyle='dotted', zorder=0.)

def plot_uncertainty(ages, lower_bound, upper_bound, **params):
    default_params = {'facecolor': '.8'}
    default_params.update(**params)

    x = np.concatenate((ages, ages[::-1]))
    y = np.concatenate((lower_bound, upper_bound[::-1]))
    pl.fill(x, y, **default_params)
                
def plot_prior(dm, type):
    # show 'zero' priors
    for prior_str in dm.get_priors(type).split(dismod3.PRIOR_SEP_STR):
        prior = prior_str.split()
        if len(prior) > 0 and prior[0] == 'zero':
            age_start = int(prior[1])
            age_end = int(prior[2]) + .5

            pl.plot([age_start, age_end], [0, 0], color='red', linewidth=15, alpha=.75)

    # write out details of priors in a friendly font as well
    if len(dm.get_estimate_age_mesh()) > 0:
        a0 = dm.get_estimate_age_mesh()[0]
        v0 = 0.
        pl.text(a0, v0, ' Priors:\n' + dm.get_priors(type).replace(dismod3.PRIOR_SEP_STR, '\n'), color='black', family='monospace', fontsize=8, alpha=.75)

    # write coeff vals for empirical priors as well, if available
#     emp_prior = dm.get_empirical_prior(type.split('+')[0])
#     alpha = emp_prior.get('alpha')
#     if alpha != None:
#         import logit_normal_model as rate_model
#         Xa, Xb = rate_model.regional_covariates(type)
        
#         coeffs = ['%.3f' % aa for ii, aa in enumerate(alpha) if Xa[ii] != 0.]
#         l,r,b,t = pl.axis()
            
#         pl.text(30, 0., ' '.join(coeffs), fontsize=10, family='monospace', alpha=.8, color='black')
        
def clear_plot(width=4*1.5, height=3*1.5):
    fig = pl.figure(figsize=(width,height))
    pl.clf()
    return fig

def label_plot(dm, type, **params):
    pl.xlabel('Age (years)', **params)
    pl.ylabel('%s %s' % (type, dm.get_units(type) or ''), **params)
    pl.title('%d: %s; %s; %s; %s' % \
                 (dm.params['id'], prettify(dm.params['condition']),
                  dm.params['sex'], prettify(dm.params['region']),
                  dm.params['year']), **params)
    #pl.legend()

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
          'case-fatality data', 'all-cause mortality data', 'duration data' or 'all'
        region : str, one of the 21 gbd regions or 'World' or 'all'
        year : int, one of 1990, 2005, 'all'
        sex : str, one of 'male', 'female', 'total', 'all'
        """
        if not self.d_hash.has_key((type, region, year, sex)):
            self.d_hash[(type, region, year, sex)] = [d for d in self.data if dismod3.relevant_to(d, type, region, year, sex)]
        return self.d_hash[(type, region, year, sex)]
