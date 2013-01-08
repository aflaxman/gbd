import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas
import matplotlib.mpl as mpl

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)

import data_model
reload(data_model)

# make all fonts bigger, etc
mpl.rcParams['axes.titlesize'] = 'xx-large'
mpl.rcParams['axes.labelsize'] = 'xx-large'
mpl.rcParams['xtick.labelsize'] = 'x-large'
mpl.rcParams['ytick.labelsize'] = 'x-large'
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['text.fontsize'] = 12

# axis def
def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])
    
# subtitle func
def subtitle(s):
    """ title where the panel names appear within each panel"""
    l,r,b,t=pl.axis()
    x = l + (r-l)*.05
    y = t - (t-b)*.05
    pl.text(x, y, s, ha='left', va='top')

def my_plot_data_bars(df, color, label, style='book'):
    """ Plot some data bars
    Input
    -----
    df : pandas.DataFrame with columns age_start, age_end, value
    """
    data_bars = zip(df['age_start'], df['age_end'], df['value'])

    # show at most 500 bars, to keep things fast
    # TODO: make 500 into an option
    if len(data_bars) > 500:
        import random
        data_bars = random.sample(data_bars, 500)

    # make lists of x and y points, faster than ploting each bar
    # individually
    x = []
    y = []
    for a_0i, a_1i, p_i in data_bars:
        x += [a_0i, a_1i, pl.nan]
        y += [p_i, p_i, pl.nan]

    pl.plot(x, y, 's-', mew=1, mec='w', ms=4, color=color, label=label)    
    
# load data to plot
age_pred = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_fruit_age_pred.csv', index_col=0)
ui_pred = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_fruit_ui_pred.csv', index_col=0)
isl = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_fruit_isl.csv', index_col=0)
grc = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_fruit_grc.csv', index_col=0)
we = pandas.read_csv('/homes/peterhm/gbd/book/applications-data_fruit_we.csv', index_col=0)

# figure with 2 subplots and legends outside of plots
labeling = dict(we_model=dict(style='k-', lab='Negative-binomial '), 
                we_log_model=dict(style='k--', lab='Lognormal '), 
                we_norm_model=dict(style='k:', lab='Normal '))

pl.figure(**book_graphics.full_page_params)
knots = list(ui_pred.index)
k=0
for i in ['ISL', 'GRC']:
    k=k+1
    pl.subplot(1,2,k)
    my_plot_data_bars(we, color='grey', label='Western Europe') 
    if i == 'ISL': my_plot_data_bars(isl, color='black', label='Iceland')
    elif i == 'GRC': my_plot_data_bars(grc, color='black', label='Greece')

    for m in ['we_model', 'we_log_model', 'we_norm_model']:
        pred = pl.array(age_pred[m + '_' + i])
        we_hpd_l = pl.array(ui_pred[m + '_' + i + '_l'])
        we_hpd_u = pl.array(ui_pred[m + '_' + i + '_u'])

        pl.plot(pred, labeling[m]['style'], linewidth=3, label=labeling[m]['lab']+'posterior mean')
        pl.plot(knots, we_hpd_l, labeling[m]['style'], linewidth=1, label=labeling[m]['lab']+'95% HPD interval')
        pl.plot(knots, we_hpd_u, labeling[m]['style'], linewidth=1)

    pl.xlabel('Age (years)')
    pl.ylabel('Consumption (kg/d)')
    pl.yticks([0, .015, .03, .045, .06], [0, 0.15, 0.30, 0.45, 0.6])
    my_axis(.075)
    if i == 'ISL': subtitle('(a)')
    elif i == 'GRC': subtitle('(b)')
    pl.legend(loc='upper center', bbox_to_anchor=(.5,-.33), fancybox=True, shadow=True)
    pl.grid()
    
pl.subplots_adjust(top=.93, bottom=.53, wspace=.35)

pl.savefig('/homes/peterhm/gbd/book/applications/fruit-isl_rate_type.pdf')
pl.savefig('/homes/peterhm/gbd/book/applications/fruit-isl_rate_type.png')