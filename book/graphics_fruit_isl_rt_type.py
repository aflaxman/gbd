import sys
sys.path += ['../gbd', '../gbd/book', '../dm3-computation_only/', '../dm3-computation_only/book']
import pylab as pl
import pymc as mc
import pandas

import dismod3
reload(dismod3)

import book_graphics
reload(book_graphics)

import data_model
reload(data_model)

# set font
book_graphics.set_font()

# axis def
def my_axis(ymax):
    pl.axis([-5,105,-ymax/10.,ymax])
    
# load data to plot
age_pred = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_age_pred.csv', index_col=0)
ui_pred = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_ui_pred.csv', index_col=0)
isl = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_isl.csv', index_col=0)
grc = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_grc.csv', index_col=0)
we = pandas.read_csv('/home/j/Project/dismod/gbd/data/applications-data_fruit_we.csv', index_col=0)

# figure with 2 subplots and legends outside of plots
labeling = dict(we_model=dict(style='k-', lab='Negative-binomial '), 
                we_log_model=dict(style='k--', lab='Lognormal '), 
                we_norm_model=dict(style='k:', lab='Normal '))

fig = pl.figure(**book_graphics.three_quarter_page_params)
knots = list(ui_pred.index)

# plot 1
ax1 = fig.add_subplot(1,2,1)
ax1.plot(dismod3.graphics.plot_data_bars(we, color='grey', label='Western Europe') )
ax1.plot(dismod3.graphics.plot_data_bars(isl, color='black', label='Iceland'))
ax1.plot(dismod3.graphics.plot_data_bars(grc, color='black', label='Greece'))
pl.xlabel('Age (years)')
pl.ylabel('Consumption (kg/d)')
pl.yticks([0, .015, .03, .045, .06], [0, 0.15, 0.30, 0.45, 0.6])
my_axis(.075)
# plot 2, 3


    
for i in ['GRC', 'ISL']: #in ['we_model', 'we_log_model', 'we_norm_model']:
    pred = pl.array(age_pred['we_model_' + i])
    we_hpd_l = pl.array(ui_pred['we_model_' + i + '_l'])
    we_hpd_u = pl.array(ui_pred['we_model_' + i + '_u'])
    # plot 1
    ax1.plot(pred, labeling['we_model']['style'], linewidth=3, label=labeling['we_model']['lab']+'posterior mean')
    ax1.plot(knots, we_hpd_l, labeling['we_model']['style'], linewidth=1, label=labeling['we_model']['lab']+'95% HPD interval')
    ax1.plot(knots, we_hpd_u, labeling['we_model']['style'], linewidth=1)  
    for m in ['we_log_model', 'we_norm_model']:
        m_pred = pl.array(age_pred[m + '_' + i])
        if i == 'GRC': 
            ax2 = fig.add_subplot(2,2,2)
            ax2.plot((pred - m_pred), labeling[m]['style'], linewidth=3, label='Negative-binomial - ' + labeling[m]['lab'])
            pl.axis([-5, 105, -.007, .017])
            pl.yticks([-.005, -.0, .005, .01, .015])
            pl.ylabel('Greece')
        elif i == 'ISL': 
            ax3 = fig.add_subplot(2,2,4, sharex=ax2)
            ax3.plot((pred - m_pred), labeling[m]['style'], linewidth=3, label='Negative-binomial - ' + labeling[m]['lab'])
            pl.axis([-5, 105, -.0022, .0008])
            pl.yticks([-.002, -.0015, -.001, -.0005, .0, .0005])
            pl.xlabel('Age (years)')
            pl.ylabel('Iceland')
        
pl.setp(ax2.get_xticklabels(), visible=False)
pl.legend(loc='upper center', bbox_to_anchor=(-.2,-.33), fancybox=True, shadow=True, ncol=3)
    
pl.subplots_adjust(top=.99, bottom=.27, wspace=.45, hspace=.1)
pl.savefig('book/graphics/fruit-isl_rate_type.pdf')
pl.savefig('book/graphics/fruit-isl_rate_type.png')

pl.show()