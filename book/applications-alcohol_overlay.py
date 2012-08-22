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

pred = pandas.read_csv('H:\\gbd\\book\\applications-alcohol_pred.csv', index_col=0)
#ui_pred = pandas.read_csv('H:\\gbd\\book\\applications-alcohol_ui_pred.csv', index_col=0)   

pl.figure(**book_graphics.full_page_params)

param_list = [dict(type='p', title='(a)', ylabel='Prevalence (%)', yticks=([0, .05, .1, .15, .2], [0, 5, 10, 15, 20]), axis=[-5,105,-0.024,.24]),
          dict(type='pf', title='(b)', ylabel='Cause-specific \n mortality \n (per 100,000 PY)', yticks=([0, .00005, .0001, .00015, .0002], [0, 5, 10, 15, 20]), axis=[-5,105,-.000024,.00024]),
          dict(type='f', title='(c)', ylabel='Excess mortality \n (per 100 PY) \n', yticks=([0, .01, .02, .03, .04], [0, 1, 2, 3, 4]), axis=[-5,105,-.005,.05]),
        ]

for i, params in enumerate(param_list):
    ax = pl.subplot(2,2,i+1)

    pl.plot(pl.arange(101), pred['csmr_'+params['type']], 'k-', linewidth=3, label='Posterior mean, $h_f'' \geq 0$')
    pl.plot(pl.arange(101), pred['pf_'+params['type']], 'k--', linewidth=3, label='Posterior mean, $h_f'' = 0$')

    pl.xlabel('Age (Years)')
    if params['type']=='pf': pl.ylabel(params['ylabel']+'\n\n\n', ha='center')
    else: pl.ylabel(params['ylabel']+'\n\n', ha='center')
    pl.axis(params.get('axis', [-5,105,-.005,.06]))
    subtitle(params['title'])
    pl.grid()
    pl.yticks(*params.get('yticks', ([0, .025, .05], [0, 2.5, 5])))
    pl.legend(loc='upper right', bbox_to_anchor=(2.1,1), fancybox=True, shadow=True) 
    
pl.subplots_adjust(hspace=.35)
pl.subplots_adjust(wspace=.45)

pl.savefig('H:\\gbd\\book\\applications\\alcohol-overlay.pdf')
pl.savefig('H:\\gbd\\book\\applications\\alcohol-overlay.png')